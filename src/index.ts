import { McpAgent } from "agents/mcp";
import { McpServer } from "@modelcontextprotocol/sdk/server/mcp.js";
import { z } from "zod";

import { findBestMatch } from 'string-similarity'

import rawDrugMap  from "./data/CIViC_therapy_name_map.json";
import rawGeneMap  from "./data/VICC_gene_alias_map.json";

export const dDrugMap: Record<string, string[]>  = rawDrugMap  as Record<string, string[]>;
export const dGeneMap: Record<string, string[]>  = rawGeneMap  as Record<string, string[]>;

function normalizeStr(s: string): string {
  return s
    .normalize('NFKD')                   // decompose accents
    .replace(/[\u0300-\u036f]/g, '')     // strip diacritics
    .toLowerCase()
    .replace(/[^a-z0-9\s]/g, ' ')        // keep only letters, digits, spaces
    .replace(/\s+/g, ' ')                // collapse runs of spaces
    .trim()
}

type AliasIndex = {
  aliasToPrimary: Record<string, string>;
  candidates: string[];
};

function buildAliasIndex(lookup: Record<string, string[]>): AliasIndex {
  const aliasToPrimary: Record<string, string> = {};

  for (const [primary, aliases] of Object.entries(lookup)) {
    aliasToPrimary[normalizeStr(primary)] = primary;
    for (const alias of aliases) {
      aliasToPrimary[normalizeStr(alias)] = primary;
    }
  }

  return { aliasToPrimary, candidates: Object.keys(aliasToPrimary) };
}

let GENE_INDEX: AliasIndex | undefined;
let DRUG_INDEX: AliasIndex | undefined;

function getGeneIndex() {
  return (GENE_INDEX ??= buildAliasIndex(dGeneMap));
}
function getDrugIndex() {
  return (DRUG_INDEX ??= buildAliasIndex(dDrugMap));
}

export function normalizeEntityFast(
  name: string | undefined | null,
  index: AliasIndex,
  threshold = 0.7
): string | undefined {
  if (!name) return undefined;

  const qNorm = normalizeStr(name);

  // Check for exact match first
  if (qNorm in index.aliasToPrimary) {
    return index.aliasToPrimary[qNorm];
  }

  // Fall back to fuzzy matching
  const { bestMatch } = findBestMatch(qNorm, index.candidates);

  return bestMatch.rating >= threshold
    ? index.aliasToPrimary[bestMatch.target]
    : undefined;
}

// export function normalizeEntityFast(
//   name: string | undefined | null,
//   index: AliasIndex,
//   threshold = 0.7
// ): string | undefined {
//   if (!name) return undefined;

//   const qNorm = normalizeStr(name);
//   const { bestMatch } = findBestMatch(qNorm, index.candidates);

//   return bestMatch.rating >= threshold
//     ? index.aliasToPrimary[bestMatch.target]
//     : undefined;
// }

// /**
//  * Map a free‑form name to its primary alias via fuzzy matching.
//  *
//  * @param name      Input string (may be undefined or null)
//  * @param lookup    Record<primary, aliases[]>
//  * @param threshold Minimum similarity (0–1) to accept a match
//  * @returns         The matched primary string, or undefined if below threshold or name missing
//  */
// export function normalizeEntity(
//   name: string | undefined | null,
//   lookup: Record<string, string[]>,
//   threshold: number = 0.7
// ): string | undefined {
//   if (!name) {
//     return undefined
//   }

//   const qNorm = normalizeStr(name)

//   // Build a map from normalized‑alias → primary
//   const aliasToPrimary: Record<string, string> = {}
//   for (const [primary, aliases] of Object.entries(lookup)) {
//     aliasToPrimary[normalizeStr(primary)] = primary
//     for (const alias of aliases) {
//       aliasToPrimary[normalizeStr(alias)] = primary
//     }
//   }

//   const candidates = Object.keys(aliasToPrimary)
//   const { bestMatch } = findBestMatch(qNorm, candidates)

//   return bestMatch.rating >= threshold
//     ? aliasToPrimary[bestMatch.target]
//     : undefined
// }

type DrugInfoNode = {
  name: string;
  conceptId?: string;
  approved?: boolean | null;
  immunotherapy?: boolean | null;
  antiNeoplastic?: boolean | null;
  drugAttributes?: { name: string; value: string }[];
  drugApprovalRatings?: {
    rating?: number | null;
    source?: {
      sourceDbName?: string;
      sourceTrustLevel?: { level?: string };
    };
  }[];
};

export type InteractionType = { type?: string; directionality?: string };
export type Interaction = {
  interactionScore?: number | null;
  interactionTypes?: InteractionType[];

  gene?: { name?: string | null } | null;
  drug?: { name?: string | null; approved?: boolean | null } | null;

  /** Only PubMed URLs — no raw pmids */
  pmidUrls?: string[];

  /** List of source database names */
  sourceDbNames?: string[];
};
export type DrugNode = { name: string; interactions: Interaction[] };

export interface GeneNode {
  /** Gene (or molecular profile) name. */
  name: string;
  /** Array of interaction objects. */
  interactions: Interaction[];
  /** Any other fields returned by the API. */
  [key: string]: unknown;
}

type NodeWithInteractions = { name: string; interactions?: Interaction[] };


/**
 * Return up to N interactions with:
 *  1) approved drugs first (drug.approved === true)
 *  2) then by highest interactionScore (desc)
 *
 * @param interactions - Array of interaction objects.
 * @param N            - Maximum number to return (default 20).
 */
export function filterInteractionScore(
  interactions: Interaction[],
  N = 20
): Interaction[] {
  return [...interactions]
    .sort((a, b) => {
      const aApproved = !!a.drug?.approved;
      const bApproved = !!b.drug?.approved;

      // Primary: approved first
      if (aApproved !== bApproved) {
        return aApproved ? -1 : 1;
      }

      // Secondary: higher interactionScore first
      const aScore = a.interactionScore ?? 0;
      const bScore = b.interactionScore ?? 0;
      return bScore - aScore;
    })
    .slice(0, N);
}

/** Exact (case-insensitive) match wins; otherwise pick the node that contains the term and has the most interactions. */
export function findBestNode<T extends NodeWithInteractions>(
  nodes: T[],
  term: string
): T | undefined {
  const t = term.toLowerCase();

  // 1) exact match first
  const exact = nodes.find(n => n.name.toLowerCase() === t);
  if (exact) return exact;

  // 2) most interactions among substring matches
  let best: T | undefined;
  let bestLen = -1;

  for (const n of nodes) {
    if (n.name.toLowerCase().includes(t)) {
      const L = n.interactions?.length ?? 0;
      if (L > bestLen) {
        best = n;
        bestLen = L;
      }
    }
  }

  return best;
}

/**
 * Fair-share allocation of a total budget across names based on available counts.
 * Ascending pass; each gets at most avg remaining (floor), capped by its available count.
 */
export function distributeNodes(
  dNodeLens: Record<string, number>,
  total = 100
): Record<string, number> {
  const entries = Object.entries(dNodeLens).sort((a, b) => a[1] - b[1]); // ascending
  const out: Record<string, number> = {};
  let i = 0;
  const num = entries.length;

  for (const [key, val] of entries) {
    const allowed = Math.floor(total / (num - i)) || 0;
    const curr = Math.min(val, Math.max(allowed, 0));
    out[key] = curr;
    total -= curr;
    i += 1;
  }
  return out;
}

/** Approved drugs first, then highest interactionScore (desc). Missing → approved=false, score=0. */
export function filterAndSort(interactions: Interaction[], N = 20): Interaction[] {
  return [...(interactions ?? [])]
    .sort((a, b) => {
      const aApproved = !!a.drug?.approved;
      const bApproved = !!b.drug?.approved;
      if (aApproved !== bApproved) return aApproved ? -1 : 1; // approved first

      const aScore = a.interactionScore ?? 0;
      const bScore = b.interactionScore ?? 0;
      return bScore - aScore; // desc
    })
    .slice(0, N);
}

/**
 * For each requested name:
 *  - choose best node (exact > longest substring)
 *  - budget interactions across names (distributeNodes)
 *  - apply filterAndSort per name with that budget
 * Returns map: requested name → filtered interactions[]
 */
export function selectNodes<T extends NodeWithInteractions>(
  nodes: T[],
  normalizedNames: string[],
  totalBudget: number
): Record<string, Interaction[]> {
  const counts: Record<string, number> = {};
  const chosen: Record<string, T> = {};

  for (const name of normalizedNames) {
    const node = findBestNode(nodes, name);
    if (!node) continue;
    chosen[name] = node;
    counts[name] = node.interactions?.length ?? 0;
  }

  const allocations = distributeNodes(counts, totalBudget);

  const out: Record<string, Interaction[]> = {};
  for (const name of Object.keys(allocations)) {
    const node = chosen[name];
    const N = allocations[name] ?? 0;
    out[name] = filterAndSort(node?.interactions ?? [], N);
  }

  return out;
}

function convertNodesPMID(
  nodes: Record<string, Interaction[]>
): Record<string, Interaction[]> {
  const base = "https://pubmed.ncbi.nlm.nih.gov/";

  const out: Record<string, Interaction[]> = {};

  for (const [key, interactions] of Object.entries(nodes)) {
    out[key] = interactions.map(interaction => {
      // treat publications/sources as untyped
      const raw = interaction as any;

      const pmidUrls =
        (raw.publications ?? [])
          .filter((p: any) => p?.pmid)
          .map((p: any) => `${base}${p.pmid}/`);

      // remove publications
      const { publications, ...rest } = raw;

      // return only allowed Interaction fields plus new URL fields
      return {
        ...rest,
        pmidUrls,
      } as Interaction & {
        pmidUrls: string[];
      };
    });
  }

  return out;
}


// ========================================
// API CONFIGURATION - Customize for your GraphQL API
// ========================================

export const API_CONFIG = {
  name:        "DGIdb_MCP",
  version:     "0.1.0",
  description: "Fixed-schema MCP tools for the DGIdb GraphQL API",
  mcpSpecVersion: "2025-06-18",
  features: {
    structuredToolOutput: true,
    metaFields:           true,
    protocolVersionHeaders: true,
    titleFields:          true,
    toolAnnotations:      true,
  },
  headers: {
    "User-Agent": "MCPDGIdbServer/0.1.0"
  },
} as const;


/** -----------------------------------------------------------------
 *  Tool definitions
 *  ---------------------------------------------------------------- */

export const tools = {
  getDrugInfo: {
    name: "get_drug_info",
    description:
      "Gets drug info including approval, if used in immunotherapy, and other drug attributes for a list of drugs.",
    inputSchema: {drugNames: z.array(z.string()).min(1).max(25).describe("Drug names only")},
    annotations: {
      destructive: false,
      idempotent:  true,
      cacheable:   false,
      world_interaction: "open",
      side_effects: ["external_api_calls"],
      resource_usage: "network_io_heavy",
    },
    async handler (
  		{ drugNames }: { drugNames: string[] } ) {

    const query = /* GraphQL */ `
          query drugs_info($names: [String!]) {
            drugs(names: $names) {
                nodes {
                    name
                    conceptId
                    approved
                    immunotherapy
                    antiNeoplastic
                    drugAttributes {
                        name
                        value
                        }
                    drugApprovalRatings {
                        rating
                        source {
                            sourceDbName
                            sourceTrustLevel {
                                level
                            }
                        }
                    }
                }
            }
        }	
      `;
    
      const rawNames = (drugNames ?? [])
        .map(s => (s ?? "").trim())
        .filter(Boolean);


      const normalizedNames = rawNames.map(raw =>
        normalizeEntityFast(raw, getDrugIndex(), 0.7) || raw
      );

      const res = await fetch("https://dgidb.org/api/graphql", {
        method: "POST",
        headers: { "Content-Type": "application/json", ...API_CONFIG.headers },
        body: JSON.stringify({ query, variables: { names: normalizedNames } })
      }).then(r => r.json()) as {
        data?: { drugs?: { nodes?: DrugInfoNode[] } };
        errors?: unknown[];
      };

      if (!res.data?.drugs?.nodes) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: JSON.stringify(res, null, 2) }]
        };
      }

      const nodes = res.data.drugs.nodes;
      if (!nodes.length) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: `No drug nodes found for: ${normalizedNames.join(", ")}` }]
        };
      }


      return {
        isError: false,
        content: [{ type: "text" as const, text: JSON.stringify(nodes, null, 2) }]
      };
    
    
    },
  },

  getGeneInfo: {
    name: "get_gene_info",
    description:
      "Gets gene category info for a list of genes.",
    inputSchema: {geneNames: z.array(z.string()).min(1).max(25).describe("HGNC gene symbols only")},
    annotations: {
      destructive: false,
      idempotent:  true,
      cacheable:   false,
      world_interaction: "open",
      side_effects: ["external_api_calls"],
      resource_usage: "network_io_heavy",
    },

    async handler(
        { geneNames }: { geneNames: string[] }  
    ) {

    const query = /* GraphQL */ `
    query gene_category($names: [String!]){
        genes(names: $names) {
            nodes {
                name
                longName
                conceptId
                geneCategoriesWithSources {
                  name
                  sourceNames
                }
            }
        }
    }	
    `;

    const rawNames = (geneNames ?? [])
      .map(s => (s ?? "").trim())
      .filter(Boolean);

    //const rawNames = geneNames.split(",").map(s => s.trim()).filter(Boolean);
    const normalizedNames = rawNames.map(raw =>
      normalizeEntityFast(raw, getGeneIndex(), 0.7) || raw
    );


    const upperNames = normalizedNames.map(n => n.toUpperCase());

    const res = await fetch("https://dgidb.org/api/graphql", {
				method: "POST",
				headers: { "Content-Type": "application/json", ...API_CONFIG.headers },
				body: JSON.stringify({ query, variables: { names: upperNames } })
			}).then(r => r.json()) as {
        data?: { genes?: { nodes?: Array<{
          name: string;
          longName?: string;
          conceptId?: string;
          geneCategoriesWithSources?: Array<{ name: string; sourceNames: string[] }>;
        }> } };
        errors?: unknown[];
      };

    if (!res.data?.genes?.nodes) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: JSON.stringify(res, null, 2) }]
        };
      }

      const nodes = res.data.genes.nodes;
      if (!nodes.length) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: `No gene nodes found for: ${upperNames.join(", ")}. DEBUG normalizedNames: ${JSON.stringify(normalizedNames)}` }]
        };
      }

      return {
        isError: false,
        content: [{ type: "text" as const, text: JSON.stringify(nodes, null, 2) }]
      };
    },
  },


  getGeneInteractionsForDrug: {
    name: "get_gene_interactions_for_drug_list",
    description:
      "Use this when the user provides drug(s) and asks for genes that interact. Only the top ranked interactions will be returned.",
    inputSchema: {drugNames: z.array(z.string()).min(1).max(25).describe("Drug names only")},
    annotations: {
      destructive: false,
      idempotent:  true,
      cacheable:   false,
      world_interaction: "open",
      side_effects: ["external_api_calls"],
      resource_usage: "network_io_heavy",
    },

    async handler(
  		{ drugNames }: { drugNames: string[] }  // <- matches schema
	) {

		const query = /* GraphQL */ `
			query drugs($names: [String!]) {
        drugs(names: $names) {
          nodes {
              name
              interactions {
                  gene {
                      name
                  }
                  interactionScore
                  interactionTypes {
                      type
                      directionality
                  }
                  publications {
                    pmid
                  }
                  sources {
                    sourceDbName
                    baseUrl
                  }
              }
          }
        }
      }`;

    const rawNames = (drugNames ?? [])
      .map(s => (s ?? "").trim())
      .filter(Boolean);


    // const rawNames = drugNames.split(",").map(s => s.trim()).filter(Boolean);
    // Call normalizeEntity(rawName, dDrugMap, 0.7). Fallback to raw if it returns false.
    const normalizedNames = rawNames.map(raw =>
      normalizeEntityFast(raw, getDrugIndex(), 0.7) || raw
    );
    const totalBudget = normalizedNames.length === 1 ? 15 : 100;

    const res = await fetch("https://dgidb.org/api/graphql", {
      method: "POST",
      headers: { "Content-Type": "application/json", ...API_CONFIG.headers },
      body: JSON.stringify({ query, variables: { names: normalizedNames } })
    }).then(r => r.json()) as {
      data?: { drugs?: { nodes?: DrugNode[] } };
      errors?: unknown[];
    };

    if (!res.data?.drugs?.nodes) {
      return {
        isError: true,
        content: [{ type: "text" as const, text: JSON.stringify(res, null, 2) }]
      };
    }

    const nodes = res.data.drugs.nodes;
    if (!nodes.length) {
      return {
        isError: true,
        content: [{ type: "text" as const, text: `No drug nodes found for: ${normalizedNames.join(", ")}` }]
      };
    }

    const instructions =
          "Interaction Score: Scoring metric based on the evidence supporting an interaction.\n" +
          "Interaction Type:  Nature of the association between a particular drug and gene.\n" +
          "Interaction Direction: Describes the effect of the interaction on the biological activity of the gene.\n" +
          "PMID url: Direct link to the PubMed publication supporting this interaction.\n" +
          "The returned list of interactions is limited to the top 15 and does not include all of them.\n" +
          "When returning information to users you MUST cite URLs used for specific information. Always cite the PubMed IDs and Database Sources as links";

    const selected_nodes = selectNodes(nodes, normalizedNames, totalBudget)

    const payload = {
        instructions,
        "API Results": convertNodesPMID(selected_nodes),
      };

    return {
      isError: false,
      content: [{ type: "text" as const, text: JSON.stringify(payload, null, 2) }]
    };
	},
}, 

  getDrugInteractionsForGene: {
    name: "get_drug_interactions_for_gene_list",
    description:
      "Use this when the user provides gene(s) and asks for drugs that interact. Only the top ranked interactions will be returned. Drugs will be first sorted by FDA approval, then Interaction Score.",
    inputSchema: {geneNames: z.array(z.string()).min(1).max(25).describe("HGNC gene symbols only")},
    annotations: {
      destructive: false,
      idempotent:  true,
      cacheable:   false,
      world_interaction: "open",
      side_effects: ["external_api_calls"],
      resource_usage: "network_io_heavy",
    },

    async handler(
  		{ geneNames }: { geneNames: string[] }  // <- matches schema
	  ) {

		const query = /* GraphQL */ `
      query genes($names: [String!]) {
          genes(names: $names) {
              nodes {
                  name
                  interactions {
                      drug {
                          name
                          approved
                      }
                      interactionScore
                      interactionTypes {
                          type
                          directionality
                      }
                      publications {
                        pmid
                      }
                      sources {
                        sourceDbName
                        baseUrl
                      }
                  }
              }
          }
      }`;

      const rawNames = (geneNames ?? [])
        .map(s => (s ?? "").trim())
        .filter(Boolean);

      //const rawNames = geneNames.split(",").map(s => s.trim()).filter(Boolean);
      const normalizedNames = rawNames.map(raw =>
        normalizeEntityFast(raw, getGeneIndex(), 0.7) || raw
      );

      const totalBudget = normalizedNames.length === 1 ? 15 : 100;

			const res = await fetch("https://dgidb.org/api/graphql", {
				method: "POST",
				headers: { "Content-Type": "application/json", ...API_CONFIG.headers },
				body: JSON.stringify({ query, variables: { names: normalizedNames } })
			}).then(r => r.json()) as {
				data?: { genes?: { nodes: { name: string; interactions: any[] }[] } };
				errors?: unknown[];
			};

    const instructions =
          "Interaction Score: Scoring metric based on the evidence supporting an interaction.\n" +
          "Interaction Type:  Nature of the association between a particular drug and gene.\n" +
          "Interaction Direction: Describes the effect of the interaction on the biological activity of the gene.\n" +
          "PMID url: Direct link to the PubMed publication supporting this interaction.\n" +
          "The returned list of interactions is limited to the top 15 and does not include all of them.\n" +
          "When returning information to users you MUST cite URLs used for specific information. Always cite the PubMed IDs and Database Sources as links";

      // If no data, wrap in a content-based return for MCP compliance
      if (!res.data?.genes) {
        return {
        isError: true,
        content: [{ type: "text" as const, text: JSON.stringify(res, null, 2) }]
        };
      }

      const nodes = res.data.genes.nodes ?? [];

      if (!nodes.length) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: `No gene nodes found for: ${normalizedNames.join(", ")}` }]
        };
      }

      const selected_nodes = selectNodes(nodes, normalizedNames, totalBudget)

      const payload = {
          instructions,
          "API Results": convertNodesPMID(selected_nodes),
        };

      return {
        isError: false,
        content: [{ type: "text" as const, text: JSON.stringify(payload, null, 2) }]
      };
    },
    },

} as const;


// -------------------------------------------------------------
// MCP SERVER (only the two fixed tools)
// -------------------------------------------------------------
class DgidbMCP extends McpAgent {
  server = new McpServer({
    name:        API_CONFIG.name,
    version:     API_CONFIG.version,
    description: API_CONFIG.description,},
    {
      instructions: `Use the tools to find drug-gene interactions with the Drug-Gene Interaction Database (DGIdb).`,
    });

  async init() {
    /* register fixed-schema tools */
    const { getGeneInteractionsForDrug, getDrugInteractionsForGene, getDrugInfo, getGeneInfo } = tools;

    this.server.tool(getDrugInfo.name, getDrugInfo.description, getDrugInfo.inputSchema, getDrugInfo.handler);
    this.server.tool(getGeneInfo.name, getGeneInfo.description, getGeneInfo.inputSchema, getGeneInfo.handler);

    this.server.tool(getGeneInteractionsForDrug.name, getGeneInteractionsForDrug.description, getGeneInteractionsForDrug.inputSchema, getGeneInteractionsForDrug.handler);
    this.server.tool(getDrugInteractionsForGene.name, getDrugInteractionsForGene.description, getDrugInteractionsForGene.inputSchema, getDrugInteractionsForGene.handler);
  }
}

// -------------------------------------------------------------
// CLOUDFLARE WORKER RUNTIME (SSE only)
// -------------------------------------------------------------
interface ExecutionContext {
  waitUntil(promise: Promise<unknown>): void;
  passThroughOnException(): void;
}


export default {
  async fetch(
    request: Request,
    env: Env,
    ctx: ExecutionContext
  ): Promise<Response> {
    const url = new URL(request.url);
 
    console.log(Object.keys(env));
 
    if (url.pathname === "/") {
      const accept = request.headers.get("Accept") ?? "";
      if (accept.includes("text/html")) {
        return new Response(`<!doctype html>
    <html>
      <head>
        <meta charset="utf-8" />
        <title>DGIdb MCP</title>
      </head>
      <body>DGIdb MCP</body>
    </html>`, { headers: { "Content-Type": "text/html; charset=utf-8" } });
      }
 
      return new Response(
        `${API_CONFIG.name} – MCP Server ${API_CONFIG.version}. Use /mcp (Streamable HTTP) or /sse (legacy).`,
        { status: 200, headers: { "Content-Type": "text/plain" } }
      );
    }
 
    if (url.pathname === "/privacy") {
      const html = `<!doctype html>
      <html><head><meta charset="utf-8"><title>Privacy Policy</title></head>
      <body style="font-family:system-ui;max-width:760px;margin:40px auto;padding:0 16px;line-height:1.5">
        <h1>Privacy Policy</h1>
        <p><strong>What this service does:</strong> This MCP server receives query terms you provide (e.g., gene names, drug names) and forwards them to the DGIdb GraphQL API to retrieve drug-gene interaction results.</p>
        <p><strong>Data we process:</strong> The text inputs you submit to the tools and standard request metadata.</p>
        <p><strong>How we use data:</strong> Only to fulfill your request and return DGIdb results.</p>
        <p><strong>Data sharing:</strong> Requests are sent to (1) Cloudflare (hosting/infrastructure) and (2) the DGIdb API at dgidb.org.</p>
        <p><strong>Data retention:</strong> We do not intentionally store your tool inputs in an application database. Cloudflare Workers Observability is enabled and may retain operational logs/telemetry for debugging and service reliability.</p>
        <p><strong>Security:</strong> Data is transmitted over HTTPS/TLS.</p>
        <p><em>Last updated: 2026-04-14</em></p>
      </body></html>`;
      return new Response(html, { headers: { "Content-Type": "text/html; charset=utf-8" } });
    }
 
    /* ────────────────────────────────────────────────
       NEW: Streamable HTTP transport (/mcp)
    ─────────────────────────────────────────────────*/
    if (url.pathname === "/mcp" || url.pathname.startsWith("/mcp/")) {
      const protocolVersion = request.headers.get("MCP-Protocol-Version");
 
      // @ts-ignore – serve helper is mixed-in by DgidbMCP (Streamable HTTP)
      const response = await DgidbMCP.serve("/mcp").fetch(request, env, ctx);
 
      if (protocolVersion && response instanceof Response) {
        const headers = new Headers(response.headers);
        headers.set("MCP-Protocol-Version", protocolVersion);
        return new Response(response.body, {
          status: response.status,
          statusText: response.statusText,
          headers
        });
      }
      return response;
    }
 
    /* ────────────────────────────────────────────────
       Legacy SSE transport (kept for now)
    ─────────────────────────────────────────────────*/
    if (url.pathname === "/sse" || url.pathname.startsWith("/sse/")) {
      const protocolVersion = request.headers.get("MCP-Protocol-Version");
 
      // @ts-ignore – serveSSE helper is mixed-in by DgidbMCP (SSE)
      const response = await DgidbMCP.serveSSE("/sse").fetch(request, env, ctx);
 
      if (protocolVersion && response instanceof Response) {
        const headers = new Headers(response.headers);
        headers.set("MCP-Protocol-Version", protocolVersion);
        return new Response(response.body, {
          status: response.status,
          statusText: response.statusText,
          headers
        });
      }
      return response;
    }
 
    return new Response(
      `${API_CONFIG.name} – MCP Server ${API_CONFIG.version}. Use /mcp (Streamable HTTP) or /sse (legacy).`,
      { status: 200, headers: { "Content-Type": "text/plain" } }
    );
  }
};

// Export the class for tests if you like
export { DgidbMCP };