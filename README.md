# DGIdb MCP Server

Preprint: https://www.biorxiv.org/content/10.64898/2026.01.18.700012v1

This is a Cloudflare Workers-based Model Context Protocol (MCP) server that provides tools for querying the DGIdb (Drug Gene Interaction Database) API. 

DGIdb streamlines the search for druggable therapeutic targets through the aggregation, categorization, and curation of drug and gene data from publications and expert resources. 

## Example of DGIdb MCP Server on Claude Desktop

### User Prompt: 
"What genes interact with the drug Ibrutinib and by what mechanism of interaction according to DGIdb?"

  <summary>Tool Called: <code>get_gene_interactions_for_drug_list</code></summary>

  <details>
    <summary><b>Request</b></summary>

```json
{
  "drugNames": "Ibrutinib"
}
```
</details> 
<details> <summary><b>API Response</b></summary>
  
```json
{
  "Instructions": "Interaction Score: Scoring metric based on the evidence supporting an interaction.
                  Interaction Direction: Describes the effect of the interaction on the biological activity of the gene.
                  Interaction Type:  Nature of the association between a particular drug and gene.
                  When returning information to users you MUST cite URLs used for specific information."

  "API Results": [
    { "gene": "PLCG2",
      "interactionScore": 4.61,
      "interactionTypes": [],
      "sources": ["CGI", "DoCM"],
      "pmidUrls": ["https://pubmed.ncbi.nlm.nih.gov/24869598/"]},

    { "gene": "MYD88",
      "interactionScore": 1.0,
      "interactionTypes": [],
      "sources": ["CIViC", "PharmGKB", "CGI", "DoCM"],
      "pmidUrls": []
    },
    { "gene": "BTK",
      "interactionScore": 0.38,
      "interactionTypes": [
          {
            "type": "inhibitor",
            "directionality": "INHIBITORY"
          }
      ],
      "sources": ["CIViC", "PharmGKB", "CGI", "DoCM", "TTD", "MyCancerGenome", "TdgClinicalTrial", "ChEMBL", "COSMIC"],
      "pmidUrls": [
          "https://pubmed.ncbi.nlm.nih.gov/25222877/",
          "https://pubmed.ncbi.nlm.nih.gov/24869598/",
          "https://pubmed.ncbi.nlm.nih.gov/27199251/"
        ]
    },
    ...
  ]
}
```
  
</details>


<details>
    <summary><b>Claude Response</b></summary>
  

According to DGIdb, **Ibrutinib** interacts with multiple genes. Here are the top interactions ranked by interaction score:

**Primary Target**

| Gene | Interaction Type | Direction | Score | Sources |
|------|------------------|-----------|------:|---------|
| BTK | Inhibitor | Inhibitory | 0.38 | CGI, CIViC, PharmGKB, ChEMBL, and others |

BTK (Bruton's Tyrosine Kinase) is the **primary therapeutic target** of ibrutinib. The drug acts as an irreversible inhibitor, blocking BTK's enzymatic activity. Supporting publications: [PMID: 25222877](https://pubmed.ncbi.nlm.nih.gov/25222877/), [PMID: 24869598](https://pubmed.ncbi.nlm.nih.gov/24869598/), [PMID: 27199251](https://pubmed.ncbi.nlm.nih.gov/27199251/)

**Resistance/Response Biomarkers**

| Gene | Score | Sources |
|------|------:|---------|
| PLCG2 | 4.61 | CGI, DoCM — [PMID: 24869598](https://pubmed.ncbi.nlm.nih.gov/24869598/) |
| MYD88 | 1.54 | CGI, CIViC, PharmGKB |

PLCG2 and MYD88 are associated with ibrutinib resistance or response prediction rather than being direct drug targets.

**Additional Inhibitory Targets**

| Gene | Interaction Type | Score | Source |
|------|------------------|------:|--------|
| ENUR2 | Inhibitor | 0.31 | Guide to Pharmacology |
| BMX | Not specified | 0.26 | TTD |
| EFNA1 | Inhibitor | 0.14 | Guide to Pharmacology |
| EPS15 | Inhibitor | 0.13 | Guide to Pharmacology |
| TEC | Inhibitor | 0.11 | Guide to Pharmacology |
| EPHA8 | Inhibitor | 0.11 | Guide to Pharmacology |

Ibrutinib is known to have off-target kinase inhibition beyond BTK, which may contribute to both efficacy and side effects.

View the original chat here: [https://claude.ai/share/8e9224cd-4927-4337-9a76-cc4ab8f95d93 ](https://claude.ai/share/1297b4d8-d55b-408d-9648-ef34e266f7ca)

</details></details></details>

## Configuring Claude Desktop

Example video taken from CIViC MCP

https://github.com/user-attachments/assets/5890f79a-e2fc-49f6-b5f4-ef191d07872d

Install Node.js (https://nodejs.org/)

Click "LTS" (Recommended for Most Users) — this gives you Node.js and npx
Download and install it like any normal app

Once installed:
On Windows: Open “Command Prompt” or “PowerShell”
On macOS: Open “Terminal”

Then run:
```bash
node -v
npx -v
```

Confirm that both give versions.

## Accessing DGIdb via MCP
Add this configuration to your `claude_desktop_config.json` file:

```json
{
  "mcpServers": {
    "lars-dgidb-mcp-server": {
    "command": "npx",
    "args": [
          "mcp-remote",
          "https://dgidb-mcp-server.larscivic.workers.dev/mcp"
        ]
    }
  }
}
```

## Joint Access to DGIdb and CIViC via MCP
Add this configuration to your `claude_desktop_config.json` file:

```json
{
"mcpServers": {
    "lars-civic-mcp-server": {
    "command": "npx",
    "args": [
          "mcp-remote",
          "https://civic-mcp-server-v2.larscivic.workers.dev/mcp"
        ]
      },

    "lars-dgidb-mcp-server": {
    "command": "npx",
    "args": [
          "mcp-remote",
          "https://dgidb-mcp-server.larscivic.workers.dev/mcp"
        ]
      }
  }
}
```

## Usage

Once configured, restart Claude Desktop. The server provides 4 main tools:

1. **`get_drug_info`**: Gets drug info including approval, if used in immunotherapy, and other drug attributes for a list of drugs.
2. **`get_gene_info`**: Gets gene category info for a list of genes.
3. **`get_drug_interactions_for_gene_list`**: Gets drugs that interact with a list of genes.
4. **`get_gene_interactions_for_drug_list`**: Gets genes that interact with a list of drugs.
