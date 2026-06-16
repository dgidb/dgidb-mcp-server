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

## Usage

Once configured, restart Claude Desktop. The server provides 4 main tools:

1. **`get_drug_info`**: Gets drug info including approval, if used in immunotherapy, and other drug attributes for a list of drugs.
2. **`get_gene_info`**: Gets gene category info for a list of genes.
3. **`get_drug_interactions_for_gene_list`**: Gets drugs that interact with a list of genes.
4. **`get_gene_interactions_for_drug_list`**: Gets genes that interact with a list of drugs.

### Example Use Cases

| Category | Example Query | Intended Output* |
| --- | --- | --- |
| Drug Discovery | Which FDA-approved drugs target the gene KIT? | A ranked list of FDA-approved KIT-interacting drugs with interaction scores, interaction types, approval status, and links to supporting sources. [Chat Link](https://claude.ai/share/62d8e788-7b0f-453a-af08-dc066d34b9b1). |
| Drug Info | What information is available about Imatinib? | A structured overview of Imatinib including its FDA approval status, year of approval, therapeutic and drug class annotations. [Chat Link](https://claude.ai/share/999dfa2e-aa4b-4a5c-bcf7-9be72886b873). |
| Gene Category Info | What categories and functional annotations are associated with the gene BRAF? | A structured summary of BRAF annotations, including clinical actionability, drug resistance relevance, druggability, enzyme class, and kinase classification, with the supporting source databases listed for each category. [Chat Link](https://claude.ai/share/91a62f7d-d278-4a4b-a70c-8f54d598d147). |
| Gene Symbol Disambiguation | Is FLT3 an unambiguous gene symbol? | An assessment of whether FLT3 maps uniquely to a single gene, including any known aliases/synonyms, overlapping symbols, and the resolved canonical gene identifier(s), with supporting source links. [Chat Link](https://claude.ai/share/5fa0883b-538b-47c8-8b68-99b8d9db4045). |
| Mechanism-specific targeting | Which inhibitors interact with EGFR? | A filtered list of EGFR inhibitors, prioritized by FDA approval and interaction score, with evidence sources cited. [Chat Link](https://claude.ai/share/0daf9c70-b6cd-4c00-b2af-ff2c91a6aff7). |
| Interaction type interpretation | What genes interact with the drug Ibrutinib and by what mechanism of interaction? | A structured summary of Ibrutinib-gene interactions, annotated with interaction direction and type (e.g., inhibitor, modulator). [Chat Link](https://claude.ai/share/1297b4d8-d55b-408d-9648-ef34e266f7ca). |

Chat links to conversations with Opus 4.5 via Claude Desktop on 1/14/2026.

*API responses are not visible in the linked chats.

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

### Joint Use Cases (DGIdb + CIViC)

| Category | Example Query | Intended Output* |
| --- | --- | --- |
| Resistance-guided drug discovery (CIViC→DGIdb) | What genes can cause resistance to ibrutinib in chronic lymphocytic leukemia, and what alternative drugs can target them? | CIViC-supported resistance genes for ibrutinib in CLL and a ranked list of interacting alternative drugs from DGIdb, with citations. [Chat Link](https://claude.ai/share/be18a9b9-0a1b-4831-9dc9-9f72dd9d5049). |
| Oncogenicity-driven functional validation (CIViC→DGIdb) | What is the evidence for oncogenicity of variants in ERBB2 for breast cancer, and what inhibitors could I use to experimentally validate their importance in vitro? | CIViC-supported oncogenic ERBB2 variants in breast cancer with summarized evidence and citations, plus a ranked list of ERBB2-targeting inhibitors from DGIdb suitable for in vitro validation. [Chat Link](https://claude.ai/share/1bd81002-fe12-44aa-9b5c-775a0fb426d2). |
| Drug target profiling with clinical evidence (DGIdb→CIViC) | What genes interact with the drug Ibrutinib, and which of those genes have specific variants in CIViC associated with sensitivity or resistance? | A list of genes interacting with ibrutinib annotated with CIViC sensitivity or resistance evidence where available. [Chat Link](https://claude.ai/share/9ad52a71-74c1-40a5-aa22-afc03e5fa940). |
| Target resistance interpretation (DGIdb→CIViC) | Are there known resistance mechanisms for genes targeted by Osimertinib? | Osimertinib target genes with CIViC-documented resistance mechanisms and supporting evidence links. [Chat Link](https://claude.ai/share/9da56313-74b4-4978-954b-43437a63e47f). |
| Gene-guided therapy matching (DGIdb→CIViC) | Which drugs interact with PIK3CA, and is there evidence of sensitivity to these drugs in breast cancer? | PIK3CA-interacting drugs prioritized by DGIdb and annotated with CIViC sensitivity evidence in breast cancer. [Chat Link](https://claude.ai/share/b73748ee-74d6-4b7f-92f1-8a8620960bcd). |
| Expression-based therapy matching (DGIdb→CIViC) | For tumors with EGFR amplification, which drugs target EGFR and is there clinical evidence supporting sensitivity? | EGFR-targeting drugs with CIViC evidence supporting sensitivity in EGFR-amplified tumors. [Chat Link](https://claude.ai/share/b2d6ffa7-a2a2-4c2f-9805-3d6aad577362). |

Chat links to conversations with Opus 4.5 via Claude Desktop on 1/14/2026 and 3/25/2026.

*API responses are not visible in the linked chats.

## GPT-5 Evaluation

GPT_eval/src details the scripts used for the published evaluation. The randomly sampled input files used for this are available. The full experiment can be run with this command:

```bash
bash run_experiments.sh
```
