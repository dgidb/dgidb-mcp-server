import pandas as pd
import sys, requests
import re
import csv
import os
from host_dgidb_civic_MCP import get_drug_interactions_for_gene_list 

civic_df = pd.read_csv('data/CIViC_evidence_extracts_curators_9_6_25.csv')

res_rows = civic_df.groupby(['molecularProfile_name','feature_name', 'disease_name', 'therapies']).filter(lambda g: set(g['significance']) == {'RESISTANCE'} and set(g['evidenceDirection']) == {'SUPPORTS'})
filtered_df = res_rows[['molecularProfile_name','feature_name', 'disease_name', 'therapies', 'evidenceType', 'significance']].drop_duplicates()

res_df = res_rows[['molecularProfile_name','feature_name', 'disease_name', 'therapies', 'evidenceType', 'significance']]

HEADERS = {"Content-Type": "application/json"}

def run_civic_query(query, variables):
    #print("Running GraphQL Query with variables:", variables, file=sys.stderr)
    r = requests.post(
        "https://civicdb.org/api/graphql",
        json={"query": query, "variables": variables},
        headers=HEADERS,
        timeout=30,
    )
    #print("GraphQL Response:", r.text, file=sys.stderr)
    r.raise_for_status()
    return r.json()

evidence_query = """
query evidenceItems($molecularProfileName: String, $diseaseName: String, $therapyName: String, $significance: EvidenceSignificance ) {
    evidenceItems(molecularProfileName: $molecularProfileName, diseaseName: $diseaseName, therapyName: $therapyName, significance: $significance) {
        nodes { 
            status 
            evidenceType
            evidenceDirection 
            significance
            molecularProfile{
                variants{
                    name
                    feature{
                        name
                    }
                }
            }
            disease{
                displayName
            }
            therapies{
                name
            }
            variantOrigin 
            description
            evidenceLevel
            evidenceRating
            id 
        }
    }
}
"""


gene_query = """
query genes($names: [String!]) {
    genes(names: $names) {
        nodes {
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
            }
        }
    }
}
"""

def run_dgidb_query(query, variables: dict) -> dict:
    #print("Running GraphQL Query with variables:", variables, file=sys.stderr)
    r = requests.post(
        "https://dgidb.org/api/graphql",
        json={"query": query, "variables": variables},
        headers=HEADERS,
        timeout=30,
    )
    #print("GraphQL Response:", r.text, file=sys.stderr)
    r.raise_for_status()
    return r.json()

def parse_list(s: str):
    # split on comma or whitespace (one or more), strip empties
    return [d.strip() for d in re.split(r'[,\s]+', s.strip()) if d.strip()]


FULL_DATASET_PATH = 'data/full_joint_task_dataset.csv'
SUBSET_DATASET_PATH = 'data/per_gene_ranked_joint_task_drug_lists_dataset.csv'

already_checked_disease_therapies = set()

if not os.path.exists(FULL_DATASET_PATH):
    with open(FULL_DATASET_PATH, 'w', newline='') as csvfile: 
        writer = csv.writer(csvfile)
        writer.writerow(['molecularProfile_name', 'disease_name', 'therapies', 'gene_list', 'etype', 'significance'])
        
        for idx, row in res_df.iterrows():
            MP, gene, disease, therapies, etype, significance = row[['molecularProfile_name', 'feature_name', 'disease_name', 'therapies', 'evidenceType','significance']].values
            #get all genes for this combo of diease + therapies

            #print(gene, disease, therapies)

            if (disease, therapies) in already_checked_disease_therapies:
                continue

            already_checked_disease_therapies.add((disease, therapies))

            if ',' in therapies or any(s in MP for s in [':', ',', 'AND', 'OR']): #dont get any disease + therapy that result in this
                continue

            resp = run_civic_query(evidence_query, {'diseaseName':disease, 'therapyName':therapies, 'significance': 'RESISTANCE'})

            nodes = resp['data']['evidenceItems']['nodes']

            if not nodes:
                print('missing nodes')
                continue

            sub_genes = set()

            if any(len(node['molecularProfile']['variants']) > 1 for node in nodes):
                #print('complex')
                continue

            for node in nodes:
                gene = node['molecularProfile']['variants'][0]['feature']['name']
                sub_genes.add(gene)
            
            halt = False
            for sub_gene in sub_genes:
                if '::' in sub_gene:
                    halt = True
            if halt:
                #print('fusion')
                continue

            writer.writerow([MP, disease, therapies, ','.join(sub_genes), etype, significance])

df = pd.read_csv(FULL_DATASET_PATH)
df = df.drop_duplicates(subset=['molecularProfile_name'], keep='first')

df = df.sample(n=50)

#This is the starting point script used for prompting
if not os.exists('data/subset_50_ranked_joint_task_drug_lists_dataset.csv'):
    df.to_csv('data/subset_50_ranked_joint_task_drug_lists_dataset.csv')

FIELDNAMES = ['molecularProfile_name','disease_name','therapies','evidenceType','significance', 'level', 'rating', 'label_CIViC_genes', 'label_DGIdb_drugs']

#Writes the order of genes that should be returned
#and then the labels for what drugs should be returned for each gene
#the molecular profile is an artifact and not used directly by the task.

if not os.path.exists(SUBSET_DATASET_PATH):
    with open(SUBSET_DATASET_PATH, "a", newline="", encoding="utf-8") as f: 
        writer = csv.writer(f)
        writer.writerow(FIELDNAMES)
        for idx, row in df.iterrows():
            molecularProfile_name, disease_name, therapies, evidenceType, significance = row[['molecularProfile_name','disease_name','therapies','evidenceType','significance']]

            print('-----------------------------------------')
            print(molecularProfile_name, disease_name, therapies)

            resp = run_civic_query(evidence_query, {'diseaseName':disease_name, 'therapyName':therapies, 'significance': 'RESISTANCE'})
            nodes = resp['data']['evidenceItems']['nodes']

            if not nodes:
                print('missing nodes')
                continue

            gene_list = []
            for node in nodes:
                gene = node['molecularProfile']['variants'][0]['feature']['name']
                level = node['evidenceLevel']
                rating = node['evidenceRating']
                gene_list.append((gene, level, rating))

            seen = set()
            deduped = []
            for gene, level, rating in gene_list:
                if gene not in seen:
                    seen.add(gene)
                    deduped.append((gene, level, rating))

            gene_list = deduped

            for gene, level, rating in gene_list:
                print(gene_list)
                drug_nodes = get_drug_interactions_for_gene_list(gene)

                if not drug_nodes:
                    print('MISSING drug nodes')
                    writer.writerow([molecularProfile_name, disease_name, therapies, evidenceType, significance, level, rating, gene, ''])
                    continue

                drug_list = []
                for drug in drug_nodes:
                    interactions = drug_nodes[drug]
                    for i in interactions:
                        d = i['drug']['name']
                        drug_list.append(d)

                drug_list = list(dict.fromkeys(drug_list))

                writer.writerow([
                    molecularProfile_name, disease_name, therapies,
                    evidenceType, significance, level, rating,
                    gene, ','.join(drug_list)
                ])