# README

## 1. Convert all Drug IDs to CHEMBL ID (./scripts/01-studyDrugs_in_ChEMBL.R)

Since each study uses a different type of id or alias for the drugs, we first need to convert all the information into one type of ID. The final format we want to have is CHEMBL IDs as DGIdb uses CHEMBL. 

### 1.1. Map drug Names to PubCHEM ID
Although the conversion between different IDs is facilitated using UniChem, there is no universal tool to map drug names to CHEMBL IDs. We choose to first convert names to pubchem IDs and then convert them to CHEMBL because PubChem API allows search with partial matches and aliases.

#### 1.1.1 Using name2CID function
This function uses PubChem pug-rest API and searches for the drug in PubChem database to retrieve its CID. For each dataset, we created drugList_CID.tsv under its own directory under ./data/processed

1. **Aliper2016**: All drugs in Aliper 2016 list are mapped to CIDs.
2. **Barardo2017**: All drugs in Barardo 2017 list are mapped to CIDs.
3. **Calvert2016**: All drugs are mapped to CIDs.
4. **Donertas2018**: Except for "15-delta prostaglandin J2", all other drugs are mapped to the CIDs.
5. **Yang2018**: We could find the CID for GW3965. Although there is no entry for the conjugated linoleic acid, we included it in the list as 'linoleic acid'. There is no CID for "bis(2-ethylhexyl)".
6. **Liu2016**: The list was already given as PubChem CIDs, however, to be consistent with the other studies we converted the compounds to their flat form.
7. **Fernandes2016**: 4 compounds were not mapped from names to CID (FK‐228, CHR‐3996, GDC‐0068, MK‐2206). One additional compound was not mapped from CID to chembl.

### 1.2 Convert IDs to CHEMBL IDs using UniChem

#### 1.2.1 unichem function

This function uses UniChem API to map compounds between different sources. We used this function to map directly identifiers from drugbank and pdb to chembl ids. For each dataset, we created drugList_CHEMBL.tsv under its own directory under ./data/processed

8. **Snell2016**: 24 out of 27 drugbank identifiers were successfully mapped to chembl using unichem. Gamma Hydroxybutyric Acid (DB01440), Nafarelin (DB00666) and Colistin (DB00803) were mapped manually.
9. **Snell2018**:  We were not able to map 4 compounds automatically from drugbank to chembl using to unichem: colistin      (DB00803), salicylate-sodium (DB01398), ivermectin (DB00602), gadopentetate dimeglumine (DB00789). We mapped these compounds manually.
10. **Mofidifar2018**: We couln't find a drug named Nacitentan in any database (probably Macitentan was named incorrectly). All the other compounds were mapped to chembl.
11. **Ziehm2017**: All compounds were mapped from pdb to chembl
12. **Fuentealba2018**: All compounds are mapped from drugbank to chembl using unichem.

#### 1.2.2 using src1src2.txt unichem conversion data

At the time we write the scripts, there was a problem with the UniChem API in mapping PubChem CIDs to CHEMBL IDs. We instead downloaded the conversion file from [ftp](link here) and used this for conversion. cid2chembl function makes use of this file. The resulting list include the rows with NAs which represent the IDs we cannot map to CHEMBL. Although there are many CIDs that are not mapped to Chembl ID, we calculated how many of the initial drug names do not have any CHEMBL id. Only three datasets have some drugs missing: 


The whole list of drugs that are missing ChEMBL IDs are given as './data/processed/noCHEMBL.tsv'

## 2. Compile Drug-Target Information Using DGIdb (./scripts/02-drugs2targets.R)

<hr>
**Drugs missing ChEMBL IDs**

1. **Barardo2017**: 5 drugs cannot be matched: mmk-1, rdp-58, ro 25-1392, gv1001, cardiolipin.
2. **Donertas2018**: 1 drug cannot be matched: quinostatin
3. **Liu2016**: 110 drugs cannot be matched
<hr>

### Retrieve Interactions

We downloaded [DGIdb Interaction File(http://www.dgidb.org/data/interactions.tsv) and mapped to studies using the CHEMBL IDs we compiled. Among 278 unique Chembl IDs from the studies, we found interactions for 163 (58.6%). In total there are 796 unique genes that are targeted by at least one drug found in at least one of the studies. Study - gene - drug interaction file is saved as './data/processed/combinedList.tsv'. For each study, the number of drugs and genes are as follows:

|study          | numDrugs| numGenes|
|:--------------|--------:|--------:|
|Aliper2016     |        7|      123|
|Barardo2017    |       11|      113|
|Calvert2016    |        7|      102|
|Donertas2018   |       17|      261|
|Fernandes2016  |       14|       92|
|Fuentealba2018 |       10|      164|
|Liu2016        |       54|      246|
|Mofidifar2018  |        4|        5|
|Snell2016      |       21|      150|
|Snell2018      |       23|      140|
|Yang2018       |        2|        8|
|Ziehm2017      |       11|      103|

### Create Adjacency Matrices for genes, drugs, and studies

Using the data we generated in the previous steps we constructed matrices for the interactions between 1) geneXdrug: ./data/processed/gene_drugMat.tsv, 2) studyXgene: ./data/processed/study_geneMat.tsv, and 3) studyXdrug: ./data/processed/study_drugMat.tsv
Files are also saved as .rds files to ease the upload to rsession.

## 3. Literature Information for Drugs & Targets (03-literatureInfo.R)

### 1. Convert drugIDs in DrugAge to CHEMBL IDs

The drugs with no PubChem CID are saved as './data/processed/drugAge_noCID.tsv' and the resulting table with CHEMBL IDs for all drugs that can be mapped is saved as './data/processed/drugAgeCHEMBL.tsv'.

### 2. Create a literature matrix for Drugs

Create a matrix where the rows are Human, Model and DrugAge and the columns are the drugs either discovered in any of the studies or in DrugAge. If a drug is in DrugAge, the matrix has a value of 1, 0 otherwise. If the genes in GenAge Human or Model data are targeted by a given drug, then it gets a value of 1 for the Human or Model row. ('./data/processed/lit_drugMat.rds')

### 3. Create a literature matrix for Genes

Create a matrix where the rows are Human, Model and DrugAge and the columns are the genes in either genage or targeted by the drugs discovered in any of the studies. If a gene is in GenAge Human or Model data, the matrix has a value of 1, 0 otherwise. If the gene is among the targets of DrugAge drugs, then it gets a value of 1 for the DrugAge row. ('./data/processed/lit_geneMat.rds')

## 4. Create matrix for pathways (04-pathwayMat.R)

Create matrices for genage/drugage and studies such that each column corresponds to a KEGG pathway and the numbers in the matrices are the % genes discovered among the druggable genes, i.e. genes targeted by at least one drug in DGIdb. There were 11 pathways with no druggable gene, thus excluded from the analysis. We continued with 319 pathways.

## 5. Druggable Genome Distribution (05-druggableDist.R)
Result: ./results/druggableGenomeDist.pdf

## 6. Study-Drug Network (06-studyDrugNet.R)
Result: ./results/st_drug_net.pdf

## 7. Study - Literature Overview as a Heatmap (07-circosPlot.R)
Result: ./results/circularHeatmap.pdf

## 8. Create lists of studies, genes, drugs, pathways, and literature as a reference (08-createMasterLists.R)

* './data/processed/study_drug_gene_literature_masterList.tsv'
* './data/processed/study_pathway_literature_masterList.tsv'
