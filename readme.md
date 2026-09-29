# TCGA-Subtypes

**Integrative Machine Learning and Network Analysis Extend TCGA Subtypes and Reveal Pan-Cancer Functional Convergence**

This repository provides the TCGA molecular subtype assignments used in the study, including subtype annotations extended by machine-learning prediction for samples without previously assigned subtypes.

- Repository: https://github.com/CBIIT-CGBB/TCGA-Subtypes
- TCGA subtype data: https://github.com/CBIIT-CGBB/TCGA-Subtypes/tree/main/TCGA_Subtype
- OmicPath R package: https://github.com/CBIIT-CGBB/OmicPath

## Overview

Molecular subtypes capture important biological heterogeneity within individual cancer types, but subtype annotations are incomplete for a substantial fraction of TCGA samples. In this study, we used transcriptomic data and supervised machine learning to extend existing TCGA subtype annotations, then integrated gene-level, pathway-level, and network-level analyses to investigate how subtype biology is organized within and across cancer types.

The analysis had three major goals:

1. **Extend TCGA subtype coverage** by predicting missing subtype assignments from RNA-seq expression profiles.
2. **Characterize subtype-associated genes and pathways** and compare subtype relationships at the gene and pathway levels.
3. **Resolve pathway-level convergence into sub-pathway structure** to determine whether tumor subtypes sharing the same pathway use the same or different internal network modules.

## Study summary

The analysis included **24 TCGA tumor types**, representing **10,468 samples** and **109 molecular subtypes**.

- **8,428 samples** had existing subtype annotations.
- **2,040 samples** across **20 tumor types** required subtype prediction.
- Four tumor types with complete subtype assignments did not require prediction but were retained for downstream pan-cancer analyses.

Five supervised machine-learning algorithms were evaluated:

- `glmnet` — regularized generalized linear model
- `knn` — k-nearest neighbors
- `rf` — random forest
- `svmLinear` — linear support vector machine
- `xgbTree` — extreme gradient boosting

Across subtype–method combinations, the study observed a mean AUC of approximately **0.887** and a median AUC of approximately **0.954**.

## Performance-Weighted Subtype Voting (PWSV)

To integrate predictions across machine-learning methods, we used a **Performance-Weighted Subtype Voting (PWSV)** framework.

For sample \(i\), candidate subtype \(k\), and model \(m\), the model-predicted probability is weighted by the corresponding predictive performance:

$$w_{km} = AUC_{km} \times I(AUC_{km} \ge 0.8)$$

$$ S_{ik} = \sum_m w_{km} p_{ikm} $$

where:

- $$p_{ikm}$$ is the probability assigned to subtype \(k\) for sample \(i\) by model \(m\);
- $$AUC_{km}$$ is the corresponding model/subtype AUC;
- models with AUC < 0.8 receive zero weight.

The subtype with the largest aggregated score is assigned as the final predicted subtype:

$$\widehat{k}_i = \arg\max_k S_{ik}$$

This strategy combines prediction confidence with model-specific discriminative performance.

## Repository contents

### `TCGA_Subtype/`

The [`TCGA_Subtype`](https://github.com/CBIIT-CGBB/TCGA-Subtypes/tree/main/TCGA_Subtype) directory contains the TCGA molecular subtype assignments distributed with this study.

These files are intended to provide a convenient subtype resource for downstream TCGA analyses, including:

- subtype-level gene-expression analysis;
- pathway and gene-set enrichment analysis;
- survival and clinical association studies;
- pan-cancer subtype comparisons;
- external validation or benchmarking studies.

Users should preserve TCGA sample identifiers when joining these subtype annotations to molecular or clinical data.

### `HALLMARK_Cluster/`

The [`HALLMARK_Cluster`](https://github.com/CBIIT-CGBB/TCGA-Subtypes/tree/main/HALLMARK_Cluster) directory contains the HALLMARK pathway network-clustering results used to define the sub-pathways analyzed in this study.

Gene-link information used to construct the pathway networks was obtained from **NeST**:

- NeST: https://idekerlab.ucsd.edu/nest/

For each HALLMARK pathway, genes were connected using the available NeST gene-link information. A connected network backbone was then generated using a **minimum spanning tree (MST)** approach. The resulting graph was partitioned into topological sub-pathways using **Louvain community detection**. Each Louvain community was treated as a pathway sub-cluster (sub-pathway) for downstream analyses.

For each HALLMARK pathway, three output file types are provided:

| File type | Description |
|---|---|
| `*_edge.csv` | Gene–gene relationships used to define the pathway network edges. |
| `*_node.csv` | Node-level information for pathway genes, including plotting coordinates and the assigned Louvain cluster/sub-pathway ID. |
| `*.pdf` | Network visualization of the HALLMARK pathway, with nodes colored according to their cluster/sub-pathway assignment. |

Together, these files provide both the machine-readable network representation and the corresponding graphical view of each HALLMARK pathway decomposition.

## Downstream subtype analysis

### Subtype-associated genes

Within each tumor type, subtype-associated genes were identified using one-subtype-versus-other-subtypes generalized linear models with covariate adjustment. Genes passing the study-defined false-discovery-rate threshold were used for downstream pathway analyses.

### HALLMARK pathway analysis

Subtype-associated genes were evaluated against MSigDB HALLMARK gene sets to identify recurrent biological programs across molecular subtypes.

The study identified pathway-level relationships involving processes such as:

- cell-cycle regulation;
- immune and interferon signaling;
- angiogenesis;
- metabolism;
- MYC/E2F-associated programs.

### Gene-level and pathway-level clustering

Subtype relationships were evaluated using hierarchical clustering based on cosine similarity.

The analyses revealed an important distinction:

- **gene-level clustering** tended to retain stronger tumor-of-origin structure;
- **pathway-level clustering** more frequently grouped subtypes from different tumor types, indicating functional convergence across cancers.

## Three recurring patterns of molecular subtype organization

The integrated analyses identified three recurring patterns.

### 1. Gene- and pathway-concordant subtype organization

Some subtype relationships were preserved at both the gene and pathway levels, indicating concordant molecular organization across analytical scales.

### 2. Pathway-convergent but sub-pathway-divergent organization

Some subtypes shared enrichment of the same HALLMARK pathway but involved different sets of genes within different internal network modules.

For example, **BRCA.Her2** and **PCPG.Corticaladmixture** both showed involvement of `HALLMARK_E2F_TARGETS`, while their subtype-associated genes preferentially mapped to different sub-pathways.

### 3. Cross-tumor pathway convergence

Subtypes originating from different tumor types sometimes clustered together at the pathway level despite distinct tissue origins and gene-level signatures.

Representative examples included:

- **ACC.CIMP.high**, **COAD_GI.GS**, and **LIHC.iCluster:2**, with shared cell-cycle/proliferative programs;
- **AML.3**, **AML.4**, **OVCA.Immunoreactive**, and **THCA.2**, with shared immune/interferon-related programs.

These observations support a model in which tissue-associated gene programs can coexist with shared functional states across cancers.

## HALLMARK sub-pathway construction with OmicPath

To investigate pathway structure beyond conventional gene-set enrichment, we used our R package **OmicPath**:

**OmicPath: an R package for gene set analysis and pathway network analysis**  
https://github.com/CBIIT-CGBB/OmicPath

OmicPath supports gene-set analysis, pathway-to-network mapping, extraction of gene relationships, and graph-based pathway/network analyses.

### Installation

```r
library(devtools)
install_github("CBIIT-CGBB/OmicPath")
```

For the present study, HALLMARK gene sets were represented as gene-interaction networks and partitioned into smaller network modules ("sub-pathways"). Subtype-associated genes were then mapped back to these modules to test whether genes associated with different tumor subtypes were distributed uniformly or preferentially concentrated in specific sub-pathways.

This network-level analysis was used to distinguish:

> **shared pathway enrichment** from **shared molecular implementation of the pathway**.

In several examples, subtypes converged at the HALLMARK pathway level while remaining distinct at the sub-pathway level.

## Conceptual workflow

```text
TCGA RNA-seq + existing subtype annotations
                    |
                    v
        Machine-learning prediction
     glmnet / knn / rf / svmLinear / xgbTree
                    |
                    v
     Performance-Weighted Subtype Voting
                    |
                    v
       Extended TCGA subtype resource
                    |
                    v
       Subtype-associated gene analysis
                    |
                    v
         HALLMARK pathway enrichment
                    |
            +-------+-------+
            |               |
            v               v
     Gene-based        Pathway-based
      clustering         clustering
            |               |
            +-------+-------+
                    |
                    v
        Pan-cancer functional convergence
                    |
                    v
      HALLMARK pathway network analysis
                    |
                    v
          Sub-pathway decomposition
                    |
                    v
   Shared pathways vs. distinct network modules
```

## Data sources

The study uses data from [**The Cancer Genome Atlas (TCGA)**](https://portal.gdc.cancer.gov/), including RNA-seq expression data and existing molecular subtype annotations.

Subtype annotations were assembled using the Bioconductor package [**TCGAbiolinks**](https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) and were extended in this study for samples lacking existing subtype assignments.

Users should consult the original TCGA/GDC resources and subtype-defining publications when interpreting individual molecular subtype labels.

## Intended use

The subtype files in this repository are provided to support research and reproducible downstream analyses.

They may be useful for:

- pan-cancer molecular studies;
- subtype-specific biomarker analyses;
- pathway and network analyses;
- machine-learning benchmarking;
- integration with TCGA clinical, genomic, epigenomic, or other molecular data.

These subtype predictions are intended for **research use** and should not be interpreted as clinical diagnostic classifications.

## Reproducibility

For analyses involving pathway-network construction and sub-pathway analysis, see:

**CBIIT-CGBB/OmicPath**  
https://github.com/CBIIT-CGBB/OmicPath

The TCGA subtype assignments associated with this study are available in:

**`TCGA_Subtype/`**  
https://github.com/CBIIT-CGBB/TCGA-Subtypes/tree/main/TCGA_Subtype

## Citation

If you use the subtype resource or methods from this repository, please cite the associated manuscript once published:

> **Integrative Machine Learning and Network Analysis Extend TCGA Subtypes and Reveal Pan-Cancer Functional Convergence.**  
> Manuscript in preparation.

Please also cite the original TCGA/GDC datasets and relevant subtype-defining publications as appropriate.

## Contact and issues

Questions, reproducibility issues, or suggestions can be submitted through the GitHub issue tracker:

https://github.com/CBIIT-CGBB/TCGA-Subtypes/issues

