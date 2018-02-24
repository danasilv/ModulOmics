# ModulOmics

## Summary
ModulOmics is a method to identify cancer driver pathways de novo by integrating multiple data types (protein-protein interactions, mutual exclusivity of mutations or copy number alterations, transcriptional co-regulation, and RNA co-expression) into a single probabilistic model. ModulOmics uses patient data, and protein-protein interaction networks, based here on the Hippie database and regulatory connections, based here on the Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining database (TRRUST).

## Installation
ModulOmics is available as ```R``` code. It has been tested on Mac and Linux, on R 3.4.3. It needs an active ```cplex``` (IBM CPLEX Optimizer) instance installed, and depends on the following R packages:
* ```TiMEx```, current version downloadable from [its github repository](https://githhub.com/csimona/TiMEx)
* ```igraph``` v. 1.1.2
* ```cplexAPI``` v. 1.3.3
* ```gtools``` v. 3.5.0

For more details on the input and output structures of ModulOmics, please see the pdf manual.


