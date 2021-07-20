# Cross-species RNA-Seq analysis

Human and mouse RNA-Seq/scRNA-Seq downstream comparison, considering baseline expresison level? ~~or homology strucutre?~~

Eventually will come to 2D pathway enrichment analysis described in [Cox and Mann, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3489530/).


## DE and pathway annotation (regular way)
- 1.getmat.R
  - RNA-seq quantification result --> matrix
- 2.setcutoff.R
  - Kernel density estimation approach to select the normalized count cutoffs for low-expressed genes on a per-species basis ([Ramsey, 2018](http://link.springer.com/10.1007/978-1-4939-7456-6_14))
  - _haven't figure out if it's necessary._
- 3.biomart.R
  - To map mouse Ensembl ID to human HGNC symbol using `biomaRt`.
- 4.DE.R
  - Using `DESeq2` to get DE result.
  - 4-x.mouseDE.R
    - Skipped cutoff step.
- 5-1.fgsea.R
  - GSEA enrichment by `fgsea`.
- 5-2.gsea_pathwayfc.R
  - Get mean of gene log2FC of certain pathway, defined it as "pathway FC".
  - _Turns out it's useless._
- 6-1.topgo_preldefinedist.R
  - 6-1-x.topgo-genescore.R
  - ORA GO enrichment by `topGO`.
  - Two way to build genelist object, same results.
- 6-2.go_pathwayfc.R
  - Get mean of gene log2FC of certain pathway, defined it as "pathway FC".
  -  _Turns out it's useless as ever._
- 6-x
  - Just different ways of annotation
  - useless
- 7.xgsa.R
  - XSGA method described in [Djordjevic, Kusumi and Ho, 2016](https://academic.oup.com/bioinformatics/article/32/17/i620/2450758)
  - XSGA takes homology structure into account, but doesn't work for this analysis.
- 8-1.overlap_gsea.R
  - Use NES from `gsea` result as axes
  - _But 2D enrichment doesn't work like this, it'a failed attempt._
- 8-2.overlap_gofc.R
- 8-3.overlap_gseafc.R
  - 8-2 and 8-3 use pathway FC as axes
  - _still, useless as ever._
- 9-1.position_score_pathfc.R
  - Calculate postion score descirbed in [Cox and Mann, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3489530/) by ranks of pathway FC from 5-2 and 6-2.
  - Pairs of positon score as axes.
  - _MESS_
- 9-2.position_score_genefc.R
  - Calculate positiob score by ranks of gene log2FC.
  - This position score is actually the wanted one, but except this, other value between two species/sample is still not comparable.
- 10-1.2D_pathfc_ps.R
  - Prove 9-1 is rubbish.
- 10-2.2D_genefc_ps.R
  - Prove 9-2 is okay but still wrong.

## 2D annotation
- 1.diapause_geneset.R
  - Generate geneset objects.
- 2.map2path.R
  - Map every gene in DE table to all genesets.
  - Massive matrix.

    |        | Path 1 | Path 2 | Path 3 | Path 4 | ... |
    | :----: | :----: | :----: | :----: | :----: | :-: |
    | gene A |    1   |    0   |    0   |    0   | ... |
    | gene B |    0   |    1   |    0   |    0   | ... |
    | gene C |    0   |    0   |    1   |    0   | ... |
    |   ...  |   ...  |   ...  |   ...  |   ...  | ... |
- 3.rank-MANOVA.R
  - Core of this 2D enrichment method.
  - MANOVA test on gene rank.
  - Genes are ranked by log2FC (Wald statistics also works).
  - MANOVA to get p-val for each pathway, p-vals then were adjusted by BH method.
  - MANOVA input table looks like:
    
    | Genes | Rank.human | Rank.mouse |   Group  |
    | :---: | :--------: | :--------: | :------: |
    |   A   |      2     |      1     |  in path |
    |   B   |      1     |      2     | the rest |
    |   C   |      3     |      3     | the rest |
 
- 4.position_score.R
  - Postion score for each pathway, calculated by gene ranks.
- 5.2Ddotplot.R
  - Final plot.
  - Pairs of positon score as axes.
- 6*_dotplot.R
  - Dotplots of all genes for each pathway.
  - To look at how genes in a certain pathway behave between two species/samples.
- mitch_*.R
  - Use Bioconductor package [`mitch`](https://doi.org/10.1186/s12864-020-06856-9) to do 2D annotation.
  - By default, it takes `stat` column form `DESeq2` result, and rank genes from the midpoint.
  - When adopted with log2FC, the results are exactly the same, except `mitch` is way much faster than the manual appraoch.
