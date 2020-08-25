# OptEnrichedSetCover.jl package

*OptEnrichedSetCover.jl* stands for "Optimal Enriched-Set Cover".
It implements the advanced *gene set enrichment analysis* (GSEA) that could be applied to
study complex experimental designs.

## Introduction

[Gene set enrichment analysis](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis) (GSEA) is a popular bioinformatic
approach for the system-level analysis of high-throughput biological data (next-gen sequencing,
proteomics etc). The current knowledge of cellular biology is formalized by *annotating*
the biological *entities* (gene, proteins etc) according to the trait of interest (biological role,
cellular localization, presence of sequence motifs etc). The aim of GSEA is then to detect, which
annotations are particularly enriched in the outcome of the experiment.

Detailed and precise annotations (e.g. [Gene Ontologies](http://geneontology.org/))
may contain thousands of *annotation terms* and be organized into hierarchies, where less specific
*annotation terms* are subdivided into more specific.
The downside of high detalization is the increased complexity of GSEA results: more terms are reported
as enriched, while many of them are *redundant* (a significant portion of genes is shared between
the terms that describe the same feature or process, just at slightly different levels of detail,
e.g. [*"Ribosome"*](https://www.ebi.ac.uk/QuickGO/term/GO:0005840) and
[*"Cytosolic Ribosome"*](https://www.ebi.ac.uk/QuickGO/term/GO:0022626)).
There are different approaches to simplify the interpretation of GSEA output in such cases:
  * using simplified annotation collections ([GOslim](http://geneontology.org/docs/go-subset-guide/))
  * automatic detection of the optimal level of detail ([topGO]())
  * grouping similar annotations into clusters by similarity ([GOzilla](), [GOSemSim]())
Often, these methods work only for the limited range of annotation collections, where the relationship between
different annotation terms is explicitly specified (topGO) or rely on the pregenerated simplified annotations.

As the technology advances and gets more accessible, data sets become more and more complex, the number of
measured experimental conditions grows and could reach 50-100 or even more. GSEA could help to grasp
the high-level differences and similarities between the conditions by highlighting condition-specific or
shared annotation traits. However, complex data multiply the redundancy problem of GSEA results as each condition
introduces its own set of redundant terms. Some redundancy-reduction approaches are simply not designed to
work with complex experimental designs (e.g. similarity graph visualization of enriched terms).
The naive application of redundancy-reducing GSEA to one condition at a time does not guarantee that the
same level of detail would be chosen for different conditions, so the combined result may contain redundant terms
and be suboptimal for conditions comparison.

*OptEnrichedSetCover* method was designed to address these challenges:
  * work with any annotation collection, including those with redundancy, but no explicit terms hierarchy
    (e.g. [CORUM protein complexes](http://mips.helmholtz-muenchen.de/corum/));
  * allow the analysis of complex experimental designs and deduction of annotation terms that best
    describe all variation and commonalities that are present in the data.

The mathematical details are provided in the ["Method description"](@ref method) section.
