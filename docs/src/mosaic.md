# [Sets collection](@id mosaic)

## Sets mosaic

Gene Ontologies provide thousands of annotation terms.
However, different genes often share the same annotations. This observation could be used to improve the
processing of large collections. The genes could be grouped into disjoint *"tiles"* so that:
* all genes of the same tile share the same annotations
* the tiles don't overlap and their *mosaic* covers all the genes 
* any annotation term could be represented as the (disjoint) union of the tiles

This *"mosaic"* representation of the sets collection is implemented by the `SetMosaic`
type. Internally it uses [`SparseMaskMatrix`](@ref OptEnrichedSetCover.SparseMaskMatrix) to maintain efficient
element-to-set, element-to-tile etc mappings.
Since there are (much) fewer tiles than individual elements, operations like sets union or intersection
could be made much faster.

```@docs
SetMosaic
nelements
nsets
ntiles
set
setsize
tile
```

## [Set relevance](@id set_relevance)

In the normal experiments, even the high-throughput ones, it's not possible to detect all annotated
entities (genes or proteins). There's e.g. *detection bias* due to the experimental protocol or
the *sensitivity limit* of the instrument.
This must be taken into account when estimating, whether a particular gene set is enriched among e.g.
upregulated genes -- only the detected genes should be considered for the enrichment scores.

It may happen that the two distinct annotation terms share the same set of observed genes.
In that case, their enrichment scores would be identical.
If the enrichment is significant, both terms could be included in the report,
but that would increase its *redundancy*.
To solve this issue, *OptEnrichedSetCover* introduces the *set relevance score*.
For example, if both *"ribosome"* and *"small ribosomal subunits"* terms are significant,
but the genes of the *large ribosomal subunit* are not detected in the data, it's natural
to prefer the "small ribosomal subunit" term over the whole ribosome.
The *relevance score* formalizes that by estimating the enrichment of *detected* entities
within each annotation term. See [cover quality](@ref cover_quality_enrichment) section for the exact
definition of the *relevance score*.

```@docs
set_relevance
```

## Masked sets mosaic

```@docs
MaskedSetMosaic
mask
unmask
```