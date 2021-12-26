var documenterSearchIndex = {"docs":
[{"location":"sparse_mask_matrix.html#sparse_mask_matrix","page":"Sparse mask matrix","title":"Sparse mask matrix","text":"","category":"section"},{"location":"sparse_mask_matrix.html","page":"Sparse mask matrix","title":"Sparse mask matrix","text":"Modules=[OptEnrichedSetCover]\nPages = [\"sparse_mask_matrix.jl\"]","category":"page"},{"location":"sparse_mask_matrix.html#OptEnrichedSetCover.SparseMaskMatrix","page":"Sparse mask matrix","title":"OptEnrichedSetCover.SparseMaskMatrix","text":"Sparse representation of Matrix{Bool}, when all \"falses\" are \"structural\". It uses compressed sparse column (CSC) representation, except no non-zero values need to be stored.\n\n\n\n\n\n","category":"type"},{"location":"sparse_mask_matrix.html#OptEnrichedSetCover.SparseMaskMatrix-Tuple{AbstractArray{Bool,2}}","page":"Sparse mask matrix","title":"OptEnrichedSetCover.SparseMaskMatrix","text":"SparseMaskMatrix(mtx::AbstractMatrix{Bool}) -> SparseMaskMatrix\n\nSparse representation of boolean matrix.\n\n\n\n\n\n","category":"method"},{"location":"sparse_mask_matrix.html#OptEnrichedSetCover.SparseMaskMatrix-Tuple{Integer,AbstractArray{Array{Int64,1},1}}","page":"Sparse mask matrix","title":"OptEnrichedSetCover.SparseMaskMatrix","text":"SparseMaskMatrix(m::Integer, rowvals_percol::Vector{Vector{Int}}) -> SparseMaskMatrix\n\nConstruct SparseMaskMatrix given the vector of true row indices per each column and the total number of rows (m).\n\n\n\n\n\n","category":"method"},{"location":"sparse_mask_matrix.html#OptEnrichedSetCover.SparseMaskMatrix-Union{Tuple{T}, Tuple{Any,AbstractDict{T,Int64}}} where T","page":"Sparse mask matrix","title":"OptEnrichedSetCover.SparseMaskMatrix","text":"SparseMaskMatrix(sets, elm2ix::Dict{T, Int}) where T -> SparseMaskMatrix\n\nConstruct SparseMaskMatrix from the collection of sets.\n\nArguments\n\nsets: a collection of sets, one set per mask matrix column\nelm2ix: maps set element to its index (mask row index)\n\n\n\n\n\n","category":"method"},{"location":"method.html#method","page":"Method description","title":"Method description","text":"","category":"section"},{"location":"method.html#Basic-Definitions","page":"Method description","title":"Basic Definitions","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"Let","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"G =  g_1 g_2 ldots g_N : the annotated entities (e.g. ENTREZ genes)\nO =  o_1 o_2 ldots o_M : the entities observed in the experiments (e.g. protein groups)\nmathrmobs G mapsto O: the function that maps annotated entities to the observed ones; forall o in O exists mathrmobs^-1(o) subset G (an observed entity may refer to several annotated ones)\nmathcalA =  A_1 A_2 ldots A_N_T : the collection of annotation terms, A_i subset G\nmathcalX =  X_1 X_2 ldots X_M_X : the collection of experiment hits (e.g. significantly regulated proteins), X_j subset O\nC =  c_1 c_2 ldots c_N_T , 0 leq c_i leq 1: weighted cover of the experiment hits, where c_i defines the probability to use the annotation term A_i for the cover (if c_i = 0, the term A_i is never used; if c_i = 1, it is always in the cover)\nP_mathrmFET(A H mathrmAll) = P(X geq N_Acap H), where X propto mathrmHypergeomtetric(N_A N_H N_mathrmAll) – the P-value of the Fisher's Exact Test that the overlap of the sets A subset mathrmAll (e.g. annotation term) and H subset mathrmAll (e.g. hits) is significant","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"<img alt=\"Figure: Basic definitions\" src=\"assets/method_defs.svg\" width=\"80%\"/>","category":"page"},{"location":"method.html#cover_quality","page":"Method description","title":"Cover quality","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"The method works by defining the cover quality q(C mathcalX) of how well the selected terms (C) cover the experiment hits (mathcalX). Since some aspects of cover quality could not be simultaneously satisfied (e.g. using as few annotation terms as possible vs covering all hits), the score q is composed of several components, and multi-objective optimization is employed to find the family of optimal covers.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"More specifically, q = (q_e q_r q_u q_c), where","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_e: enrichment score, measures how much the cover terms are enriched within the experiment hits (across all experiments)\nq_r: redundancy score, penalizes for the presence of redundant terms in the cover\nq_u: hits cover score, measures how well individual hits are covered by C\nq_c: non-hits cover score, penalizes for covering observed entities that are not hits","category":"page"},{"location":"method.html#cover_quality_enrichment","page":"Method description","title":"Enrichment score","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_e(C mathcalX) = sum_i=1^N_T c_i L(A_i O) E(mathrmobs(A_i) mathcalX)","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"where L(A_i O) is the relevance of the term A_i with respect to the observable elements O, and E(mathrmobs(A_i) mathcalX) is the overall enrichment of A_i in the experiment hits mathcalX.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"Relevance component allows to prioritize the annotation terms that are better represented among the observed entities, but otherwise have the same enrichment as less represented ones:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"L(A O) = max(1 - P_mathrmFET(A mathrmobs^-1(O) G) L_min)^beta_R","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"The relevance score quickly approaches 1, when the representation of the term among all observed genes is significant. L_min in 0 1 controls that the terms don't get overpenalized due to low representation. For more practical discussion, see \"Set relevance\" section.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"The enrichment component is defined as","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"E(Y mathcalX) = epsilon - sum_j^M_X alpha^(j-1) left(- log P_mathrmFET(Y X_i_j O)right)^beta","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"where epsilon ge 0 is the penalty for including arbitrary term into a cover; i_1, i_2, ldots, i_M_X is the permutation of 1, 2, ldots, M_X that orders the hits mathcalX by the significance of their overlap with Y (from most to least significant); alpha in (0 1 is the enrichment discount, and beta  0 is the enrichment shape.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"The discount parameter alpha is important for complex data (M_X gg 1). It allows terms that are very significant, but only in a few conditions, be preferred over the terms that are moderately represented in most conditions. The lower is alpha, the more enrichment score favors condition-specific terms.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"The enrichment shape parameter beta affects the size of the selected annotation terms. For small beta the method would prefer several smaller and more specific terms over the single one that combines multiple small terms.","category":"page"},{"location":"method.html#Redundancy-score","page":"Method description","title":"Redundancy score","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_r(C) = sum_i=1^N_T sum_j=1^N_T min(c_i c_j) R(mathrmobs(A_i) mathrmobs(A_j))","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"where R(mathrmobs(A_i) mathrmobs(A_j)) measures the redundancy of observed entities from A_i and A_j annotation terms as the enrichment of one term in their union:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"R(B_i B_j) = -left( -log min left(\n  P_mathrmFET(B_i B_i cup B_j N_O (1 + Delta))\n  P_mathrmFET(B_j B_i cup B_j N_O (1 + Delta)) right)right)^beta","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"where Delta N_O ge 0 is added to keep the overlap significant for very large annotation terms that are comparable in size to the set of all observed genes O.","category":"page"},{"location":"method.html#Hits-cover-scores","page":"Method description","title":"Hits cover scores","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_u counts how many hit entities were not covered by the selected terms:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_u(C mathcalX) = sum_o in cup_i X_i left(1 - max_j in С(o) c_jright)","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"where С(o) is the set of all indices j, s.t. o in mathrmobs(A_j).","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"Similarly, q_c counts how many non-hit observations in each experiment are covered by the selected terms:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_c(C mathcalX) = sum_i=1^M_X sum_o in Osetminus X_i max_j in С(o) c_j","category":"page"},{"location":"method.html#Optimal-covers","page":"Method description","title":"Optimal covers","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"With q = (q_e q_u q_c q_r) cover quality defined, one can consider the multi-objective optimization problem of finding the Pareto-optimal family of covers mathcalC_* =  C_* 1 C_* 2 ldots :","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"forall C notin mathcalC_* exists C_* in mathcalC_* \nq_e(C_*) leq q_e(C) q_r(C_*) leq q_r(C) q_u(C_*) leq q_u(C_*) q_c(C_*) leq q_c(C)","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"and at least for one of the q components the inequality is strict.","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"The resulting Pareto front mathcalC_* will contain: the covers with minimally redundant but less enriched terms, very redundant covers of more enriched terms, as well as intermediate solutions between these extremes. It's possible to weight the importance of individual quality components (w_e = 1, w_r  0, w_u ge 0, w_c ge 0) to identify the unique optimal cover C_* Sigma in mathcalC_* that minimizes the components weighted sum:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"q_Sigma(C) = q_e(C) + w_r q_r(C) + w_u q_u(C) + w_c q_c(C) to min quad C in mathcalC_*","category":"page"},{"location":"method.html#cover_score_convolution","page":"Method description","title":"Components convolution","text":"","category":"section"},{"location":"method.html","page":"Method description","title":"Method description","text":"The number of points on the approximated 4-objective Pareto front is O(frac1delta^(4-1)), where delta  0 is the grid step in the multi-objective space. The maintenance of such front creates significant overhead during optimization. However, as the components q_u and q_c correlate with the enrichment component q_e, they could be convoluted with the q_e:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"tildeq_e = q_e + w_u q_u + w_c q_c","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"or omitted altogether by setting w_u = w_c = 0. The optimization of the simplified 2-objective (tildeq_e q_r) problem would require much less resources. The corresponding Pareto front would essentially show the trade-offs between the annotation terms enrichment tildeq_e and their redundancy q_r. One can define the redunancy index k_r as","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"k_r(C) = fracq_r(C)-tildeq_e(C)","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"High k_r(C) means the terms in the cover C are redundant. In practice, covers with k_r(C) gg 1, although marginally improving q_e, do not reveal any new patterns in the data. To improve the optimization efficiency, it's possible to use the modified cover score tildeq, so that the covers with high k_r become dominated by less redundant covers, and thus keep only the relevant solutions on the Pareto front:","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"tildeq = left(tildeq_e + w_r t q_r (1 - t) q_r right)textwhere\nt = frac1k_r left(w_r + expleft(s^-1 - sright)right)^-1 quad\ns = alpha_k max(k_r - k_max 0)","category":"page"},{"location":"method.html","page":"Method description","title":"Method description","text":"k_max is the maximal accepted redundancy index (by default the method uses k_max = 1), and alpha_k  0 defines the transformation strength. With this definition, when k_r(C) leq k_max, then s = 0, t = 0, and both score components are unchanged: tildeq = (tildeq_e q_r). When k_r(C) gg k_max, i.e. when the redundancy of the cover is very high, t approx frac1k_r w_r = frac-tildeq_ew_r q_r, and tildeq approx (0 q_r + fractildeq_ew_r). The latter solution is dominated by any less redundant cover.  Also note that, for any C, tildeq_1 + w_r tildeq_2 = q_Sigma, so this  transformation preserves the optimal cover C_*Sigma that minimizes q_Sigma(C).","category":"page"},{"location":"mosaic.html#mosaic","page":"Sets collection","title":"Sets collection","text":"","category":"section"},{"location":"mosaic.html#Sets-mosaic","page":"Sets collection","title":"Sets mosaic","text":"","category":"section"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"Gene Ontologies provide thousands of annotation terms. However, different genes often share the same annotations. This observation could be used to improve the processing of large collections. The genes could be grouped into disjoint \"tiles\" so that:","category":"page"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"all genes of the same tile share the same annotations\nthe tiles don't overlap and their mosaic covers all the genes \nany annotation term could be represented as the (disjoint) union of the tiles","category":"page"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"This \"mosaic\" representation of the sets collection is implemented by the SetMosaic type. Internally it uses SparseMaskMatrix to maintain efficient element-to-set, element-to-tile etc mappings. Since there are (much) fewer tiles than individual elements, operations like sets union or intersection could be made much faster.","category":"page"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"SetMosaic\nnelements\nnsets\nntiles\nset\nsetsize\ntile","category":"page"},{"location":"mosaic.html#OptEnrichedSetCover.SetMosaic","page":"Sets collection","title":"OptEnrichedSetCover.SetMosaic","text":"SetMosaic{T,S}\n\nRepresents a collection of (potentially overlapping) sets as a \"mosaic\" of non-overlapping \"tiles\".\n\nType parameters\n\nT: type of set elements\nS: type of set keys\n\n\n\n\n\n","category":"type"},{"location":"mosaic.html#OptEnrichedSetCover.nelements","page":"Sets collection","title":"OptEnrichedSetCover.nelements","text":"nelements(mosaic::SetMosaic) -> Int\n\nThe number of distinct elements (i.e. genes) in the sets collection.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#OptEnrichedSetCover.nsets","page":"Sets collection","title":"OptEnrichedSetCover.nsets","text":"nsets(mosaic::SetMosaic) -> Int\n\nThe number of sets in the collection.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#OptEnrichedSetCover.ntiles","page":"Sets collection","title":"OptEnrichedSetCover.ntiles","text":"nsets(mosaic::SetMosaic) -> Int\n\nThe number of tiles (pairwise disjoint sets) in the mosaic representation of the set collection.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#OptEnrichedSetCover.set","page":"Sets collection","title":"OptEnrichedSetCover.set","text":"set(mosaic::SetMosaic, i::Integer) -> AbstractVector{Int}\n\nGet the i-th set as the vector of indices of its elements.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#OptEnrichedSetCover.setsize","page":"Sets collection","title":"OptEnrichedSetCover.setsize","text":"set(mosaic::SetMosaic, i::Integer) -> Int\n\nGet the size of the i-th set.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#OptEnrichedSetCover.tile","page":"Sets collection","title":"OptEnrichedSetCover.tile","text":"tile(mosaic::SetMosaic, i::Integer) -> AbstractVector{Int}\n\nGet the i-th mosaic tile as the vector of indices of its elements.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#set_relevance","page":"Sets collection","title":"Set relevance","text":"","category":"section"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"In the normal experiments, even the high-throughput ones, it's not possible to detect all annotated entities (genes or proteins). There's e.g. detection bias due to the experimental protocol or the sensitivity limit of the instrument. This must be taken into account when estimating, whether a particular gene set is enriched among e.g. upregulated genes – only the detected genes should be considered for the enrichment scores.","category":"page"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"It may happen that the two distinct annotation terms share the same set of observed genes. In that case, their enrichment scores would be identical. If the enrichment is significant, both terms could be included in the report, but that would increase its redundancy. To solve this issue, OptEnrichedSetCover introduces the set relevance score. For example, if both \"ribosome\" and \"small ribosomal subunits\" terms are significant, but the genes of the large ribosomal subunit are not detected in the data, it's natural to prefer the \"small ribosomal subunit\" term over the whole ribosome. The relevance score formalizes that by estimating the enrichment of detected entities within each annotation term. See cover quality section for the exact definition of the relevance score.","category":"page"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"set_relevance","category":"page"},{"location":"mosaic.html#OptEnrichedSetCover.set_relevance","page":"Sets collection","title":"OptEnrichedSetCover.set_relevance","text":"set_relevance(nset_observed::Integer, nset::Integer,\n              nobserved::Integer, ntotal::Integer) -> Float64\n\nCalculates the relevance weight of the set that contains nset elements, nset_observed of which were present (not necessarily enriched) in the data that identified nobserved elements out of all known (ntotal). It is used by SetMosaic to penalize the sets, which could not be observed in the data (e.g. biological processes or pathways that involve proteins not expressed by the cells used in the experiments).\n\nWhile for MaskedSetMosaic it's recommended to use the IDs of data entities (e.g. protein group IDs for proteomic data) to correctly count the set sizes and estimate enrichment; set_relevance() should use the counts derived from the original IDs of the annotation database (e.g. UniProt accession codes). Otherwise it's not possible to correctly estimate the number of elements that belong to the given annotated set, but were not observed in the data.\n\nThe returned value is the probability that no more than nset_observed elements were observed at random.\n\n\n\n\n\n","category":"function"},{"location":"mosaic.html#Masked-sets-mosaic","page":"Sets collection","title":"Masked sets mosaic","text":"","category":"section"},{"location":"mosaic.html","page":"Sets collection","title":"Sets collection","text":"MaskedSetMosaic\nmask\nunmask","category":"page"},{"location":"mosaic.html#OptEnrichedSetCover.mask","page":"Sets collection","title":"OptEnrichedSetCover.mask","text":"mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};\n     [experiment_ids::Union{AbstractVector, AbstractSet, Nothing} = nothing],\n     [min_nmasked=1], [max_setsize=nothing],\n     [max_overlap_logpvalue=0.0]) -> MaskedSetMosaic\n\nConstruct MaskedSetMosaic from the SetMosaic and the collection of element masks.\n\nArguments\n\nmin_nmasked: the minimal number of masked elements in a set to include in the mosaic\nmax_setsize (optional): ignore the annotation sets bigger than the specified size\nmax_overlap_logpvalue: the threshold of Fisher's Exact Test log P-value of the overlap  between the set and the mask for the inclusion of the set into the mosaic.  0 accepts all sets.\n\n\n\n\n\n","category":"function"},{"location":"cover_problem.html#cover_problem","page":"Cover problem","title":"Cover Problem","text":"","category":"section"},{"location":"cover_problem.html","page":"Cover problem","title":"Cover problem","text":"OptEnrichedSetCover.CoverParams\nAbstractCoverProblem\nMultiobjCoverProblem\nOptEnrichedSetCover.MultiobjOptimizerParams\noptimize\nOptEnrichedSetCover.MultiobjProblemSoftFold2d\nscore\nOptEnrichedSetCover.MultiobjCoverProblemResult","category":"page"},{"location":"cover_problem.html#OptEnrichedSetCover.CoverParams","page":"Cover problem","title":"OptEnrichedSetCover.CoverParams","text":"CoverParams(; [sel_tax=0.0], [set_shape=1.0],\n              [min_weight=1E-2], [mask_discount=0.9], [setXset_factor=1.0],\n              [uncovered_factor=0.1], [covered_factor=0.001],\n              [set_relevance_shape=0.5], [set_relevance_min=0.5]) -> CoverParams\n\nSpecify parameters for the AbstractCoverProblem. One can specify a non-default value for a particular parameter.\n\nSee cover quality score for the detailed description of the parameters.\n\nArguments\n\nsel_tax: epsilon, the constant added to the set score of each selected set to penalize  the number of sets in the cover\nset_shape: beta, applied to set or set×set scores\nmin_weight: minimal non-zero set probability ???\nmask_discount: alpha, how much the overlap score of each subsequent mask (from most to less enriched) is discounted\nsetXset_factor: w_r, the weight of redundancy cover quality component (setXset_score scale), 0 = no redunancy penalty\nuncovered_factor: w_u, the weight of uncovered hits component in cover quality\ncovered_factor: w_c, the weight of covered non-hits component in cover quality\nset_relevance_shape: beta_L, how much set relevance affects set score, 0 = no effect\nset_relevance_min: L_min, if shaped relevance is below, it's set to set_relevance_min\n\n\n\n\n\n","category":"type"},{"location":"cover_problem.html#OptEnrichedSetCover.AbstractCoverProblem","page":"Cover problem","title":"OptEnrichedSetCover.AbstractCoverProblem","text":"AbstractCoverProblem{T}\n\nOptimal Enriched-Set Cover problem – choose the sets from the collection 𝒞 to cover the masked(selected) elements M. The optimal sets cover C = {c₁, c₂, ..., cₙ} ⊂ 𝒞 has to deliver 3 goals:\n\nbe relevant (i.e. minimize the P-values of M and cᵢ sets overlap)\nbe minimal (i.e. minimize the number of sets in C)\nbe non-redundant (i.e. minimize the P-values of the pairwise non-overlap of C sets with each other).\n\nFuzzy set selection is possible – each set is assigned a weight from [0, 1] range.\n\n\n\n\n\n","category":"type"},{"location":"cover_problem.html#OptEnrichedSetCover.MultiobjCoverProblem","page":"Cover problem","title":"OptEnrichedSetCover.MultiobjCoverProblem","text":"Abstract multi-objective optimal Enriched-Set Cover problem.\n\nSee \"Method Description\" for more details.\n\n\n\n\n\n","category":"type"},{"location":"cover_problem.html#OptEnrichedSetCover.MultiobjOptimizerParams","page":"Cover problem","title":"OptEnrichedSetCover.MultiobjOptimizerParams","text":"Parameters for the Borg-based optimization of MultiobjCoverProblem.\n\nSee optimize.\n\n\n\n\n\n","category":"type"},{"location":"cover_problem.html#OptEnrichedSetCover.optimize","page":"Cover problem","title":"OptEnrichedSetCover.optimize","text":"optimize(problem::MultiobjCoverProblem,\n         [opt_params::MultiobjOptimizerParams]) -> MultiobjCoverProblemResult\n\nOptimize MultiobjCoverProblem and return the result. Uses Borf multi-objective optimization method from BlackBoxOptim.jl package.\n\n\n\n\n\n","category":"function"},{"location":"cover_problem.html#OptEnrichedSetCover.MultiobjProblemSoftFold2d","page":"Cover problem","title":"OptEnrichedSetCover.MultiobjProblemSoftFold2d","text":"Transforms the 4-component cover quality score into 2-component score that makes the highly redundant solutions dominated by any less redundant ones.\n\nSee \"Cover score convolution\" section for the discussion.\n\nArguments\n\nsetXset_factor: w_r, same as in CoverParams\nuncovered_factor: w_u, same as in CoverParams\ncovered_factor: w_c, same as in CoverParams\nratio_threshold: k_max, defaults to 1\nshape: alpha_k, defaults to 0.5.\n\n\n\n\n\n","category":"type"},{"location":"cover_problem.html#OptEnrichedSetCover.score","page":"Cover problem","title":"OptEnrichedSetCover.score","text":"score(w::AbstractVector{Float64}, problem) -> NTuple{4, Float64}\n\nUnfolded multiobjective score (fitness) of the OESC coverage.\n\nw: probabilities of the sets being covered\n\nSee \"Cover quality\".\n\n\n\n\n\n","category":"function"},{"location":"cover_problem.html#OptEnrichedSetCover.MultiobjCoverProblemResult","page":"Cover problem","title":"OptEnrichedSetCover.MultiobjCoverProblemResult","text":"The result of optimize. Contains the solutions on the Pareto front: weights of the annotation terms and corresponding cover scores.\n\n\n\n\n\n","category":"type"},{"location":"index.html#OptEnrichedSetCover.jl-package","page":"Introduction","title":"OptEnrichedSetCover.jl package","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"OptEnrichedSetCover.jl stands for \"Optimal Enriched-Set Cover\". It implements the advanced gene set enrichment analysis (GSEA) that could be applied to study complex experimental designs.","category":"page"},{"location":"index.html#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Gene set enrichment analysis (GSEA) is a popular bioinformatic approach for the system-level analysis of high-throughput biological data (next-gen sequencing, proteomics etc). The current knowledge of cellular biology is formalized by annotating the biological entities (gene, proteins etc) according to the trait of interest (biological role, cellular localization, presence of sequence motifs etc). The aim of GSEA is then to detect, which annotations are particularly enriched in the outcome of the experiment.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Detailed and precise annotations (e.g. Gene Ontologies) may contain thousands of annotation terms and be organized into hierarchies, where less specific annotation terms are subdivided into more specific. The downside of high detalization is the increased complexity of GSEA results: more terms are reported as enriched, while many of them are redundant (a significant portion of genes is shared between the terms that describe the same feature or process, just at slightly different levels of detail, e.g. \"Ribosome\" and \"Cytosolic Ribosome\"). There are different approaches to simplify the interpretation of GSEA output in such cases:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"using simplified annotation collections (GOslim)\nautomatic detection of the optimal level of detail (topGO)\ngrouping similar annotations into clusters by similarity (GOzilla, GOSemSim)","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Often, these methods work only for the limited range of annotation collections, where the relationship between different annotation terms is explicitly specified (topGO) or rely on the pregenerated simplified annotations.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"As the technology advances and gets more accessible, data sets become more and more complex, the number of measured experimental conditions grows and could reach 50-100 or even more. GSEA could help to grasp the high-level differences and similarities between the conditions by highlighting condition-specific or shared annotation traits. However, complex data multiply the redundancy problem of GSEA results as each condition introduces its own set of redundant terms. Some redundancy-reduction approaches are simply not designed to work with complex experimental designs (e.g. similarity graph visualization of enriched terms). The naive application of redundancy-reducing GSEA to one condition at a time does not guarantee that the same level of detail would be chosen for different conditions, so the combined result may contain redundant terms and be suboptimal for conditions comparison.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"OptEnrichedSetCover method was designed to address these challenges:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"work with any annotation collection, including those with redundancy, but no explicit terms hierarchy (e.g. CORUM protein complexes);\nallow the analysis of complex experimental designs and deduction of annotation terms that best describe all variation and commonalities that are present in the data.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The mathematical details are provided in the \"Method description\" section.","category":"page"},{"location":"array_pool.html#array_pool","page":"Array pool","title":"Array pool","text":"","category":"section"},{"location":"array_pool.html","page":"Array pool","title":"Array pool","text":"The ArrayPool type is meant to reduce the stress on Julia Garbage Collector by maintaining the pool of preallocated arrays that the program can temporarily borrow and then return back to the pool when not needed anymore.","category":"page"},{"location":"array_pool.html","page":"Array pool","title":"Array pool","text":"Modules=[OptEnrichedSetCover]\nPages = [\"array_pool.jl\"]","category":"page"},{"location":"array_pool.html#OptEnrichedSetCover.ArrayPool","page":"Array pool","title":"OptEnrichedSetCover.ArrayPool","text":"ArrayPool{T}\n\nHelps to maintain the pool of reusable arrays of different sizes and reduce the burden on garbage collection.\n\nType parameters\n\nT: the type of array elements\n\n\n\n\n\n","category":"type"},{"location":"array_pool.html#OptEnrichedSetCover.acquire!-Union{Tuple{T}, Tuple{OptEnrichedSetCover.ArrayPool{T},Integer}} where T","page":"Array pool","title":"OptEnrichedSetCover.acquire!","text":"acquire!(pool::ArrayPool{T}, size) where T -> Array{T}\n\nGets an array of a specific size from the pool. The acquired array should be returned back to the pool using release!. The size could be either an integer or a tuple of integers.\n\n\n\n\n\n","category":"method"},{"location":"array_pool.html#OptEnrichedSetCover.release!-Union{Tuple{T}, Tuple{OptEnrichedSetCover.ArrayPool{T},Array{T,N} where N}} where T","page":"Array pool","title":"OptEnrichedSetCover.release!","text":"release!(pool::ArrayPool{T}, arr::Array{T}) where T\n\nReleases the array previously obtained by acquire! back into the pool.\n\n\n\n\n\n","category":"method"}]
}
