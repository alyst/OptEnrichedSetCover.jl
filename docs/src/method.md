# [Method description](@id method)

## Basic Definitions

Let
  * ``G = \{ g_1, g_2, \ldots, g_N \}``: the *annotated entities* (e.g. ENTREZ genes)
  * ``O = \{ o_1, o_2, \ldots, o_M \}``: the entities *observed* in the experiments (e.g. *protein groups*)
  * ``\mathrm{obs}: G \mapsto O``: the function that *maps* annotated entities to the observed ones;
    ``\forall o \in O`` ``\exists \mathrm{obs}^{-1}(o) \subset G`` (an *observed* entity may refer to several *annotated* ones)
  * ``\mathcal{A} = \{ A_1, A_2, \ldots, A_{N_T} \}``: the collection of *annotation terms*, ``A_i \subset G``
  * ``\mathcal{X} = \{ X_1, X_2, \ldots, X_{M_X} \}``: the collection of *experiment hits*
    (e.g. significantly regulated proteins), ``X_j \subset O``
  * ``C = \{ c_1, c_2, \ldots, c_{N_T} \}``, ``0 \leq c_i \leq 1``: weighted *cover* of the experiment hits,
    where ``c_i`` defines the probability to use the annotation term ``A_i`` for the cover
    (if ``c_i = 0``, the term ``A_i`` is never used; if ``c_i = 1``, it is always in the cover)
  * ``P_{\mathrm{FET}}(A, H, \mathrm{All}) = P(X \geq N_{A\cap H})``, where
    ``X \propto \mathrm{Hypergeomtetric}(N_{A}, N_{H}, N_{\mathrm{All}})`` -- the ``P``-value of the
    [Fisher's Exact Test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) that the overlap
    of the sets ``A \subset \mathrm{All}`` (e.g. annotation term) and ``H \subset \mathrm{All}`` (e.g. hits) is significant

## [Cover quality](@id cover_quality)

The method works by defining the *cover quality* ``q(C, \mathcal{X})`` of how well the selected
terms (``C``) *cover* the experiment hits (``\mathcal{X}``). Since some aspects of cover
quality could not be simultaneously satisfied (e.g. using as few annotation terms as possible vs covering all hits),
the score ``q`` is composed of several components, and multi-objective optimization is employed to find the family of
optimal covers.

More specifically, ``q = (q_e, q_r, q_u, q_c)``, where
  * ``q_e``: *enrichment* score, measures how much the *cover terms* are enriched within the experiment hits (across all experiments)
  * ``q_r``: *redundancy* score, penalizes for the presence of redundant terms in the cover
  * ``q_u``: *hits cover* score, measures how well individual hits are covered by ``C``
  * ``q_c``: *non-hits cover* score, penalizes for covering observed entities that are not hits

### [Enrichment score](@id cover_quality_enrichment)

```math
q_e(C, \mathcal{X}) = \sum_{i=1}^{N_T} c_i L(A_i, O) E(\mathrm{obs}(A_i), \mathcal{X}),
```
where ``L(A_i, O)`` is the *relevance* of the term ``A_i`` with respect to the observable elements ``O``,
and ``E(\mathrm{obs}(A_i), \mathcal{X})`` is the *overall enrichment* of ``A_i`` in the experiment hits ``\mathcal{X}``.

*Relevance* component allows to prioritize the annotation terms that are better represented among the
observed entities, but otherwise have the same enrichment as less represented ones:
```math
L(A, O) = \max(1 - P_{\mathrm{FET}}(A, \mathrm{obs}^{-1}(O), G), L_{\min})^{\beta_R}.
```
The *relevance score* quickly approaches ``1``, when the representation of the term among all observed genes is significant.
``L_{\min} \in [0, 1]`` controls that the terms don't get overpenalized due to low representation.
For more practical discussion, see ["Set relevance"](@ref set_relevance) section.

The *enrichment* component is defined as
```math
E(Y, \mathcal{X}) = \epsilon - \sum_j^{M_X} \alpha^{(j-1)} \left(- \log P_{\mathrm{FET}}(Y, X_{i_j}, O)\right)^{\beta},
```
where ``\epsilon \ge 0`` is the penalty for including arbitrary term into a cover;
``i_1``, ``i_2``, ``\ldots``, ``i_{M_X}`` is the permutation of ``1``, ``2``, ``\ldots``, ``M_X``
that orders the hits ``\mathcal{X}``
by the significance of their overlap with ``Y`` (from most to least significant);
``\alpha \in (0, 1]`` is the *enrichment discount*,
and ``\beta > 0`` is the *enrichment shape*.

The discount parameter ``\alpha`` is important for complex data (``M_X \gg 1``).
It allows terms that are very significant, but only in a few conditions,
be preferred over the terms that are moderately represented in most conditions.
The lower is ``\alpha``, the more enrichment score favors condition-specific terms.

The *enrichment shape* parameter ``\beta`` affects the size of the selected annotation terms.
For small ``\beta`` the method would prefer several smaller and more specific terms over
the single one that combines multiple small terms.

### Redundancy score
```math
q_r(C) = \sum_{i=1}^{N_T} \sum_{j=1}^{N_T} \min(c_i, c_j) R(\mathrm{obs}(A_i), \mathrm{obs}(A_j)),
```
where ``R(\mathrm{obs}(A_i), \mathrm{obs}(A_j))`` measures the redundancy of observed entities from
``A_i`` and ``A_j`` annotation terms as the enrichment of one term in their union:
```math
R(B_i, B_j) = -\left( -\log \min \left(
  P_{\mathrm{FET}}(B_i, B_i \cup B_j, N_O (1 + \Delta)),
  P_{\mathrm{FET}}(B_j, B_i \cup B_j, N_O (1 + \Delta)) \right)\right)^\beta,
```
where ``\Delta N_O \ge 0`` is added to keep the overlap significant for very large annotation terms that are
comparable in size to the set of all observed genes ``O``.

### Hits cover scores
``q_u`` counts how many hit entities were *not covered* by the selected terms:
```math
q_u(C, \mathcal{X}) = \sum_{o \in \cup_i X_i} \left(1 - \max_{j \in С(o)} c_j\right),
```
where ``С(o)`` is the set of all indices ``j``, s.t. ``o \in \mathrm{obs}(A_j)``.

Similarly, ``q_c`` counts how many non-hit observations in *each* experiment are *covered* by
the selected terms:
```math
q_c(C, \mathcal{X}) = \sum_{i=1}^{M_X} \sum_{o \in O\setminus X_i} \max_{j \in С(o)} c_j.
```

### Optimal covers
With ``q = (q_e, q_u, q_c, q_r)`` cover quality defined, one can consider
the *multi-objective optimization problem* of finding the
[Pareto-optimal](https://en.wikipedia.org/wiki/Pareto_efficiency) family of
covers ``\mathcal{C}_{*} = \{ C_{*, 1}, C_{*, 2}, \ldots \}``:
```math
\forall C \not\in \mathcal{C}_{*} \;\exists C_{*} \in \mathcal{C}_{*}: 
q_e(C_*) \leq q_e(C), q_r(C_*) \leq q_r(C), q_u(C_*) \leq q_u(C_*), q_c(C_*) \leq q_c(C),
```
and at least for one of the ``q`` components the inequality is strict.

The resulting Pareto front ``\mathcal{C}_{*}`` will contain: the covers with minimally redundant but less enriched terms,
very redundant covers of more enriched terms, as well as intermediate solutions between these extremes.
It's possible to weight the importance of individual quality components
(``w_e = 1``, ``w_r > 0``, ``w_u \ge 0``, ``w_c \ge 0``) to identify the unique optimal
cover ``C_{*, \Sigma} \in \mathcal{C}_{*}`` that minimizes the components weighted sum:
```math
q_{\Sigma}(C) = q_e(C) + w_r q_r(C) + w_u q_u(C) + w_c q_c(C) \to \min, \quad C \in \mathcal{C}_{*}.
```

### [Components convolution](@id cover_score_convolution)
The number of points on the approximated 4-objective Pareto front is ``O(\frac{1}{\delta^{(4-1)}})``,
where ``\delta > 0`` is the grid step in the multi-objective space. The maintenance of such
front creates significant overhead during optimization.
However, as the components ``q_u`` and ``q_c`` correlate with the enrichment component ``q_e``,
they could be convoluted with the ``q_e``:
```math
\tilde{q}_e = q_e + w_u q_u + w_c q_c,
```
or omitted altogether by setting ``w_u = w_c = 0``.
The optimization of the simplified 2-objective ``(\tilde{q}_e, q_r)`` problem would require much less resources.
The corresponding Pareto front would essentially show the trade-offs between
the annotation terms enrichment ``\tilde{q}_e`` and their redundancy ``q_r``.
One can define the *redunancy index* ``k_r`` as
```math
k_r(C) = \frac{q_r(C)}{-\tilde{q}_e(C)}.
```
High ``k_r(C)`` means the terms in the cover ``C`` are redundant. In practice, covers with
``k_r(C) \gg 1``, although marginally improving ``q_e``, do not reveal any new patterns in the data.
To improve the optimization efficiency, it's possible to use the modified cover score
``\tilde{q}``, so that the covers with high ``k_r`` become dominated by less redundant covers,
and thus keep only the relevant solutions on the Pareto front:
```math
\tilde{q} = \left(\tilde{q}_e + w_r t q_r, (1 - t) q_r \right),\;\text{where}\\
t = \frac{1}{k_r} \left(w_r + \exp\left(s^{-1} - s\right)\right)^{-1}, \quad
s = \alpha_k \max(k_r - k_{\max}, 0),
```
``k_{\max}`` is the maximal accepted redundancy index (by default the method uses ``k_{\max} = 1``),
and ``\alpha_k > 0`` defines the transformation strength.
With this definition, when ``k_r(C) \leq k_{\max}``,
then ``s = 0``, ``t = 0``, and both score components are unchanged: ``\tilde{q} = (\tilde{q}_e, q_r)``.
When ``k_r(C) \gg k_{\max}``, i.e. when the redundancy
of the cover is very high, ``t \approx \frac{1}{k_r w_r} = \frac{-\tilde{q}_e}{w_r q_r}``,
and ``\tilde{q} \approx (0, q_r + \frac{\tilde{q}_e}{w_r})``.
The latter solution is dominated by any less redundant cover. 
Also note that, for any ``C``, ``\tilde{q}_1 + w_r \tilde{q}_2 = q_{\Sigma}``, so this 
transformation preserves the optimal cover ``C_{*,\Sigma}`` that minimizes ``q_{\Sigma}(C)``.
