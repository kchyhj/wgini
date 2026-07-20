{smcl}
{* *! version 1.1.0  17jul2026}{...}
{vieweralsosee "[R] inequality" "help inequality"}{...}
{vieweralsosee "ineqdeco (if installed)" "help ineqdeco"}{...}
{viewerjumpto "Syntax" "wgini##syntax"}{...}
{viewerjumpto "Description" "wgini##description"}{...}
{viewerjumpto "Options" "wgini##options"}{...}
{viewerjumpto "Method" "wgini##method"}{...}
{viewerjumpto "Stored results" "wgini##results"}{...}
{viewerjumpto "Examples" "wgini##examples"}{...}
{viewerjumpto "References" "wgini##references"}{...}
{viewerjumpto "Author" "wgini##author"}{...}
{title:Title}

{phang}
{bf:wgini} {hline 2} Weighted Gini coefficient with Lerman-Yitzhaki source decomposition


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:wgini} {varname} {ifin} {weight} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt source(varlist)}}decompose the Gini into the additive
contributions of these sources; their sum must equal {varname}{p_end}
{synopt:{opt gi(newvar)}}store the observation-level Gini contribution,
which sums to the Gini{p_end}
{synopt:{opt top(numlist)}}top-share diagnostics: for each {it:p}, the top
{it:p}%'s population share, value share, Gini share, and the Gini without
them{p_end}

{syntab:Details}
{synopt:{opt tol(#)}}tolerance, in the units of the data, for the check that
the sources sum to {varname}; default {cmd:tol(1)}; {cmd:tol(0)} skips it{p_end}
{synopt:{opt nopr:int}}suppress the display{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{cmd:aweight}s, {cmd:pweight}s, and {cmd:fweight}s are allowed; see
{help weight}. With no weight every observation is weighted equally.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:wgini} computes the Gini coefficient of {varname} using the
covariance form of {help wgini##LY1984:Lerman and Yitzhaki (1984)},
G = 2 cov(x, F(x)) / mu, applied to individual records with weights as in
{help wgini##LY1989:Lerman and Yitzhaki (1989)}, which admits sampling or
population weights directly. The variable may take negative values (for
example net worth, assets minus debt): unlike the Lorenz-curve or
mean-absolute-difference constructions, which implicitly assume
nonnegative values, the covariance form needs only each observation's
deviation from the mean and its fractional rank, both defined for any
real values — the only requirement is a positive mean. Negative
observations are therefore used as-is, with no truncation at zero; the
Gini is then no longer bounded by 1 and may exceed it. This is a
documented feature, not an error.

{pstd}
With {opt source()}, when {varname} is the sum of {it:K} sources
{it:x} = {it:y}{sub:1} + {it:...} + {it:y}{sub:K}, {cmd:wgini} decomposes the
Gini additively into each source's contribution and factors each contribution
into share, own Gini, and Gini correlation with the rank of the total. A
source entered with a negative mean (for example {bf:-}debt) contributes
negatively; no absolute values are taken.

{pstd}
With {opt gi()}, {cmd:wgini} stores each observation's own contribution to the
Gini, which sums exactly to the coefficient. This is useful for asking what a
single unit contributes, or what the Gini becomes when a set of units is
removed.

{pstd}
{bf:Ties.} Tied values of the ranking variable all receive the mid-rank of
their tie group, ({it:cumulative weight before the group} + 0.5 {c 215}
{it:group weight}) / {it:W}. Without this, tied observations would be ranked
in whatever order the sort happens to leave them, which is random in Stata:
the total Gini is provably invariant to that order, but source contributions
and observation-level contributions are not, and results would change from
run to run. Mid-ranking makes the decomposition well defined and exactly
reproducible; the total Gini is unchanged. {cmd:wgini} also restores the
original sort order of the data before exiting.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt source(varlist)} decomposes the Gini into the contributions of the listed
sources. The sources must add up to {varname} at every observation; this is
checked (see {opt tol()}). The decomposition and its factors are returned in
{cmd:r(decomp)} and displayed.

{phang}
{opt gi(newvar)} creates {it:newvar} holding the observation-level Gini
contribution

{p 12 12 2}
{it:g_i} = {it:w_i} {c 215} 2({it:x_i} {c 45} {c 956})({it:p_i} {c 45} 0.5) /
({c 956} {it:W}),

{pmore}
where {it:p_i} is the weighted fractional rank, {c 956} the weighted mean, and
{it:W} the total weight. The sum of {it:newvar} equals {cmd:r(gini)}.

{dlgtab:Details}

{phang}
{opt tol(#)} sets the tolerance, in the units of {varname}, for the check that
the sources sum to the total. The default {cmd:tol(1)} allows rounding at the
last recorded unit of the data. Increase it when the sources come from coarser
rounding; {cmd:tol(0)} skips the check, which is safe only when the sources
sum to the total by construction.

{phang}
{opt top(numlist)} reports top-share diagnostics for each {it:p} in the
list (0 < {it:p} < 100): the top group is defined by the weighted
fractional mid-rank exceeding 1 {c 45} {it:p}/100. Results are returned in
{cmd:r(top)} and as per-{it:p} scalars (see {help wgini##results:stored results});
the scalars exist so that {helpb statsby} can collect them by subgroup.

{phang}
{opt noprint} suppresses the display. The stored results are produced either
way.


{marker method}{...}
{title:Method}

{pstd}
Let {it:p_i} be the weighted fractional rank of {varname},

{p 8 8 2}
{it:p_i} = (running weight up to {it:i} {c 45} 0.5{it:w_i}) / {it:W},
{space 3}{it:W} = total weight,

{pstd}
after sorting by {varname}. The weighted mean of {it:p} is exactly 0.5, so
the Gini is

{p 8 8 2}
{it:G} = 2 cov{sub:w}({it:x}, {it:p}) / {c 956}
= (2 / ({c 956} {it:W})) {c 8721}{sub:i} {it:w_i}({it:x_i} {c 45} {c 956})({it:p_i} {c 45} 0.5).

{pstd}
{bf:Source decomposition.} Because covariance is linear in its first argument,
the Gini of a sum of sources is the sum of the sources' contributions:

{p 8 8 2}
{it:G} = {c 8721}{sub:k} 2 cov{sub:w}({it:y_k}, {it:p}) / {c 956}
= {c 8721}{sub:k} {it:S_k G_k R_k},

{pstd}
where the rank {it:p} and mean {c 956} are those of the {it:total}, and

{p 8 12 2}
{it:S_k} = {c 956}{sub:k} / {c 956} {space 8}(share of the total){break}
{it:G_k} = 2 cov{sub:w}({it:y_k}, {it:p_k}) / {c 956}{sub:k} {space 3}(Gini of the source, on its own rank {it:p_k}){break}
{it:R_k} = cov{sub:w}({it:y_k}, {it:p}) / cov{sub:w}({it:y_k}, {it:p_k}) {space 3}(Gini correlation with the total's rank).

{pstd}
A large contribution can come from a large share {it:S_k}, from a highly
concentrated source (large {it:G_k}), or from a source that tracks the ranking
of the total (large {it:R_k}); the three factors separate these.

{pstd}
{bf:Validation.} {cmd:descogini} (Lopez-Feldman 2006) performs the same
decomposition but does not accept weights ({cmd:weights not allowed}), and
{cmd:ineqrbd} (Jenkins 1999) implements a different, regression-based
(Fields, Shorrocks) decomposition. Run with all weights equal to one,
{cmd:wgini}'s {it:S_k}, {it:G_k}, {it:R_k}, and share reproduce
{cmd:descogini} to four decimals.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wgini} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(gini)}}weighted Gini coefficient{p_end}
{synopt:{cmd:r(N)}}number of observations used{p_end}
{synopt:{cmd:r(mean)}}weighted mean of {varname}{p_end}
{synopt:{cmd:r(sumdev)}}sum of contributions minus {cmd:r(gini)} (identity
check); with {opt source()} only{p_end}
{synopt:{cmd:r(actual_}{it:p}{cmd:)}, {cmd:r(vshare_}{it:p}{cmd:)}, {cmd:r(gshare_}{it:p}{cmd:)}, {cmd:r(gexcl_}{it:p}{cmd:)}}per-{it:p}
scalars from {opt top()}, named by the value ({cmd:.} {c 45}> {cmd:_}){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(sources)}}the source {varlist}; with {opt source()} only{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(decomp)}}{it:K} {c 215} 5 matrix, one row per source, columns
{cmd:contrib share Sk Gk Rk}; with {opt source()} only{p_end}
{synopt:{cmd:r(top)}}{it:K} {c 215} 5 matrix, one row per {it:p}, columns
{cmd:top_pct actual_pct value_share gini_share gini_excl}; with {opt top()} only{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
{bf:1. Weighted Gini of one variable.} Computes the Gini of {cmd:networth}
with the survey weight applied and stores {cmd:r(gini)}, {cmd:r(N)}, and
{cmd:r(mean)}.{p_end}

{phang2}{cmd:. wgini networth [aw=weight]}{p_end}

{pstd}
{bf:2. Source decomposition.} When net worth is the sum of asset
components — debt entered as a negative variable
({cmd:gen negdebt = -debt}) — this splits the Gini into one additive
contribution per component and factors each into
{it:S_k} {c 215} {it:G_k} {c 215} {it:R_k}. The full table is also returned
in {cmd:r(decomp)}.{p_end}

{phang2}{cmd:. wgini networth [aw=weight], source(home other_re other_real stockbondfund othersaving negdebt)}{p_end}
{phang2}{cmd:. matrix list r(decomp)}{p_end}

{pstd}
{bf:3. Top-share diagnostics — who drives the Gini?} This answers
questions like: {it:how much of measured inequality comes from the top
1%? And what would the Gini be without them?} One option does it for any
list of top shares:{p_end}

{phang2}{cmd:. wgini networth [aw=weight], top(1 5 10)}{p_end}
{phang2}{cmd:. matrix list r(top)}{p_end}

{pstd}
For each {it:p} the command reports one row with four numbers.
{cmd:actual_pct} is the top group's actual weighted population share:
the top group is every observation whose weighted fractional mid-rank
exceeds 1 {c 45} {it:p}/100 (a rank cut), and because the weighted data
are discrete and a tie group — sharing one mid-rank — stays together on
one side, this can differ slightly from {it:p}. {cmd:value_share} is the
share of the weighted total held by the top group. {cmd:gini_share} is
the share of the Gini contributed by the top group: each observation has
an additive contribution {it:g_i} that sums exactly to the Gini, and this
column sums the top group's {it:g_i} and divides by {it:G}.
{cmd:gini_excl} is the Gini {it:recomputed without} the top group, on the
remaining sample with its own mean and its own ranks.{p_end}

{pstd}
Reading {cmd:value_share} against {cmd:gini_share} shows how
disproportionate the top is; comparing {cmd:gini_excl} with the full Gini
shows how much of measured inequality hinges on a thin top slice. Two
warnings: {cmd:gini_excl} is {bf:not} the Gini minus {cmd:gini_share} —
deleting observations changes the mean and every rank, so the exclusion
is a genuine recomputation. And contributions are positive in {it:both}
tails — for the rich both ({it:x_i} {c 45} {c 956}) and
({it:F_i} {c 45} 0.5) are positive, for the poor both are negative — so
the {it:bottom} 1% also contributes positively; {opt top()} looks at the
top by construction.{p_end}

{pstd}
{bf:4. Custom sets beyond top shares.} When the group of interest is not
a top share — one specific household, a region, an occupation —
{cmd:gi()} stores each observation's contribution {it:g_i} as a variable,
and the sum over any set divided by {cmd:r(gini)} is that set's share of
the Gini. For the single richest household ({helpb gsort} with a minus
sign sorts in descending order, so observation 1 is the top):{p_end}

{phang2}{cmd:. wgini networth [aw=weight], gi(gcon) noprint}{p_end}
{phang2}{cmd:. scalar G = r(gini)}{p_end}
{phang2}{cmd:. gsort -networth}{p_end}
{phang2}{cmd:. display gcon[1]/G}{p_end}

{pstd}
To {it:exclude} an arbitrary set, restrict the sample with {cmd:if}
({cmd:... if hhid != "<top id>"}).{p_end}

{pstd}
{bf:5. Gini by subgroup, with statsby.} {cmd:wgini} computes one set of
results for the sample it is given. To repeat it per subgroup, use
{helpb statsby} — generic Stata machinery that works with any command and
is not part of {cmd:wgini}. It runs the command once per group and
collects the returned scalars, one row per group. The {opt top()} results
are returned as scalars ({cmd:r(gshare_1)} etc.) precisely so that
{cmd:statsby} can collect them. Note that {cmd:statsby, clear} replaces
the data in memory — save your data first.{p_end}

{phang2}{cmd:. statsby gini=r(gini) n=r(N), by(year agegrp) clear: wgini networth [aw=weight]}{p_end}
{phang2}{cmd:. statsby gini=r(gini) gsh=r(gshare_1) gex=r(gexcl_1), by(year) clear: wgini networth [aw=weight], top(1)}{p_end}

{pstd}
{it:Caveat.} This {it:computes} a separate Gini within each subgroup; it
does not {it:decompose} the overall Gini into subgroups. The Gini
decomposes exactly by income or asset source (example 2), but not by
population subgroup: when subgroup distributions overlap,
"within + between" does not add up to the total — a residual overlap term
remains. For an exact within/between subgroup decomposition use a
generalized entropy index (e.g. Theil), for example with {cmd:ineqdeco}
(Jenkins, SSC).{p_end}


{marker references}{...}
{title:References}

{marker LY1984}{...}
{phang}
Lerman, R. I., and S. Yitzhaki. 1984. A note on the calculation and
interpretation of the Gini index. {it:Economics Letters} 15: 363-368.

{marker LY1985}{...}
{phang}
Lerman, R. I., and S. Yitzhaki. 1985. Income inequality effects by income
source: A new approach and applications to the United States.
{it:Review of Economics and Statistics} 67: 151-156.

{marker LY1989}{...}
{phang}
Lerman, R. I., and S. Yitzhaki. 1989. Improving the accuracy of estimates of
Gini coefficients. {it:Journal of Econometrics} 42: 43-47.

{phang}
Lopez-Feldman, A. 2006. Decomposing inequality and obtaining marginal effects.
{it:Stata Journal} 6: 106-111.

{phang}
Stark, O., J. E. Taylor, and S. Yitzhaki. 1986. Remittances and inequality.
{it:Economic Journal} 96: 722-740.


{marker author}{...}
{title:Author}

{pstd}ChangHwan Kim, University of Kansas.{p_end}

{pstd}License: PolyForm Noncommercial 1.0.0 — noncommercial use, modification,
and distribution are free; commercial use requires a separate license.
If you use {cmd:wgini} in published work, please cite the repository:
{browse "https://github.com/kchyhj/wgini"}.{p_end}
