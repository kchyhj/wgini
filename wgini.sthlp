{smcl}
{* *! version 1.0.0  16jul2026}{...}
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
Lerman-Yitzhaki (1984, 1989) covariance form, which admits sampling or
population weights directly. The variable may take negative values (for
example net worth, assets minus debt); the Gini is then still defined and
may exceed 1. This is a documented feature, not an error.

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

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(sources)}}the source {varlist}; with {opt source()} only{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(decomp)}}{it:K} {c 215} 5 matrix, one row per source, columns
{cmd:contrib share Sk Gk Rk}; with {opt source()} only{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
{bf:1. Weighted Gini of one variable.} Computes the Gini of {cmd:networth}
with the survey weight applied and stores {cmd:r(gini)}, {cmd:r(N)}, and
{cmd:r(mean)}.{p_end}

{phang2}{cmd:. wgini networth [aw=weight]}{p_end}

{pstd}
{bf:2. Gini computed by subgroup.} {cmd:wgini} returns one Gini for the
sample it is given. To get one Gini per subgroup, wrap it in
{helpb statsby}, which runs the command once for every subgroup defined by
{cmd:by()} and collects the returned scalars into a new dataset with one
row per subgroup. Note that {cmd:statsby, clear} replaces the data in
memory — save your data first.{p_end}

{phang2}{cmd:. statsby gini=r(gini) n=r(N), by(year agegrp) clear: wgini networth [aw=weight]}{p_end}
{phang2}{cmd:. list year agegrp gini n}{p_end}

{pstd}
{it:Caveat.} This {it:computes} a separate Gini within each subgroup; it
does not {it:decompose} the overall Gini into subgroups. The reason is in
the covariance form {it:G} = 2 cov{sub:w}({it:x}, {it:F}({it:x}))/{c 956}.
Sources enter the {it:first} argument: with
{it:x} = {it:y}{sub:1} + ... + {it:y}{sub:K}, linearity of the covariance
in its first argument, with {it:F}({it:x}) and {c 956} of the {it:total}
held fixed, gives cov{sub:w}({c 931}{sub:k} {it:y_k}, {it:F}({it:x})) =
{c 931}{sub:k} cov{sub:w}({it:y_k}, {it:F}({it:x})), so {it:G} =
{c 931}{sub:k} 2 cov{sub:w}({it:y_k}, {it:F}({it:x}))/{c 956} exactly —
that is what {opt source()} computes. Subgroups instead act on the
{it:second} argument, the ranking: a household's within-group rank
{it:F_g}({it:x}) differs from its population rank {it:F}({it:x}) whenever
group distributions overlap, and the covariance is not separable in the
ranking argument, so "within + between" does not add up to {it:G} — a
residual overlap term remains. For an exact within/between subgroup
decomposition use a generalized entropy index (e.g. Theil), for example
with {cmd:ineqdeco} (Jenkins, SSC).{p_end}

{pstd}
{bf:3. Source decomposition.} When net worth is the sum of asset
components — debt entered as a negative variable
({cmd:gen negdebt = -debt}) — this splits the Gini into one additive
contribution per component and factors each into
{it:S_k} {c 215} {it:G_k} {c 215} {it:R_k}. The full table is also returned
in {cmd:r(decomp)}.{p_end}

{phang2}{cmd:. wgini networth [aw=weight], source(home other_re other_real stockbondfund othersaving negdebt)}{p_end}
{phang2}{cmd:. matrix list r(decomp)}{p_end}

{pstd}
{bf:4. Observation-level contributions.} {cmd:gi(gcon)} creates the
variable {cmd:gcon} holding each observation's additive contribution to
the Gini; the contributions sum exactly to {cmd:r(gini)}, so
{cmd:gcon/r(gini)} is the share of overall inequality attributable to that
observation. {helpb gsort} with a minus sign sorts in descending order, so
after {cmd:gsort -networth} observation 1 is the richest household, and
{cmd:in 1} restricts {cmd:list} to that first row.{p_end}

{phang2}{cmd:. wgini networth [aw=weight], gi(gcon) noprint}{p_end}
{phang2}{cmd:. scalar G = r(gini)}{p_end}
{phang2}{cmd:. gsort -networth}{p_end}
{phang2}{cmd:. display gcon[1]/G}{p_end}
{phang2}{cmd:. list networth gcon in 1}{p_end}

{pstd}
{bf:5. Recomputing the Gini without the top unit(s).} The contribution in
{cmd:gcon} describes the observation's role {it:within the current sample};
it is {bf:not} the amount the Gini would fall by if the observation were
deleted, because deleting it changes the mean and every rank. To get the
Gini {it:without} some units, rerun {cmd:wgini} on the restricted sample
with {cmd:if}. Without the single largest household ({cmd:r(max)} from
{cmd:summarize} is the sample maximum):{p_end}

{phang2}{cmd:. quietly sum networth}{p_end}
{phang2}{cmd:. wgini networth if networth < r(max) [aw=weight]}{p_end}

{pstd}
Without the top 1%: {helpb _pctile} with the weight computes the weighted
99th percentile, returned in {cmd:r(r1)}, and {cmd:if} keeps the households
at or below that cutoff.{p_end}

{phang2}{cmd:. _pctile networth [aw=weight], p(99)}{p_end}
{phang2}{cmd:. wgini networth if networth <= r(r1) [aw=weight]}{p_end}

{pstd}
If several households tie exactly at the maximum, the first {cmd:if} drops
all of them; to pin down one specific household, condition on its id
instead ({cmd:... if hhid != "<top id>"}).{p_end}


{marker references}{...}
{title:References}

{phang}
Lerman, R. I., and S. Yitzhaki. 1984. A note on the calculation and
interpretation of the Gini index. {it:Economics Letters} 15: 363-368.

{phang}
Lerman, R. I., and S. Yitzhaki. 1985. Income inequality effects by income
source: A new approach and applications to the United States.
{it:Review of Economics and Statistics} 67: 151-156.

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
