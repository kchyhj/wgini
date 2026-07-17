# wgini — Weighted Gini coefficient with Lerman–Yitzhaki source decomposition

`wgini` is a Stata command that computes the Gini coefficient under sampling
or population weights and, optionally, decomposes it additively by income or
asset **source** (Lerman & Yitzhaki 1985) or by **observation**. It exists
because the standard tools cover only parts of this:

| | weights | LY source decomposition | observation-level contribution |
|---|---|---|---|
| `descogini` (López-Feldman 2006, *Stata Journal*) | **no** (`weights not allowed`) | yes | no |
| `ineqrbd` (Jenkins) | yes | no — regression-based (Fields/Shorrocks) CV decomposition, a different method | no |
| `wgini` | yes | yes | yes |

Survey microdata (CPS, LIS, SHFLC, EU-SILC, …) require weights, which is what
motivated this command.

## Installation

```stata
net install wgini, from("https://raw.githubusercontent.com/kchyhj/wgini/main/")
```

or copy `wgini.ado` and `wgini.sthlp` into your personal ado directory
(`adopath` shows where). Requires Stata 14 or later. After installing, see
`help wgini`.

## Quick start

**1. Weighted Gini of one variable.** Computes the Gini of `networth` with
the survey weight applied, and leaves `r(gini)`, `r(N)`, and `r(mean)`
behind for later use.

```stata
wgini networth [aw=weight]
```

**2. Source decomposition.** When net worth is the sum of asset components
(here housing + other real estate + vehicles + savings + debt entered as a
*negative* variable, e.g. `gen negdebt = -debt`), this splits the Gini into
one additive contribution per component and factors each into
share × own Gini × Gini correlation. The decomposition is an exact
identity — see [Method](#method) for why.

```stata
wgini networth [aw=weight], source(home other_re vehicles saving negdebt)
matrix list r(decomp)     // contrib share Sk Gk Rk, one row per source
```

**3. Top-share diagnostics — who drives the Gini?** This answers
questions like: *how much of measured inequality comes from the top 1%?
And what would the Gini be without them?* One option does it for any list
of top shares:

```stata
wgini networth [aw=weight], top(1 5 10)
matrix list r(top)
```

For each `p` the command reports one row with four numbers:

| column | meaning |
|---|---|
| `actual_pct` | the top group's actual weighted population share. The top group is everyone whose weighted fractional mid-rank exceeds 1−p/100 (a rank cut); because the weighted data are discrete and a tie group — sharing one mid-rank — stays together on one side, this can differ slightly from `p`. |
| `value_share` | the share of total wealth (or income) held by the top group. |
| `gini_share` | the share of the Gini contributed by the top group: each observation has an additive contribution $g_i$ that sums exactly to the Gini, and this column sums the top group's $g_i$ and divides by $G$. |
| `gini_excl` | the Gini *recomputed without* the top group, on the remaining sample with its own mean and its own ranks. |

Reading `value_share` against `gini_share` shows how disproportionate the
top is: in one application to a national wealth survey, the richest
household among 20-something householders held 12% of the group's assets
but produced 17% of the group's Gini. Comparing `gini_excl` to the full
Gini shows how much of measured inequality hinges on a thin top slice.

Two warnings. `gini_excl` is **not** the Gini minus `gini_share`:
deleting observations changes the mean and every rank, so the exclusion
is a genuine recomputation. And contributions are positive in *both*
tails — for the rich, $(x_i-\mu)>0$ and $(F_i-\tfrac12)>0$; for the poor,
both factors are negative — so the *bottom* 1% also contributes
positively to the Gini; `top()` looks at the top by construction.

**4. Custom sets beyond top shares.** When the group of interest is not
a top share — one specific household, a region, an occupation — `gi()`
stores each observation's contribution $g_i$ as a variable, and the sum
over any set divided by `r(gini)` is that set's share of the Gini. For
the single richest household:

```stata
wgini networth [aw=weight], gi(gcon) noprint
scalar G = r(gini)
gsort -networth              // descending sort: observation 1 is the top
display gcon[1]/G            // its share of the Gini
```

To *exclude* an arbitrary set, restrict the sample with `if`
(e.g. `wgini networth if hhid != "<top id>" [aw=weight]`).

## Computing the Gini by subgroup

`wgini` computes one set of results for the sample it is given. To repeat
it for every subgroup (say, year × age group), use Stata's `statsby` —
generic Stata machinery that works with any command and is **not part of
`wgini`**. It runs the command once per group and collects the returned
scalars into a new dataset, one row per group:

```stata
statsby gini=r(gini) n=r(N), by(year agegrp) clear: ///
    wgini networth [aw=weight]
list year agegrp gini n
```

The `top()` results are returned as scalars precisely so `statsby` can
collect them by subgroup — for example, the top-1% Gini share and the
Gini without the top 1%, per year:

```stata
statsby gini=r(gini) gsh=r(gshare_1) gex=r(gexcl_1), by(year) clear: ///
    wgini networth [aw=weight], top(1)
```

(`statsby ..., clear` replaces the data in memory with the collected
results — save your data first.)

> **Caveat.** This *computes* a separate Gini within each subgroup; it
> does not *decompose* the overall Gini into subgroups. The Gini
> decomposes exactly by income/asset source (example 2), but not by
> population subgroup: when subgroup distributions overlap,
> "within + between" does not add up to the total — a residual overlap
> term remains. For an exact within/between subgroup decomposition use a
> generalized entropy index (e.g. Theil), for example with `ineqdeco`
> (Jenkins, SSC).

## Method

The Gini is computed in the Lerman–Yitzhaki (1984) covariance form

$$G \;=\; \frac{2\,\mathrm{cov}_w\!\big(x,\,F(x)\big)}{\mu},$$

where $F(x)$ is the weighted fractional rank and $\mu$ the weighted mean.

**Why the source decomposition is exact.** Sources enter the *first*
argument of the covariance. If $x = \sum_k y_k$, then because the
covariance is linear in its first argument while $F(x)$ and $\mu$ — the
ranking and the mean of the *total* — stay fixed,

$$\mathrm{cov}_w\!\Big(\textstyle\sum_k y_k,\,F(x)\Big)
= \sum_k \mathrm{cov}_w\big(y_k,\,F(x)\big)
\quad\Longrightarrow\quad
G = \sum_k \frac{2\,\mathrm{cov}_w\big(y_k,\,F(x)\big)}{\mu},$$

an exact identity with one term per source. Each term factors as

$$\frac{2\,\mathrm{cov}_w\big(y_k,\,F(x)\big)}{\mu} = S_k\,G_k\,R_k,$$

with $S_k$ the source's share of the total, $G_k$ the source's own Gini,
and $R_k$ its Gini correlation with the rank of the total (Lerman &
Yitzhaki 1985; Stark, Taylor & Yitzhaki 1986). The identity
$\sum_k \text{contrib}_k = G$ is returned as a check (`r(sumdev)`); in the
test suite it holds to ~1e-16. A large contribution can come from a large
share, from a concentrated source (large $G_k$), or from a source that
tracks the overall ranking (large $R_k$) — the three factors separate
these.

**Negative values are allowed throughout — a property of the
Lerman–Yitzhaki form.** Textbook constructions of the Gini (the Lorenz
curve, or the mean of absolute differences divided by twice the mean)
implicitly assume nonnegative values. The covariance form does not: it
needs only each observation's deviation from the mean, $x_i - \mu$, and
its fractional rank, $F(x_i)$, both of which are perfectly well defined
for negative values. The only requirement is a positive mean, $\mu > 0$.
So a variable like net worth, where indebted households are negative, is
computed as-is — no truncation at zero, no dropping of negative
observations. One consequence to be aware of: with negative values the
Gini is no longer bounded by 1 and can exceed it (documented, not an
error). Likewise a source with a negative mean (enter debt as `-debt`)
has negative $S_k$, $G_k$, $R_k$ whose product is still its correct
negative contribution. No absolute values are taken anywhere.

### Ties and reproducibility

Tied values receive the **mid-rank** of their tie group. This matters more
than it sounds: Stata's `sort` orders tied observations randomly, and while
the total Gini is provably invariant to that order, source contributions and
observation-level contributions are not — with arbitrary tie ranks they change
from run to run (in our survey application, by up to ~1e-4 in small groups
with many zero-wealth ties). Mid-ranking makes the decomposition well defined
and exactly reproducible: same data in, bit-identical results out. `wgini`
also restores the sort order of your data before exiting.

## Validation

The test suite (`test/wgini_test.do`, self-contained synthetic data) checks:

1. **Identity**: source contributions sum to the Gini (~1e-16).
2. **Factorization**: `contrib_k = S_k * G_k * R_k` holds row by row,
   including for a negative-mean (debt) source.
3. **Against `descogini`**: with all weights equal to 1 — the only case
   `descogini` accepts — `wgini` reproduces its `Sk`, `Gk`, `Rk`, and share
   columns to four decimals. (On a national wealth survey with 18,664
   households and six sources the match was exact to four decimals in every
   cell.) With actual weights the two commands are not comparable, since
   `descogini` refuses weights.
4. **Tie invariance**: permuting the physical order of tied observations
   leaves every result unchanged.
5. **Reproducibility**: two runs on the same data return bit-identical
   `r(decomp)`.
6. **Sort restoration**: the command does not change the order of your data.
7. **Error handling**: sources that do not sum to the total are rejected
   (r 459); an all-zero source gets contribution 0 with `G_k`, `R_k` missing.
8. **`top()` equals the manual recipes**: its Gini share matches the
   `gi()`-based sum, its excluded Gini matches a rerun with `if`, and its
   value share matches direct aggregation.
9. **`statsby` collection**: the `top()` scalars collected by subgroup
   reproduce direct by-group calls.

## Stored results

| | |
|---|---|
| `r(gini)` | weighted Gini |
| `r(N)`, `r(mean)` | observations used, weighted mean |
| `r(decomp)` | K×5 matrix: `contrib share Sk Gk Rk` (with `source()`) |
| `r(sumdev)` | sum of contributions − Gini, identity check (with `source()`) |
| `r(sources)` | the source varlist (with `source()`) |
| `r(top)` | K×5 matrix: `top_pct actual_pct value_share gini_share gini_excl` (with `top()`) |
| `r(actual_p)`, `r(vshare_p)`, `r(gshare_p)`, `r(gexcl_p)` | scalar versions per `p` in `top()`, named by the value (`.` → `_`), for `statsby` |

## References

- Lerman, R. I., and S. Yitzhaki. 1984. A note on the calculation and
  interpretation of the Gini index. *Economics Letters* 15: 363–368.
- Lerman, R. I., and S. Yitzhaki. 1985. Income inequality effects by income
  source: A new approach and applications to the United States. *Review of
  Economics and Statistics* 67: 151–156.
- López-Feldman, A. 2006. Decomposing inequality and obtaining marginal
  effects. *Stata Journal* 6: 106–111.
- Stark, O., J. E. Taylor, and S. Yitzhaki. 1986. Remittances and inequality.
  *Economic Journal* 96: 722–740.

## Citing wgini

If you use `wgini` in published work, please cite it (see
[CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button gives
BibTeX and APA):

> Kim, ChangHwan. 2026. *wgini: Weighted Gini coefficient with
> Lerman–Yitzhaki source decomposition* (Version 1.1.0) [Stata command].
> https://github.com/kchyhj/wgini

## Author

ChangHwan Kim, University of Kansas
(ORCID [0000-0001-7149-1386](https://orcid.org/0000-0001-7149-1386))

**AI disclosure.** The code and documentation were drafted with the
assistance of Claude (Anthropic; Claude Code, Opus 4.8), working under the
author's direction. The statistical design, the choice of estimator and
tie treatment, and all validation decisions are the author's. Every release
is verified by the test suite in `test/` and against `descogini` before
publication.

한국어 안내는 [README.ko.md](README.ko.md)를 보세요.

## License

[PolyForm Noncommercial 1.0.0](LICENSE.md) — free to use, modify, and share
for any **noncommercial** purpose (research, teaching, government, and
nonprofit use included, regardless of funding source). Commercial use
requires a separate license from the author.

Required Notice: Copyright ChangHwan Kim (https://github.com/kchyhj/wgini)
