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

```stata
* weighted Gini
wgini networth [aw=weight]

* by cell
statsby gini=r(gini) n=r(N), by(year agegrp): wgini networth [aw=weight]

* source decomposition: sources must sum to the total
wgini networth [aw=weight], source(home other_re vehicles saving negdebt)

* observation-level contributions (sum to the Gini)
wgini networth [aw=weight], gi(gcon) noprint
gsort -networth
list networth gcon in 1   // what does the top household contribute?
```

## Method

The Gini is computed in the Lerman–Yitzhaki (1984) covariance form

```
G = 2 cov_w(x, F) / mu
```

where `F` is the weighted fractional rank and `mu` the weighted mean. With
`source(y1 ... yK)` and `x = y1 + ... + yK`, linearity of the covariance gives
an exact additive decomposition,

```
G = sum_k 2 cov_w(y_k, F) / mu = sum_k S_k * G_k * R_k
```

with `S_k` the source's share of the total, `G_k` the source's own Gini, and
`R_k` its Gini correlation with the rank of the total (Lerman & Yitzhaki
1985; Stark, Taylor & Yitzhaki 1986). The identity `sum_k contrib_k = G` is
returned as a check (`r(sumdev)`); in the test suite it holds to ~1e-16.

Negative values are allowed throughout. A variable like net worth can be
negative (the Gini may then exceed 1 — documented, not an error), and a source
with a negative mean (enter debt as `-debt`) has negative `S_k`, `G_k`, `R_k`
whose product is still its correct negative contribution. No absolute values
are taken.

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

## Stored results

| | |
|---|---|
| `r(gini)` | weighted Gini |
| `r(N)`, `r(mean)` | observations used, weighted mean |
| `r(decomp)` | K×5 matrix: `contrib share Sk Gk Rk` (with `source()`) |
| `r(sumdev)` | sum of contributions − Gini, identity check (with `source()`) |
| `r(sources)` | the source varlist (with `source()`) |

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
> Lerman–Yitzhaki source decomposition* (Version 1.0.0) [Stata command].
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
