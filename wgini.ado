*! wgini 1.0.0  Weighted Gini + Lerman-Yitzhaki source decomposition  2026-07-16
*! Author: ChangHwan Kim
*! License: PolyForm Noncommercial 1.0.0 (see LICENSE.md; noncommercial use only)
*! Required Notice: Copyright ChangHwan Kim (https://github.com/kchyhj/wgini)
*! Drafted with the assistance of Claude (Anthropic) under the author's direction;
*! design and validation by the author. See README.md for the full disclosure.
*----------------------------------------------------------------------*
* Weighted Gini via the Lerman-Yitzhaki covariance formula:
*     G = 2 * cov_w(x, F) / mean_w(x)
* where F is the weighted fractional rank
*     p_i = (running weight - 0.5*w_i) / W ,   W = total weight.
* Unit-free (Gini is a ratio). Applies the supplied weight.
*
* NOTE: with net worth (assets minus debt) negative values are possible,
*       so the Gini can exceed 1 (documented caveat, not an error).
*
* TIES. Tied values all receive the mid-rank of their tie group,
*     p = (cum. weight before the group + 0.5 * group weight) / W.
* Without this, tied observations get ranks in whatever order sort leaves
* them, which is random in Stata: the total Gini is provably invariant to
* that order, but source contributions and observation-level contributions
* are not, so results would change across runs. Mid-rank makes the
* decomposition well-defined and reproducible; the total Gini is unchanged.
*
*----------------------------------------------------------------------*
* SOURCE DECOMPOSITION  -- source(varlist)
*
* When the total x is the sum of K sources, x = sum_k y_k, the Gini is
* additive in the sources by linearity of covariance:
*
*     G = sum_k  2*cov_w(y_k, F) / mu        [F and mu are those of x]
*
* Each contribution factors into the familiar three terms
*
*     contrib_k = S_k * G_k * R_k
*        S_k = mu_k / mu                     share of the total
*        G_k = 2*cov_w(y_k, F_k) / mu_k      Gini of the source itself
*        R_k = cov_w(y_k, F) / cov_w(y_k,F_k)  Gini correlation with x's rank
*
* The identity holds for any signs: a source with a negative mean (e.g. debt
* entered as -debt) yields negative S_k, G_k and R_k whose product is still
* the correct (negative) contribution. No absolute values are taken.
*
* WHY THIS EXISTS. descogini (Lopez-Feldman 2006, Stata Journal 6(1):106-111)
* performs the same decomposition but refuses weights ("weights not allowed",
* rc 101), so it cannot be used with survey data. ineqrbd (Jenkins) is a
* different, regression-based (Fields/Shorrocks) decomposition. Validated
* against descogini: run with weights all equal to 1, wgini's Share column
* reproduces descogini to four decimals on SHFLC 2025 (six sources).
*
*----------------------------------------------------------------------*
* Returns:  r(gini)   weighted Gini
*           r(N)      number of observations used
*           r(mean)   weighted mean of the variable
*         with source():
*           r(decomp) K x 4 matrix, rows = sources, columns:
*                     contrib  share(=contrib/G)  Sk  Gk  Rk
*           r(sumdev) sum of contributions minus r(gini)  [identity check]
*           r(sources) the source varlist
*
* Options:
*   source(varlist)  decompose the Gini by these sources. Their sum must
*                    equal the total variable; checked to within tol().
*   tol(#)           tolerance for that check, in the units of the data.
*                    Default 1 (i.e. rounding at the last unit). Set larger
*                    for coarse data, 0 to skip the check.
*   gi(newvar)       store the observation-level Gini contribution
*                        g_i = w_i * 2*(x_i-mu)*(p_i-0.5) / (mu*W)
*                    which sums to r(gini). Useful for asking what a single
*                    household contributes, or what the Gini becomes without
*                    a given set of units.
*   noprint          suppress the display table
*
* Examples:
*   wgini networth [aw=w]
*   wgini networth [aw=w], source(s1 s2 s3 s4 s5 s6)
*   wgini networth [aw=w], gi(g_contrib)
*   statsby gini=r(gini) n=r(N), by(year agegrp): wgini networth [aw=weight]
*----------------------------------------------------------------------*
program define wgini, rclass
    version 14
    syntax varname(numeric) [if] [in] [aweight pweight fweight] ///
        [, SOURCE(varlist numeric) TOL(real 1) GI(name) noPRINT]

    marksample touse
    markout `touse' `varlist'
    if "`source'" != "" markout `touse' `source'

    * the program sorts internally; put the data back afterwards
    tempvar order
    quietly gen long `order' = _n

    tempvar w x cw p wg hi term
    tempname W mu G n
    quietly {
        if "`weight'" == "" {
            gen double `w' = 1        if `touse'
        }
        else {
            gen double `w' `exp'      if `touse'
        }
        replace `touse' = 0 if missing(`w') | `w' <= 0

        count if `touse'
        if r(N) == 0 {
            di as error "no observations"
            exit 2000
        }

        gen double `x' = `varlist' if `touse'
        * stable: deterministic tie order -> bit-reproducible sums
        sort `touse' `x' `order'

        gen double `cw' = sum(`w') if `touse'
        summ `cw' if `touse', meanonly
        scalar `W' = r(max)

        * weighted fractional mid-rank of the TOTAL. Ties all get the rank of
        * the middle of their tie group, (cum. weight before + 0.5*group
        * weight)/W, so the result does not depend on the arbitrary order in
        * which sort leaves tied observations. The weighted mean of p is
        * exactly 0.5 (sum w_i*p_i = W^2/2), so (p - 0.5) is (p - pbar).
        by `touse' `x': gen double `wg' = sum(`w') if `touse'
        by `touse' `x': gen double `hi' = `cw'[_N]  if `touse'
        by `touse' `x': replace    `wg' = `wg'[_N]  if `touse'
        gen double `p' = (`hi' - 0.5*`wg') / `W' if `touse'

        summ `x' [aw=`w'] if `touse', meanonly
        scalar `mu' = r(mean)

        gen double `term' = `w' * (`x' - `mu') * (`p' - 0.5) if `touse'
        summ `term' if `touse'
        scalar `G' = 2 * (r(sum) / `W') / `mu'

        count if `touse'
        scalar `n' = r(N)
    }

    return scalar gini = `G'
    return scalar N    = `n'
    return scalar mean = `mu'

    *------------------------------------------------------------------*
    * observation-level contribution (sums to G)
    *------------------------------------------------------------------*
    if "`gi'" != "" {
        confirm new variable `gi'
        quietly gen double `gi' = `w' * 2*(`x' - `mu')*(`p' - 0.5) / (`mu'*`W') ///
            if `touse'
        label var `gi' "Gini contribution of the observation"
    }

    *------------------------------------------------------------------*
    * source decomposition
    *------------------------------------------------------------------*
    if "`source'" != "" {
        local K : word count `source'

        * -- the sources must add up to the total --
        if `tol' > 0 {
            tempvar sum_ dev_
            quietly {
                gen double `sum_' = 0 if `touse'
                foreach v of local source {
                    replace `sum_' = `sum_' + `v' if `touse'
                }
                gen double `dev_' = abs(`varlist' - `sum_') if `touse'
                summ `dev_' if `touse', meanonly
            }
            if r(max) > `tol' {
                di as error "source() does not sum to `varlist'"
                di as error "  max |`varlist' - sum of sources| = " r(max) ///
                    " exceeds tol(" `tol' ")"
                exit 459
            }
        }

        tempname D
        matrix `D' = J(`K', 5, .)
        local rn ""
        local gsum = 0
        local k = 0
        foreach v of local source {
            local ++k
            local rn "`rn' `v'"

            tempvar yk cwk pk wgk hik tk tkk
            quietly {
                gen double `yk' = `v' if `touse'

                * mu_k and S_k
                summ `yk' [aw=`w'] if `touse', meanonly
                local muk = r(mean)
                local Sk  = `muk'/`mu'

                * cov(y_k, F) with the rank of the TOTAL -> the contribution
                gen double `tk' = `w' * (`yk' - `muk') * (`p' - 0.5) if `touse'
                summ `tk' if `touse'
                local covkF = r(sum)/`W'
                local con   = 2*`covkF'/`mu'

                * cov(y_k, F_k) with the source's OWN mid-rank -> G_k
                sort `touse' `yk' `order'
                gen double `cwk' = sum(`w') if `touse'
                by `touse' `yk': gen double `wgk' = sum(`w')  if `touse'
                by `touse' `yk': gen double `hik' = `cwk'[_N] if `touse'
                by `touse' `yk': replace    `wgk' = `wgk'[_N] if `touse'
                gen double `pk'  = (`hik' - 0.5*`wgk') / `W' if `touse'
                gen double `tkk' = `w' * (`yk' - `muk') * (`pk' - 0.5) if `touse'
                summ `tkk' if `touse'
                local covkFk = r(sum)/`W'

                * restore the total's ordering for the next source
                sort `touse' `x' `order'
            }
            * G_k and R_k are undefined for an all-zero source (mu_k = 0);
            * its contribution is 0 and is reported as such.
            local Gk = cond(`muk' != 0, 2*`covkFk'/`muk', .)
            local Rk = cond(`covkFk' != 0, `covkF'/`covkFk', .)

            matrix `D'[`k', 1] = `con'
            matrix `D'[`k', 2] = `con'/`G'
            matrix `D'[`k', 3] = `Sk'
            matrix `D'[`k', 4] = `Gk'
            matrix `D'[`k', 5] = `Rk'
            local gsum = `gsum' + `con'
        }
        matrix rownames `D' = `rn'
        matrix colnames `D' = contrib share Sk Gk Rk

        * copy: without it return matrix *moves* `D', leaving nothing for the
        * display block below.
        return matrix decomp = `D', copy
        return scalar sumdev = `gsum' - `G'
        return local sources "`source'"
    }

    * put the data back in the order the user had
    quietly sort `order'

    *------------------------------------------------------------------*
    * display
    *------------------------------------------------------------------*
    if "`print'" != "noprint" {
        di as txt _n "Weighted Gini (" as res "`varlist'" as txt ") = " ///
            as res %8.4f `G' as txt "   (N = " as res `n' as txt ")"

        if "`source'" != "" {
            di as txt _n "Lerman-Yitzhaki decomposition by source" _n ///
                "{hline 70}"
            di as txt %-22s "Source" %10s "contrib" %10s "share" ///
                %9s "Sk" %9s "Gk" %9s "Rk"
            di as txt "{hline 70}"
            local k = 0
            foreach v of local source {
                local ++k
                di as txt %-22s abbrev("`v'",21) ///
                    as res %10.4f `D'[`k',1] %10.4f `D'[`k',2] ///
                    %9.4f `D'[`k',3] %9.4f `D'[`k',4] %9.4f `D'[`k',5]
            }
            di as txt "{hline 70}"
            di as txt %-22s "Total" as res %10.4f `gsum' %10.4f `gsum'/`G'
            di as txt "{hline 70}"
            di as txt "contrib = 2*cov(y_k,F)/mu, summing to the Gini above." _n ///
                "share = contrib/Gini.  contrib = Sk*Gk*Rk." _n ///
                "Sk = mu_k/mu.  Gk = Gini of the source.  Rk = Gini correlation" _n ///
                "with the rank of `varlist'."
            local dev = `gsum' - `G'
            di as txt "Sum of contributions - Gini = " as res %10.2e `dev'
        }
    }
end
