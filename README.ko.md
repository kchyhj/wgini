# wgini — 가중 지니계수와 Lerman–Yitzhaki 원천 분해

`wgini`는 표본 가중치 아래에서 지니계수를 계산하고, 옵션으로 지니를 소득·자산
**원천별**(Lerman & Yitzhaki 1985)로, 또는 **관측치별**로 가법 분해하는 Stata
명령어입니다. 기존 도구들이 이 조합을 지원하지 않아 만들었습니다:

| | 가중치 | LY 원천 분해 | 관측 단위 기여 |
|---|---|---|---|
| `descogini` (López-Feldman 2006, *Stata Journal*) | **불가** (`weights not allowed`) | 가능 | 불가 |
| `ineqrbd` (Jenkins) | 가능 | 불가 — 회귀 기반(Fields/Shorrocks) CV 분해로 다른 방법 | 불가 |
| `wgini` | 가능 | 가능 | 가능 |

조사 마이크로데이터(가계금융복지조사, CPS, LIS, EU-SILC 등)는 가중치가
필수라서, 가중치를 받는 원천 분해 명령이 필요했습니다.

## 설치

```stata
net install wgini, from("https://raw.githubusercontent.com/kchyhj/wgini/main/")
```

또는 `wgini.ado`와 `wgini.sthlp`를 개인 ado 폴더(`adopath`로 확인)에 복사하면
됩니다. Stata 14 이상. 설치 후 `help wgini`.

## 빠른 사용법

**1. 한 변수의 가중 지니.** `networth`의 지니를 표본 가중치를 적용해
계산하고, 이후에 쓸 수 있도록 `r(gini)`, `r(N)`, `r(mean)`을 남깁니다.

```stata
wgini networth [aw=weight]
```

**2. 원천 분해.** 순자산이 자산 구성요소의 합일 때(여기서는 거주주택 +
기타 부동산 + 자동차 + 저축 + 부채를 *음수 변수*로 넣은 것, 예:
`gen negdebt = -debt`), 지니를 구성요소별 가법 기여로 가르고 각 기여를
점유율 × 원천 내 지니 × 지니 상관으로 인수분해합니다. 이 분해는 정확한
항등식입니다 — 이유는 [방법](#방법) 참조.

```stata
wgini networth [aw=weight], source(home other_re vehicles saving negdebt)
matrix list r(decomp)     // contrib share Sk Gk Rk, 원천당 한 행
```

**3. 상위 몫 진단 — 누가 지니를 만드는가?** *측정된 불평등의 몇 %가 상위
1%에게서 오는가? 그들을 빼면 지니가 얼마가 되는가?* 같은 질문에 옵션
하나로 답합니다. 상위 몫 목록을 자유롭게 줄 수 있습니다:

```stata
wgini networth [aw=weight], top(1 5 10)
matrix list r(top)
```

각 `p`마다 네 개의 수치가 한 행으로 나옵니다:

| 열 | 의미 |
|---|---|
| `actual_pct` | 상위 집단의 실제 가중 인구 비율. 상위 집단은 가중 분위 mid-rank가 1−p/100을 넘는 관측치 전체(순위 기준 절단)인데, 가중 자료는 이산적이고 하나의 mid-rank를 공유하는 동점 그룹이 통째로 한쪽에 남으므로 `p`와 약간 다를 수 있습니다. |
| `value_share` | 상위 집단이 보유한 총액(자산 또는 소득) 점유율. |
| `gini_share` | 상위 집단이 지니에서 만드는 몫. 각 관측치는 합이 정확히 지니가 되는 가법 기여 $g_i$를 갖는데, 이 열은 상위 집단의 $g_i$를 합해 $G$로 나눈 값입니다. |
| `gini_excl` | 상위 집단을 *뺀* 표본에서 그 표본의 평균과 순위로 *재계산한* 지니. |

`value_share`와 `gini_share`를 나란히 읽으면 상위가 얼마나 불비례적인지
보입니다: 실제 전국 자산조사 응용에서 20대 가구주 집단의 최상위 한 가구가
집단 자산의 12%를 보유하며 집단 지니의 17%를 만들고 있었습니다.
`gini_excl`을 전체 지니와 비교하면 측정된 불평등이 얇은 상위층에 얼마나
얹혀 있는지 보입니다.

유의점 둘. `gini_excl`은 "지니 − `gini_share`"가 **아닙니다** — 관측치를
지우면 평균과 모든 순위가 다시 정해지므로 진짜 재계산입니다. 그리고
기여는 *양쪽* 꼬리에서 모두 양수입니다 — 부유층은 $(x_i-\mu)>0$이고
$(F_i-\tfrac12)>0$, 빈곤층은 두 인자가 모두 음수라 곱이 양수 — 그래서
*하위* 1%도 지니에 양의 기여를 하며, `top()`은 정의상 상위만 봅니다.

**4. 상위 몫이 아닌 임의 집합.** 관심 집단이 상위 몫이 아닐 때 — 특정 한
가구, 특정 지역, 특정 직업 — 는 `gi()`가 각 관측치의 기여 $g_i$를 변수로
저장하며, 어떤 집합이든 그 합을 `r(gini)`로 나누면 그 집합의 지니 몫이
됩니다. 최상위 한 가구라면:

```stata
wgini networth [aw=weight], gi(gcon) noprint
scalar G = r(gini)
gsort -networth              // 내림차순 정렬: 관측치 1번이 최상위
display gcon[1]/G            // 그 가구의 지니 몫
```

임의 집합을 *제외*하려면 `if`로 표본을 제한하면 됩니다
(예: `wgini networth if hhid != "<최상위 id>" [aw=weight]`).

## 하위 그룹별로 지니 계산하기

`wgini`는 주어진 표본에 대해 한 벌의 결과를 냅니다. 하위 그룹마다(예:
연도 × 연령대) 반복하려면 Stata의 `statsby`를 씁니다 — 어떤 명령과도
결합되는 Stata의 범용 장치이지 **`wgini`에 프로그램된 기능이 아닙니다**.
`statsby`는 그룹마다 명령을 한 번씩 실행하고 반환된 스칼라를 모아 그룹당
한 행짜리 새 데이터셋을 만듭니다:

```stata
statsby gini=r(gini) n=r(N), by(year agegrp) clear: ///
    wgini networth [aw=weight]
list year agegrp gini n
```

`top()`의 결과를 스칼라로도 반환하는 것은 바로 `statsby`가 그룹별로
수집할 수 있게 하기 위해서입니다 — 예컨대 연도별 상위 1%의 지니 몫과
상위 1% 제외 지니:

```stata
statsby gini=r(gini) gsh=r(gshare_1) gex=r(gexcl_1), by(year) clear: ///
    wgini networth [aw=weight], top(1)
```

(`statsby ..., clear`는 메모리의 데이터를 수집 결과로 바꿔치기하므로
원자료를 먼저 저장해 두세요.)

> **주의.** 이것은 그룹 안에서 지니를 따로따로 *계산*하는 것이지, 전체
> 지니를 그룹들로 *분해*하는 것이 아닙니다. 지니는 소득·자산 원천으로는
> 정확히 분해되지만(예시 2), 인구 하위 그룹으로는 분해되지 않습니다 —
> 그룹 분포가 겹치면 "집단 내 + 집단 간"의 합이 전체와 달라지고 겹침
> 잔차항이 남습니다. 집단 내/집단 간의 정확한 분해가 필요하면 일반화
> 엔트로피 지수(예: Theil)를 쓰세요. Stata에서는 `ineqdeco`(Jenkins,
> SSC)가 있습니다.

## 방법

지니는 Lerman–Yitzhaki (1984) 공분산 형식으로 계산합니다.

$$G \;=\; \frac{2\,\mathrm{cov}_w\!\big(x,\,F(x)\big)}{\mu}$$

$F(x)$는 가중 분위 순위, $\mu$는 가중 평균입니다.

**원천 분해가 정확한 이유.** 원천은 공분산의 *첫째* 인자로 들어갑니다.
$x = \sum_k y_k$이면, 공분산이 첫째 인자에 선형이고 *총계*의 순위
$F(x)$와 평균 $\mu$는 그대로 고정되므로

$$\mathrm{cov}_w\!\Big(\textstyle\sum_k y_k,\,F(x)\Big)
= \sum_k \mathrm{cov}_w\big(y_k,\,F(x)\big)
\quad\Longrightarrow\quad
G = \sum_k \frac{2\,\mathrm{cov}_w\big(y_k,\,F(x)\big)}{\mu}$$

가 성립합니다 — 원천당 한 항씩의 정확한 항등식입니다. 각 항은

$$\frac{2\,\mathrm{cov}_w\big(y_k,\,F(x)\big)}{\mu} = S_k\,G_k\,R_k$$

로 인수분해되며, $S_k$는 원천의 점유율, $G_k$는 원천 자체의 지니,
$R_k$는 총계 순위와의 지니 상관입니다 (Lerman & Yitzhaki 1985; Stark,
Taylor & Yitzhaki 1986). 항등식 $\sum_k \text{기여} = G$를 `r(sumdev)`로
돌려주며, 테스트에서 ~1e-16 수준으로 성립합니다. 기여가 큰 것은 점유율이
커서일 수도, 원천이 집중되어서($G_k$가 커서)일 수도, 원천이 전체 순위를
따라가서($R_k$가 커서)일 수도 있습니다 — 세 인자가 이를 갈라 보여줍니다.

**음수를 그대로 허용합니다 — Lerman–Yitzhaki 형식의 성질입니다.** 교과서적
지니 구성(로렌츠 곡선, 또는 절대차 평균을 평균의 2배로 나누는 방식)은
암묵적으로 음수가 없다고 가정합니다. 공분산 형식은 그렇지 않습니다: 필요한
것은 각 관측치의 평균 편차 $x_i - \mu$와 분위 순위 $F(x_i)$뿐인데, 둘 다
음수 값에서도 완벽하게 정의됩니다. 유일한 요건은 평균이 양수($\mu > 0$)라는
것입니다. 따라서 부채 초과 가구가 음수로 나타나는 순자산 같은 변수도 0에서
자르거나 음수 관측치를 버리지 않고 그대로 계산합니다. 한 가지 유의할 결과:
음수가 있으면 지니가 더 이상 1에 묶이지 않고 1을 넘을 수 있습니다(문서화된
성질이지 오류가 아님). 마찬가지로 평균이 음수인 원천(부채를 `-debt`로
입력)은 $S_k$, $G_k$, $R_k$가 모두 음수가 되면서 그 곱이 올바른 음의 기여가
됩니다. 어디에서도 절대값을 취하지 않습니다.

### 동점 처리와 재현성

동점 값에는 동점 그룹의 **mid-rank**(그룹 직전 누적 가중 + 그룹 가중의
절반)를 일괄 배정합니다. 이것이 생각보다 중요합니다: Stata의 `sort`는 동점
관측치의 순서를 실행마다 임의로 바꾸는데, 총지니는 그 순서에 수학적으로
불변이지만 **원천별 기여와 관측 단위 기여는 불변이 아니어서** 임의 순위를
쓰면 실행할 때마다 값이 달라집니다(실제 자산조사 자료의 순자산 0 동점이 많은
소집단에서 최대 ~1e-4 수준). mid-rank는 분해를 유일하게 정의하고 완전한
재현성을 줍니다 — 같은 자료를 넣으면 bit 단위로 같은 결과가 나옵니다.
실행 후 데이터의 정렬 순서도 원래대로 복원합니다.

## 검증

테스트(`test/wgini_test.do`, 인공 자료만 사용하는 자체완결 스크립트) 항목:

1. **항등식**: 원천 기여의 합 = 지니 (~1e-16).
2. **인수분해**: `기여_k = S_k × G_k × R_k`가 행마다 성립. 음의 평균을 갖는
   부채 원천 포함.
3. **`descogini` 대조**: 가중치를 모두 1로 두면 — `descogini`가 받는 유일한
   경우 — `wgini`의 `Sk`, `Gk`, `Rk`, share가 `descogini`와 소수점 넷째
   자리까지 일치. (전국 자산조사 18,664가구, 6원천에서 전 칸 일치를
   확인했습니다.) 실제 가중치를 쓰면 `descogini`가 가중치를 거부하므로 비교
   자체가 불가능합니다.
4. **동점 불변**: 동점 관측치의 물리적 순서를 뒤섞어도 모든 결과 불변.
5. **재현성**: 같은 자료에서 두 번 실행하면 `r(decomp)`가 bit 단위로 동일.
6. **정렬 복원**: 명령이 사용자 데이터의 정렬을 바꾸지 않음.
7. **오류 처리**: 원천의 합이 총계와 다르면 거부(r 459); 전부 0인 원천은
   기여 0, `G_k`·`R_k` 결측.
8. **`top()` = 수동 레시피**: 지니 몫이 `gi()` 합산과, 제외 지니가 `if`
   재실행과, 총액 점유가 직접 합산과 일치.
9. **`statsby` 수집**: 그룹별로 수집한 `top()` 스칼라가 그룹별 직접 호출을
   재현.

## 반환값

| | |
|---|---|
| `r(gini)` | 가중 지니 |
| `r(N)`, `r(mean)` | 사용 관측수, 가중 평균 |
| `r(decomp)` | K×5 행렬: `contrib share Sk Gk Rk` (`source()` 사용 시) |
| `r(sumdev)` | Σ기여 − 지니, 항등식 점검 (`source()` 사용 시) |
| `r(sources)` | 원천 varlist (`source()` 사용 시) |
| `r(top)` | K×5 행렬: `top_pct actual_pct value_share gini_share gini_excl` (`top()` 사용 시) |
| `r(actual_p)`, `r(vshare_p)`, `r(gshare_p)`, `r(gexcl_p)` | `top()`의 각 `p`별 스칼라 (값으로 명명, `.` → `_`) — `statsby` 수집용 |

## 참고문헌

- Lerman, R. I., and S. Yitzhaki. 1984. A note on the calculation and
  interpretation of the Gini index. *Economics Letters* 15: 363–368.
- Lerman, R. I., and S. Yitzhaki. 1985. Income inequality effects by income
  source: A new approach and applications to the United States. *Review of
  Economics and Statistics* 67: 151–156.
- Lerman, R. I., and S. Yitzhaki. 1989. Improving the accuracy of estimates of
  Gini coefficients. *Journal of Econometrics* 42: 43–47.
- López-Feldman, A. 2006. Decomposing inequality and obtaining marginal
  effects. *Stata Journal* 6: 106–111.
- Stark, O., J. E. Taylor, and S. Yitzhaki. 1986. Remittances and inequality.
  *Economic Journal* 96: 722–740.

## 인용

`wgini`를 출판물에 사용하시면 다음과 같이 인용해 주세요
([CITATION.cff](CITATION.cff) 참조; GitHub의 "Cite this repository" 버튼이
BibTeX·APA 형식을 제공합니다):

> Kim, ChangHwan. 2026. *wgini: Weighted Gini coefficient with
> Lerman–Yitzhaki source decomposition* (Version 1.1.0) [Stata command].
> https://github.com/kchyhj/wgini

## 저자

김창환 (ChangHwan Kim), University of Kansas
(ORCID [0000-0001-7149-1386](https://orcid.org/0000-0001-7149-1386))

**AI 사용 공개.** 코드와 문서 초안 작성에 Claude(Anthropic; Claude Code,
Opus 4.8)의 도움을 받았으며, 작업은 저자의 지휘 아래 이루어졌습니다. 통계적
설계, 추정량과 동점 처리 방식의 선택, 모든 검증 결정은 저자의 것입니다. 매
릴리스는 공개 전에 `test/`의 테스트 스위트와 `descogini` 대조로 검증합니다.

English documentation: [README.md](README.md).

## 라이선스

[PolyForm Noncommercial 1.0.0](LICENSE.md) — 연구·교육·정부·비영리 등 모든
**비상업** 목적의 사용·수정·배포가 자유롭습니다(재원 출처와 무관). 상업적
이용은 저자와의 별도 라이선스가 필요합니다.

Required Notice: Copyright ChangHwan Kim (https://github.com/kchyhj/wgini)
