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

```stata
* 가중 지니
wgini networth [aw=weight]

* 셀별 계산
statsby gini=r(gini) n=r(N), by(year agegrp): wgini networth [aw=weight]

* 원천 분해: 원천의 합 = 총계 변수여야 함
wgini networth [aw=weight], source(home other_re vehicles saving negdebt)

* 관측 단위 기여 (합 = 지니)
wgini networth [aw=weight], gi(gcon) noprint
gsort -networth
list networth gcon in 1   // 최상위 가구의 기여는?
```

## 방법

지니는 Lerman–Yitzhaki (1984) 공분산 형식으로 계산합니다.

```
G = 2 cov_w(x, F) / mu
```

`F`는 가중 분위 순위, `mu`는 가중 평균입니다. `source(y1 ... yK)`에서
`x = y1 + ... + yK`이면 공분산의 선형성에 의해 정확한 가법 분해가 성립합니다.

```
G = Σ_k 2 cov_w(y_k, F) / mu = Σ_k S_k × G_k × R_k
```

`S_k`는 원천의 점유율, `G_k`는 원천 자체의 지니, `R_k`는 총계 순위와의 지니
상관입니다 (Lerman & Yitzhaki 1985; Stark, Taylor & Yitzhaki 1986). 항등식
`Σ_k 기여 = G`를 `r(sumdev)`로 돌려주며, 테스트에서 ~1e-16 수준으로
성립합니다.

음수를 그대로 허용합니다. 순자산처럼 음수가 있는 변수도 계산되고(이때 지니가
1을 넘을 수 있음 — 문서화된 성질이지 오류가 아님), 부채를 `-debt`로 넣으면
`S_k`, `G_k`, `R_k`가 모두 음수가 되면서 그 곱이 올바른 음의 기여가 됩니다.
절대값을 취하지 않습니다.

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

## 반환값

| | |
|---|---|
| `r(gini)` | 가중 지니 |
| `r(N)`, `r(mean)` | 사용 관측수, 가중 평균 |
| `r(decomp)` | K×5 행렬: `contrib share Sk Gk Rk` (`source()` 사용 시) |
| `r(sumdev)` | Σ기여 − 지니, 항등식 점검 (`source()` 사용 시) |
| `r(sources)` | 원천 varlist (`source()` 사용 시) |

## 참고문헌

- Lerman, R. I., and S. Yitzhaki. 1984. A note on the calculation and
  interpretation of the Gini index. *Economics Letters* 15: 363–368.
- Lerman, R. I., and S. Yitzhaki. 1985. Income inequality effects by income
  source: A new approach and applications to the United States. *Review of
  Economics and Statistics* 67: 151–156.
- López-Feldman, A. 2006. Decomposing inequality and obtaining marginal
  effects. *Stata Journal* 6: 106–111.
- Stark, O., J. E. Taylor, and S. Yitzhaki. 1986. Remittances and inequality.
  *Economic Journal* 96: 722–740.

## 인용

`wgini`를 출판물에 사용하시면 다음과 같이 인용해 주세요
([CITATION.cff](CITATION.cff) 참조; GitHub의 "Cite this repository" 버튼이
BibTeX·APA 형식을 제공합니다):

> Kim, ChangHwan. 2026. *wgini: Weighted Gini coefficient with
> Lerman–Yitzhaki source decomposition* (Version 1.0.0) [Stata command].
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
