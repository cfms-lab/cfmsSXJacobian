# 조건 함수 → 텐서 꼴 (시제품 출력 미리 보기)

## 1. 인장

### BaraffStretch — 텐서 꼴

Baraff-Witkin 1998 stretch condition C = sqrt(a) (|w_u| - 1, |w_v| - 1) on triangle (0, 1, 2)

**조건 함수** (꼭짓점 3 개, 성분 2 개)

$$C_{1} = \sqrt{a} \left(|\mathbf{w}_u| - 1\right)$$

$$C_{2} = \sqrt{a} \left(|\mathbf{w}_v| - 1\right)$$

**중간 벡터** — 좌표의 1 차식이다. $m$ 은 꼭짓점 번호(0 부터).

$$\mathbf{w}_u = \sum_m c^{u}_{m}\,\mathbf{x}_m,\qquad c^{u}_{0} = \frac{v_{1} - v_{2}}{\Delta},\quad c^{u}_{1} = \frac{- v_{0} + v_{2}}{\Delta},\quad c^{u}_{2} = \frac{v_{0} - v_{1}}{\Delta}$$

$$\mathbf{w}_v = \sum_m c^{v}_{m}\,\mathbf{x}_m,\qquad c^{v}_{0} = \frac{- u_{1} + u_{2}}{\Delta},\quad c^{v}_{1} = \frac{u_{0} - u_{2}}{\Delta},\quad c^{v}_{2} = \frac{- u_{0} + u_{1}}{\Delta}$$

여기서 $\Delta = \left(- u_{0} + u_{1}\right) \left(- v_{0} + v_{2}\right) - \left(- u_{0} + u_{2}\right) \left(- v_{0} + v_{1}\right),\quad a = \tfrac12\,|\Delta|$ 이고, 계수는 쉬는 패턴만으로 정해진다($\sum_m c_m = 0$ — 병진 불변).

**1 차 미분** — 꼭짓점 $m$ 에 대한 3-벡터

$$\frac{\partial C_{1}}{\partial \mathbf{x}_m} = \sqrt{a}\,c^{u}_{m}\,\hat{\mathbf{w}}_u$$

$$\frac{\partial C_{2}}{\partial \mathbf{x}_m} = \sqrt{a}\,c^{v}_{m}\,\hat{\mathbf{w}}_v$$

**2 차 미분** — 꼭짓점 쌍 $(m, n)$ 의 3 × 3 블록

$$\frac{\partial^2 C_{1}}{\partial \mathbf{x}_m\,\partial \mathbf{x}_n} = \frac{\sqrt{a}}{|\mathbf{w}_u|}\,c^{u}_{m}\,c^{u}_{n}\left(\mathbf{I} - \hat{\mathbf{w}}_u\hat{\mathbf{w}}_u^{\mathsf T}\right)$$

$$\frac{\partial^2 C_{2}}{\partial \mathbf{x}_m\,\partial \mathbf{x}_n} = \frac{\sqrt{a}}{|\mathbf{w}_v|}\,c^{v}_{m}\,c^{v}_{n}\left(\mathbf{I} - \hat{\mathbf{w}}_v\hat{\mathbf{w}}_v^{\mathsf T}\right)$$

**지배방정식에 들어가는 항** — 꼭짓점 블록으로. $\nabla_m C \equiv \partial C/\partial\mathbf{x}_m$, $\nabla^2_{mn} C \equiv \partial^2 C/\partial\mathbf{x}_m\partial\mathbf{x}_n$, 합은 조건 성분 $c$ 에 대해.

$$\dot C_c = \sum_n \nabla_n C_c\cdot\mathbf{v}_n$$

$$\mathbf{f}_m = -k\sum_c C_c\,\nabla_m C_c,\qquad \mathbf{d}_m = -k_d\sum_c \dot C_c\,\nabla_m C_c$$

$$\mathbf{K}_{mn} = \frac{\partial \mathbf{f}_m}{\partial \mathbf{x}_n} = -k\sum_c\Big[\nabla_m C_c\,\nabla_n C_c^{\mathsf T} + \beta_c\,\nabla^2_{mn} C_c\Big],\qquad \mathbf{V}_{mn} = \frac{\partial \mathbf{d}_m}{\partial \mathbf{v}_n} = -k_d\sum_c \nabla_m C_c\,\nabla_n C_c^{\mathsf T}$$

$$\mathbf{A} = \mathbf{M} - h\,\mathbf{V} - h^2\,\mathbf{K},\qquad \mathbf{b} = h\,(\mathbf{f} + \mathbf{d}) + h^2\,\mathbf{K}\,\mathbf{v}$$

| 강성 판 | $\beta_c$ | 생성 함수 |
|---|---|---|
| 정확 | $C_c + (k_d/k)\,\dot C_c$ | `EvalImplicit` · `EvalEngine` |
| Gauss–Newton | $0$ | `EvalImplicitGN` · `EvalEngineGN` |
| 인장 C > 0 | $\max(C_c, 0)$ | `EvalImplicitTension` · `EvalEngineTension` |
| 엔진 처방 | $\max(C_c, 0) + (k_d/k)\max(\dot C_c, 0)$ | `EvalEngineTensionRate` |


검사: {'samples': 50.0, 'max_rel_C': 2.9729135545965e-15, 'max_rel_grad': 5.666314672112641e-16, 'max_rel_hess': 6.567850685691378e-16}

## 2. Provot 전단

### ProvotShear — 텐서 꼴

Provot-type shear (diagonal stretch) C = sqrt(a) (|w_u + w_v| - b_plus, |w_u - w_v| - b_minus) on triangle (0, 1, 2)

**조건 함수** (꼭짓점 3 개, 성분 2 개)

$$C_{1} = \sqrt{a} \left(- b_{+} + |\mathbf{d}_{+}|\right)$$

$$C_{2} = \sqrt{a} \left(- b_{-} + |\mathbf{d}_{-}|\right)$$

**중간 벡터** — 좌표의 1 차식이다. $m$ 은 꼭짓점 번호(0 부터).

$$\mathbf{d}_{+} = \sum_m c^{{+}}_{m}\,\mathbf{x}_m,\qquad c^{{+}}_{0} = \frac{- u_{1} + u_{2} + v_{1} - v_{2}}{\Delta},\quad c^{{+}}_{1} = \frac{u_{0} - u_{2} - v_{0} + v_{2}}{\Delta},\quad c^{{+}}_{2} = \frac{- u_{0} + u_{1} + v_{0} - v_{1}}{\Delta}$$

$$\mathbf{d}_{-} = \sum_m c^{{-}}_{m}\,\mathbf{x}_m,\qquad c^{{-}}_{0} = \frac{u_{1} - u_{2} + v_{1} - v_{2}}{\Delta},\quad c^{{-}}_{1} = \frac{- u_{0} + u_{2} - v_{0} + v_{2}}{\Delta},\quad c^{{-}}_{2} = \frac{u_{0} - u_{1} + v_{0} - v_{1}}{\Delta}$$

여기서 $\Delta = \left(- u_{0} + u_{1}\right) \left(- v_{0} + v_{2}\right) - \left(- u_{0} + u_{2}\right) \left(- v_{0} + v_{1}\right),\quad a = \tfrac12\,|\Delta|$ 이고, 계수는 쉬는 패턴만으로 정해진다($\sum_m c_m = 0$ — 병진 불변).

**1 차 미분** — 꼭짓점 $m$ 에 대한 3-벡터

$$\frac{\partial C_{1}}{\partial \mathbf{x}_m} = \sqrt{a}\,c^{{+}}_{m}\,\hat{\mathbf{d}}_{+}$$

$$\frac{\partial C_{2}}{\partial \mathbf{x}_m} = \sqrt{a}\,c^{{-}}_{m}\,\hat{\mathbf{d}}_{-}$$

**2 차 미분** — 꼭짓점 쌍 $(m, n)$ 의 3 × 3 블록

$$\frac{\partial^2 C_{1}}{\partial \mathbf{x}_m\,\partial \mathbf{x}_n} = \frac{\sqrt{a}}{|\mathbf{d}_{+}|}\,c^{{+}}_{m}\,c^{{+}}_{n}\left(\mathbf{I} - \hat{\mathbf{d}}_{+}\hat{\mathbf{d}}_{+}^{\mathsf T}\right)$$

$$\frac{\partial^2 C_{2}}{\partial \mathbf{x}_m\,\partial \mathbf{x}_n} = \frac{\sqrt{a}}{|\mathbf{d}_{-}|}\,c^{{-}}_{m}\,c^{{-}}_{n}\left(\mathbf{I} - \hat{\mathbf{d}}_{-}\hat{\mathbf{d}}_{-}^{\mathsf T}\right)$$

**지배방정식에 들어가는 항** — 꼭짓점 블록으로. $\nabla_m C \equiv \partial C/\partial\mathbf{x}_m$, $\nabla^2_{mn} C \equiv \partial^2 C/\partial\mathbf{x}_m\partial\mathbf{x}_n$, 합은 조건 성분 $c$ 에 대해.

$$\dot C_c = \sum_n \nabla_n C_c\cdot\mathbf{v}_n$$

$$\mathbf{f}_m = -k\sum_c C_c\,\nabla_m C_c,\qquad \mathbf{d}_m = -k_d\sum_c \dot C_c\,\nabla_m C_c$$

$$\mathbf{K}_{mn} = \frac{\partial \mathbf{f}_m}{\partial \mathbf{x}_n} = -k\sum_c\Big[\nabla_m C_c\,\nabla_n C_c^{\mathsf T} + \beta_c\,\nabla^2_{mn} C_c\Big],\qquad \mathbf{V}_{mn} = \frac{\partial \mathbf{d}_m}{\partial \mathbf{v}_n} = -k_d\sum_c \nabla_m C_c\,\nabla_n C_c^{\mathsf T}$$

$$\mathbf{A} = \mathbf{M} - h\,\mathbf{V} - h^2\,\mathbf{K},\qquad \mathbf{b} = h\,(\mathbf{f} + \mathbf{d}) + h^2\,\mathbf{K}\,\mathbf{v}$$

| 강성 판 | $\beta_c$ | 생성 함수 |
|---|---|---|
| 정확 | $C_c + (k_d/k)\,\dot C_c$ | `EvalImplicit` · `EvalEngine` |
| Gauss–Newton | $0$ | `EvalImplicitGN` · `EvalEngineGN` |
| 인장 C > 0 | $\max(C_c, 0)$ | `EvalImplicitTension` · `EvalEngineTension` |
| 엔진 처방 | $\max(C_c, 0) + (k_d/k)\max(\dot C_c, 0)$ | `EvalEngineTensionRate` |


검사: {'samples': 50.0, 'max_rel_C': 4.3696602676284256e-14, 'max_rel_grad': 6.923776399048371e-16, 'max_rel_hess': 8.110498570973231e-16}

## 3. 굽힘 — 원자 수준의 전개

**원자의 정의** — 힌지 $(h, i, j, k)$

$$\mathbf{N}_1 = (\mathbf{p}_i - \mathbf{p}_h)\times(\mathbf{p}_j - \mathbf{p}_h),\qquad \mathbf{N}_2 = (\mathbf{p}_j - \mathbf{p}_k)\times(\mathbf{p}_i - \mathbf{p}_k),\qquad \tilde{e} = [\mathbf{p}_j - \mathbf{p}_i]_\times$$

$$n_1 = \mathbf{N}_1/|\mathbf{N}_1|,\qquad n_2 = \mathbf{N}_2/|\mathbf{N}_2|$$

$[\mathbf{a}]_\times$ 는 $\mathbf{a}\times$ 의 3 × 3 반대칭 행렬이다.

**원자의 1 차 미분** (3 × 3 블록)

$$\frac{\partial \mathbf{N}_1}{\partial p_h} = [\mathbf{p}_j - \mathbf{p}_i]_\times = \tilde{e},\qquad \frac{\partial \mathbf{N}_1}{\partial p_i} = [\mathbf{p}_h - \mathbf{p}_j]_\times,\qquad \frac{\partial \mathbf{N}_1}{\partial p_j} = [\mathbf{p}_i - \mathbf{p}_h]_\times$$

$$\frac{\partial \mathbf{N}_2}{\partial p_k} = [\mathbf{p}_i - \mathbf{p}_j]_\times = -\tilde{e},\qquad \frac{\partial \mathbf{N}_2}{\partial p_j} = [\mathbf{p}_k - \mathbf{p}_i]_\times,\qquad \frac{\partial \mathbf{N}_2}{\partial p_i} = [\mathbf{p}_j - \mathbf{p}_k]_\times$$

$$\frac{\partial n}{\partial p_m} = \frac{1}{|\mathbf{N}|}\left(\mathbf{I} - n\,n^{\mathsf T}\right)\frac{\partial \mathbf{N}}{\partial p_m}$$

$$\frac{\partial \tilde{e}}{\partial p_j}[\boldsymbol\delta] = [\boldsymbol\delta]_\times,\qquad \frac{\partial \tilde{e}}{\partial p_i}[\boldsymbol\delta] = -[\boldsymbol\delta]_\times$$

**원자의 2 차 미분** — 방향 $\mathbf{a}$(첫 꼭짓점) · $\mathbf{b}$(둘째 꼭짓점)에 대한 값. $\mathbf{N} = (P - H)\times(Q - H)$ 로 쓰면($\mathbf{N}_1$: $H{=}h, P{=}i, Q{=}j$; $\mathbf{N}_2$: $H{=}k, P{=}j, Q{=}i$)

$$D^2\mathbf{N}[\mathbf{a}, \mathbf{b}] = \sigma\;\mathbf{a}\times\mathbf{b},\qquad \sigma_{PQ} = \sigma_{HP} = \sigma_{QH} = +1,\quad \sigma_{QP} = \sigma_{PH} = \sigma_{HQ} = -1,\quad \sigma_{PP} = \sigma_{QQ} = \sigma_{HH} = 0$$

$$D^2 n[\mathbf{a}, \mathbf{b}] = \frac{\mathbf{P}\,D^2\mathbf{N}[\mathbf{a}, \mathbf{b}]}{|\mathbf{N}|} - \frac{(n\cdot \mathbf{N}_b)\,\mathbf{P}\mathbf{N}_a + (n\cdot \mathbf{N}_a)\,\mathbf{P}\mathbf{N}_b + (\mathbf{N}_a\cdot\mathbf{P}\mathbf{N}_b)\,n}{|\mathbf{N}|^2},\qquad \mathbf{P} = \mathbf{I} - n\,n^{\mathsf T},\quad \mathbf{N}_a = D\mathbf{N}[\mathbf{a}]$$

$\tilde{e}$ 는 좌표의 1 차식이라 2 차 미분이 0 이다.

**단위 변** $\hat{\mathbf{e}} = \mathbf{e}/|\mathbf{e}|$, $\tilde{\hat{e}} = [\hat{\mathbf{e}}]_\times$ — 단위 법선과 같은 꼴이다($\mathbf{e}_a = D\mathbf{e}[\mathbf{a}]$: 꼭짓점 $j$ 이면 $\mathbf{a}$, $i$ 이면 $-\mathbf{a}$).

$$D\hat{\mathbf{e}}[\mathbf{a}] = \frac{\mathbf{P}_e\,\mathbf{e}_a}{|\mathbf{e}|},\qquad D^2\hat{\mathbf{e}}[\mathbf{a}, \mathbf{b}] = -\frac{(\hat{\mathbf{e}}\cdot\mathbf{e}_b)\,\mathbf{P}_e\mathbf{e}_a + (\hat{\mathbf{e}}\cdot\mathbf{e}_a)\,\mathbf{P}_e\mathbf{e}_b + (\mathbf{e}_a\cdot\mathbf{P}_e\mathbf{e}_b)\,\hat{\mathbf{e}}}{|\mathbf{e}|^2},\qquad \mathbf{P}_e = \mathbf{I} - \hat{\mathbf{e}}\hat{\mathbf{e}}^{\mathsf T}$$


$$C_{\mathrm{bend}} = n_1 \cdot \tilde{e} \cdot \left(n_1 + n_2\right)$$

**1 차 미분**

$$\frac{\partial C_{\mathrm{bend}}}{\partial p_{h}} = n_1 \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}} + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot n_1 + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot n_2$$

$$\frac{\partial C_{\mathrm{bend}}}{\partial p_{i}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1 + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + n_1 \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot n_1 + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot n_2$$

$$\frac{\partial C_{\mathrm{bend}}}{\partial p_{j}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + n_1 \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot n_2$$

$$\frac{\partial C_{\mathrm{bend}}}{\partial p_{k}} = n_1 \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

**2 차 미분** — 꼭짓점 쌍마다(행 = 첫 미분의 꼭짓점)

(h, ·)

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{h}\,\partial p_{h}} = n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{h}} + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{h}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{h}} \cdot \tilde{e} \cdot n_2 + 2 \cdot \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{h}\,\partial p_{i}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_1}{\partial p_{h}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{h}\,\partial p_{j}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_1}{\partial p_{h}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{h}\,\partial p_{k}} = \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(i, ·)

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{i}\,\partial p_{h}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_1}{\partial p_{h}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{i}\,\partial p_{i}} = n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{i}} \cdot n_1 + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{i}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{i}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{i}} + 2 \cdot n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_1}{\partial p_{i}} + 2 \cdot n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{i}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{i}} \cdot \tilde{e} \cdot n_2 + 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1 + 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{i}\,\partial p_{j}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_1}{\partial p_{j}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{j}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_1}{\partial p_{i}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{i}} + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{j}} \cdot n_1 + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{j}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{j}} + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2 + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2 + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{i}\,\partial p_{k}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{k}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(j, ·)

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{j}\,\partial p_{h}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_1}{\partial p_{h}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2 + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{h}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{j}\,\partial p_{i}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_1}{\partial p_{j}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{j}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_1}{\partial p_{i}} + n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{i}} + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{j}} \cdot n_1 + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{i}\,\partial p_{j}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{j}} + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_2 + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot \tilde{e} \cdot n_2 + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2 + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot n_1$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{j}\,\partial p_{j}} = n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{j}\,\partial p_{j}} \cdot n_1 + n_1 \cdot \frac{\partial^{2} \tilde{e}}{\partial p_{j}\,\partial p_{j}} \cdot n_2 + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_1}{\partial p_{j}\,\partial p_{j}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{j}} + 2 \cdot n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_1}{\partial p_{j}} + 2 \cdot n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial^{2} n_1}{\partial p_{j}\,\partial p_{j}} \cdot \tilde{e} \cdot n_1 + \frac{\partial^{2} n_1}{\partial p_{j}\,\partial p_{j}} \cdot \tilde{e} \cdot n_2 + 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_1}{\partial p_{j}} + 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{j}} + 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_1 + 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot n_2$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{j}\,\partial p_{k}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{k}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(k, ·)

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{k}\,\partial p_{h}} = \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{k}\,\partial p_{i}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{k}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{k}\,\partial p_{j}} = n_1 \cdot \frac{\partial \tilde{e}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{k}} + n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{e} \cdot \frac{\partial n_2}{\partial p_{k}}$$

$$\frac{\partial^2 C_{\mathrm{bend}}}{\partial p_{k}\,\partial p_{k}} = n_1 \cdot \tilde{e} \cdot \frac{\partial^{2} n_2}{\partial p_{k}\,\partial p_{k}}$$


검사(비정규 법선): {'samples': 5.0, 'max_rel_C': 1.3867724806200045e-16, 'max_rel_grad': 2.8424026979849545e-16, 'max_rel_hess': 3.791230428185586e-16}

검사(단위 법선): {'samples': 3.0, 'max_rel_C': 1.6789260308285602e-16, 'max_rel_grad': 4.556924105470135e-16, 'max_rel_hess': 7.648764501444884e-16}

## 4. 굽힘 명세 (atan2)

### BaraffBend — 텐서 꼴

Baraff-Witkin 1998 bend condition C = theta - theta0 on hinge (h, i, j, k)

$$C = F(\theta) = - \theta_0 + \theta,\qquad F^\prime = 1,\qquad F^{\prime\prime} = 0$$

$\theta$ 는 두 삼각형의 법선이 이루는 각이다. 단위 원자로 쓰면($n_1 = \hat{\mathbf{n}}_A$, $n_2 = \hat{\mathbf{n}}_B$; $(\mathbf{a}\times\mathbf{b})\cdot\hat{\mathbf{e}} = -\mathbf{a}\cdot(\hat{\mathbf{e}}\times\mathbf{b})$)

$$s \equiv \sin\theta = - n_1 \cdot \tilde{\hat{e}} \cdot n_2,\qquad c \equiv \cos\theta = n_1 \cdot n_2,\qquad \theta = \operatorname{atan2}(s, c)$$

명세가 정규화하지 않은 법선 · 변으로 atan2 를 쓰더라도 atan2 는 양의 배율에 불변이라 같은 $\theta$ 다(`check()` 가 확인한다).

**연쇄 법칙** — $\nabla_m \equiv \partial/\partial p_m$, $\nabla^2_{mn} \equiv \partial^2/\partial p_m\partial p_n$. 단위 벡터에서는 $s^2 + c^2 = 1$ 이 항등식이다.

$$\nabla_m \theta = c\,\nabla_m s - s\,\nabla_m c$$

$$\nabla^2_{mn} \theta = c\,\nabla^2_{mn} s - s\,\nabla^2_{mn} c - 2sc\left(\nabla_m s\,\nabla_n s^{\mathsf T} - \nabla_m c\,\nabla_n c^{\mathsf T}\right) + (s^2 - c^2)\left(\nabla_m s\,\nabla_n c^{\mathsf T} + \nabla_m c\,\nabla_n s^{\mathsf T}\right)$$

$$\nabla_m C = F^\prime\,\nabla_m\theta,\qquad \nabla^2_{mn} C = F^\prime\,\nabla^2_{mn}\theta + F^{\prime\prime}\,\nabla_m\theta\,\nabla_n\theta^{\mathsf T}$$

**원자의 정의** — 힌지 $(h, i, j, k)$

$$\mathbf{N}_1 = (\mathbf{p}_i - \mathbf{p}_h)\times(\mathbf{p}_j - \mathbf{p}_h),\qquad \mathbf{N}_2 = (\mathbf{p}_j - \mathbf{p}_k)\times(\mathbf{p}_i - \mathbf{p}_k),\qquad \tilde{e} = [\mathbf{p}_j - \mathbf{p}_i]_\times$$

$$n_1 = \mathbf{N}_1/|\mathbf{N}_1|,\qquad n_2 = \mathbf{N}_2/|\mathbf{N}_2|$$

$[\mathbf{a}]_\times$ 는 $\mathbf{a}\times$ 의 3 × 3 반대칭 행렬이다.

**원자의 1 차 미분** (3 × 3 블록)

$$\frac{\partial \mathbf{N}_1}{\partial p_h} = [\mathbf{p}_j - \mathbf{p}_i]_\times = \tilde{e},\qquad \frac{\partial \mathbf{N}_1}{\partial p_i} = [\mathbf{p}_h - \mathbf{p}_j]_\times,\qquad \frac{\partial \mathbf{N}_1}{\partial p_j} = [\mathbf{p}_i - \mathbf{p}_h]_\times$$

$$\frac{\partial \mathbf{N}_2}{\partial p_k} = [\mathbf{p}_i - \mathbf{p}_j]_\times = -\tilde{e},\qquad \frac{\partial \mathbf{N}_2}{\partial p_j} = [\mathbf{p}_k - \mathbf{p}_i]_\times,\qquad \frac{\partial \mathbf{N}_2}{\partial p_i} = [\mathbf{p}_j - \mathbf{p}_k]_\times$$

$$\frac{\partial n}{\partial p_m} = \frac{1}{|\mathbf{N}|}\left(\mathbf{I} - n\,n^{\mathsf T}\right)\frac{\partial \mathbf{N}}{\partial p_m}$$

$$\frac{\partial \tilde{e}}{\partial p_j}[\boldsymbol\delta] = [\boldsymbol\delta]_\times,\qquad \frac{\partial \tilde{e}}{\partial p_i}[\boldsymbol\delta] = -[\boldsymbol\delta]_\times$$

**원자의 2 차 미분** — 방향 $\mathbf{a}$(첫 꼭짓점) · $\mathbf{b}$(둘째 꼭짓점)에 대한 값. $\mathbf{N} = (P - H)\times(Q - H)$ 로 쓰면($\mathbf{N}_1$: $H{=}h, P{=}i, Q{=}j$; $\mathbf{N}_2$: $H{=}k, P{=}j, Q{=}i$)

$$D^2\mathbf{N}[\mathbf{a}, \mathbf{b}] = \sigma\;\mathbf{a}\times\mathbf{b},\qquad \sigma_{PQ} = \sigma_{HP} = \sigma_{QH} = +1,\quad \sigma_{QP} = \sigma_{PH} = \sigma_{HQ} = -1,\quad \sigma_{PP} = \sigma_{QQ} = \sigma_{HH} = 0$$

$$D^2 n[\mathbf{a}, \mathbf{b}] = \frac{\mathbf{P}\,D^2\mathbf{N}[\mathbf{a}, \mathbf{b}]}{|\mathbf{N}|} - \frac{(n\cdot \mathbf{N}_b)\,\mathbf{P}\mathbf{N}_a + (n\cdot \mathbf{N}_a)\,\mathbf{P}\mathbf{N}_b + (\mathbf{N}_a\cdot\mathbf{P}\mathbf{N}_b)\,n}{|\mathbf{N}|^2},\qquad \mathbf{P} = \mathbf{I} - n\,n^{\mathsf T},\quad \mathbf{N}_a = D\mathbf{N}[\mathbf{a}]$$

$\tilde{e}$ 는 좌표의 1 차식이라 2 차 미분이 0 이다.

**단위 변** $\hat{\mathbf{e}} = \mathbf{e}/|\mathbf{e}|$, $\tilde{\hat{e}} = [\hat{\mathbf{e}}]_\times$ — 단위 법선과 같은 꼴이다($\mathbf{e}_a = D\mathbf{e}[\mathbf{a}]$: 꼭짓점 $j$ 이면 $\mathbf{a}$, $i$ 이면 $-\mathbf{a}$).

$$D\hat{\mathbf{e}}[\mathbf{a}] = \frac{\mathbf{P}_e\,\mathbf{e}_a}{|\mathbf{e}|},\qquad D^2\hat{\mathbf{e}}[\mathbf{a}, \mathbf{b}] = -\frac{(\hat{\mathbf{e}}\cdot\mathbf{e}_b)\,\mathbf{P}_e\mathbf{e}_a + (\hat{\mathbf{e}}\cdot\mathbf{e}_a)\,\mathbf{P}_e\mathbf{e}_b + (\mathbf{e}_a\cdot\mathbf{P}_e\mathbf{e}_b)\,\hat{\mathbf{e}}}{|\mathbf{e}|^2},\qquad \mathbf{P}_e = \mathbf{I} - \hat{\mathbf{e}}\hat{\mathbf{e}}^{\mathsf T}$$


#### $s$ 의 전개

$$s = - n_1 \cdot \tilde{\hat{e}} \cdot n_2$$

**1 차 미분**

$$\frac{\partial s}{\partial p_{h}} = - \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial s}{\partial p_{i}} = - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot n_2 - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{i}} - \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial s}{\partial p_{j}} = - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot n_2 - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{j}} - \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial s}{\partial p_{k}} = - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

**2 차 미분** — 꼭짓점 쌍마다(행 = 첫 미분의 꼭짓점)

(h, ·)

$$\frac{\partial^2 s}{\partial p_{h}\,\partial p_{h}} = - \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{h}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial^2 s}{\partial p_{h}\,\partial p_{i}} = - \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot n_2 - \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{i}} - \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial^2 s}{\partial p_{h}\,\partial p_{j}} = - \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot n_2 - \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{j}} - \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial^2 s}{\partial p_{h}\,\partial p_{k}} = - \frac{\partial n_1}{\partial p_{h}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(i, ·)

$$\frac{\partial^2 s}{\partial p_{i}\,\partial p_{i}} = - n_1 \cdot \frac{\partial^{2} \tilde{\hat{e}}}{\partial p_{i}\,\partial p_{i}} \cdot n_2 - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{i}} - 2 \cdot n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{i}} - \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{i}} \cdot \tilde{\hat{e}} \cdot n_2 - 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot n_2 - 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{i}}$$

$$\frac{\partial^2 s}{\partial p_{i}\,\partial p_{j}} = - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{i}} - n_1 \cdot \frac{\partial^{2} \tilde{\hat{e}}}{\partial p_{i}\,\partial p_{j}} \cdot n_2 - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{j}} - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{j}} - \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot n_2 - \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{j}} - \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot n_2 - \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{i}} - \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot \tilde{\hat{e}} \cdot n_2$$

$$\frac{\partial^2 s}{\partial p_{i}\,\partial p_{k}} = - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{k}} - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{k}} - \frac{\partial n_1}{\partial p_{i}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(j, ·)

$$\frac{\partial^2 s}{\partial p_{j}\,\partial p_{j}} = - n_1 \cdot \frac{\partial^{2} \tilde{\hat{e}}}{\partial p_{j}\,\partial p_{j}} \cdot n_2 - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{j}} - 2 \cdot n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{j}} - \frac{\partial^{2} n_1}{\partial p_{j}\,\partial p_{j}} \cdot \tilde{\hat{e}} \cdot n_2 - 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot n_2 - 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{j}}$$

$$\frac{\partial^2 s}{\partial p_{j}\,\partial p_{k}} = - n_1 \cdot \frac{\partial \tilde{\hat{e}}}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{k}} - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{k}} - \frac{\partial n_1}{\partial p_{j}} \cdot \tilde{\hat{e}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(k, ·)

$$\frac{\partial^2 s}{\partial p_{k}\,\partial p_{k}} = - n_1 \cdot \tilde{\hat{e}} \cdot \frac{\partial^{2} n_2}{\partial p_{k}\,\partial p_{k}}$$


#### $c$ 의 전개

$$c = n_1 \cdot n_2$$

**1 차 미분**

$$\frac{\partial c}{\partial p_{h}} = \frac{\partial n_1}{\partial p_{h}} \cdot n_2$$

$$\frac{\partial c}{\partial p_{i}} = n_1 \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial n_1}{\partial p_{i}} \cdot n_2$$

$$\frac{\partial c}{\partial p_{j}} = n_1 \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot n_2$$

$$\frac{\partial c}{\partial p_{k}} = n_1 \cdot \frac{\partial n_2}{\partial p_{k}}$$

**2 차 미분** — 꼭짓점 쌍마다(행 = 첫 미분의 꼭짓점)

(h, ·)

$$\frac{\partial^2 c}{\partial p_{h}\,\partial p_{h}} = \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{h}} \cdot n_2$$

$$\frac{\partial^2 c}{\partial p_{h}\,\partial p_{i}} = \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{i}} \cdot n_2$$

$$\frac{\partial^2 c}{\partial p_{h}\,\partial p_{j}} = \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial^{2} n_1}{\partial p_{h}\,\partial p_{j}} \cdot n_2$$

$$\frac{\partial^2 c}{\partial p_{h}\,\partial p_{k}} = \frac{\partial n_1}{\partial p_{h}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(i, ·)

$$\frac{\partial^2 c}{\partial p_{i}\,\partial p_{i}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{i}} \cdot n_2 + 2 \cdot \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{i}}$$

$$\frac{\partial^2 c}{\partial p_{i}\,\partial p_{j}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{j}} + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{j}} + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{i}} + \frac{\partial^{2} n_1}{\partial p_{i}\,\partial p_{j}} \cdot n_2$$

$$\frac{\partial^2 c}{\partial p_{i}\,\partial p_{k}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{i}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{i}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(j, ·)

$$\frac{\partial^2 c}{\partial p_{j}\,\partial p_{j}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{j}} + \frac{\partial^{2} n_1}{\partial p_{j}\,\partial p_{j}} \cdot n_2 + 2 \cdot \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{j}}$$

$$\frac{\partial^2 c}{\partial p_{j}\,\partial p_{k}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{j}\,\partial p_{k}} + \frac{\partial n_1}{\partial p_{j}} \cdot \frac{\partial n_2}{\partial p_{k}}$$

(k, ·)

$$\frac{\partial^2 c}{\partial p_{k}\,\partial p_{k}} = n_1 \cdot \frac{\partial^{2} n_2}{\partial p_{k}\,\partial p_{k}}$$


2 차 미분은 위 삼각(행 ≤ 열)만 적었다. 아래 삼각은 그 전치다.


검사(명세의 성분 미분과): {'samples': 3.0, 'max_rel_C': 3.4317465307146135e-16, 'max_rel_grad': 2.3476209013512296e-16, 'max_rel_hess': 8.888009462324232e-16, 'max_abs_s2_c2_minus_1': 4.440892098500626e-16}

## 5. Baraff–Witkin 전단

### BaraffShear — 텐서 꼴

Baraff-Witkin 1998 shear condition C = sqrt(a) w_u . w_v on triangle (0, 1, 2)

**조건 함수** (꼭짓점 3 개, 성분 1 개)

$$C_{1} = \sqrt{a} \mathbf{w}_u\cdot\mathbf{w}_v$$

**중간 벡터** — 좌표의 1 차식이다. $m$ 은 꼭짓점 번호(0 부터).

$$\mathbf{w}_u = \sum_m c^{u}_{m}\,\mathbf{x}_m,\qquad c^{u}_{0} = \frac{v_{1} - v_{2}}{\Delta},\quad c^{u}_{1} = \frac{- v_{0} + v_{2}}{\Delta},\quad c^{u}_{2} = \frac{v_{0} - v_{1}}{\Delta}$$

$$\mathbf{w}_v = \sum_m c^{v}_{m}\,\mathbf{x}_m,\qquad c^{v}_{0} = \frac{- u_{1} + u_{2}}{\Delta},\quad c^{v}_{1} = \frac{u_{0} - u_{2}}{\Delta},\quad c^{v}_{2} = \frac{- u_{0} + u_{1}}{\Delta}$$

여기서 $\Delta = \left(- u_{0} + u_{1}\right) \left(- v_{0} + v_{2}\right) - \left(- u_{0} + u_{2}\right) \left(- v_{0} + v_{1}\right),\quad a = \tfrac12\,|\Delta|$ 이고, 계수는 쉬는 패턴만으로 정해진다($\sum_m c_m = 0$ — 병진 불변).

**1 차 미분** — 꼭짓점 $m$ 에 대한 3-벡터

$$\frac{\partial C_{1}}{\partial \mathbf{x}_m} = \sqrt{a}\,\left(c^{u}_{m}\,\mathbf{w}_v + c^{v}_{m}\,\mathbf{w}_u\right)$$

**2 차 미분** — 꼭짓점 쌍 $(m, n)$ 의 3 × 3 블록

$$\frac{\partial^2 C_{1}}{\partial \mathbf{x}_m\,\partial \mathbf{x}_n} = \sqrt{a}\,\left(c^{u}_{m}\,c^{v}_{n} + c^{v}_{m}\,c^{u}_{n}\right)\mathbf{I}$$

**지배방정식에 들어가는 항** — 꼭짓점 블록으로. $\nabla_m C \equiv \partial C/\partial\mathbf{x}_m$, $\nabla^2_{mn} C \equiv \partial^2 C/\partial\mathbf{x}_m\partial\mathbf{x}_n$, 합은 조건 성분 $c$ 에 대해.

$$\dot C_c = \sum_n \nabla_n C_c\cdot\mathbf{v}_n$$

$$\mathbf{f}_m = -k\sum_c C_c\,\nabla_m C_c,\qquad \mathbf{d}_m = -k_d\sum_c \dot C_c\,\nabla_m C_c$$

$$\mathbf{K}_{mn} = \frac{\partial \mathbf{f}_m}{\partial \mathbf{x}_n} = -k\sum_c\Big[\nabla_m C_c\,\nabla_n C_c^{\mathsf T} + \beta_c\,\nabla^2_{mn} C_c\Big],\qquad \mathbf{V}_{mn} = \frac{\partial \mathbf{d}_m}{\partial \mathbf{v}_n} = -k_d\sum_c \nabla_m C_c\,\nabla_n C_c^{\mathsf T}$$

$$\mathbf{A} = \mathbf{M} - h\,\mathbf{V} - h^2\,\mathbf{K},\qquad \mathbf{b} = h\,(\mathbf{f} + \mathbf{d}) + h^2\,\mathbf{K}\,\mathbf{v}$$

| 강성 판 | $\beta_c$ | 생성 함수 |
|---|---|---|
| 정확 | $C_c + (k_d/k)\,\dot C_c$ | `EvalImplicit` · `EvalEngine` |
| Gauss–Newton | $0$ | `EvalImplicitGN` · `EvalEngineGN` |
| 인장 C > 0 | $\max(C_c, 0)$ | `EvalImplicitTension` · `EvalEngineTension` |
| 엔진 처방 | $\max(C_c, 0) + (k_d/k)\max(\dot C_c, 0)$ | `EvalEngineTensionRate` |


검사: {'samples': 50.0, 'max_rel_C': 2.3995821735651806e-15, 'max_rel_grad': 6.502648230039918e-16, 'max_rel_hess': 2.876152390357333e-16}

## 6. Bridson 굽힘

### BridsonBend — 텐서 꼴

Bridson 2003 bending as a condition C = 2 sqrt(2 w) sin(theta/4), theta = atan2((nA x nB).e, (nA.nB)|e|)

$$C = F(\theta) = 2 \sqrt{2} \sqrt{\bar{w}} \sin{\left(\frac{\theta}{4} \right)},\qquad F^\prime = \frac{\sqrt{2} \sqrt{\bar{w}} \cos{\left(\frac{\theta}{4} \right)}}{2},\qquad F^{\prime\prime} = - \frac{\sqrt{2} \sqrt{\bar{w}} \sin{\left(\frac{\theta}{4} \right)}}{8}$$

$\theta$ 는 두 삼각형의 법선이 이루는 각이다. 단위 원자로 쓰면($n_1 = \hat{\mathbf{n}}_A$, $n_2 = \hat{\mathbf{n}}_B$; $(\mathbf{a}\times\mathbf{b})\cdot\hat{\mathbf{e}} = -\mathbf{a}\cdot(\hat{\mathbf{e}}\times\mathbf{b})$)

$$s \equiv \sin\theta = - n_1 \cdot \tilde{\hat{e}} \cdot n_2,\qquad c \equiv \cos\theta = n_1 \cdot n_2,\qquad \theta = \operatorname{atan2}(s, c)$$

명세가 정규화하지 않은 법선 · 변으로 atan2 를 쓰더라도 atan2 는 양의 배율에 불변이라 같은 $\theta$ 다(`check()` 가 확인한다).

**연쇄 법칙** — $\nabla_m \equiv \partial/\partial p_m$, $\nabla^2_{mn} \equiv \partial^2/\partial p_m\partial p_n$. 단위 벡터에서는 $s^2 + c^2 = 1$ 이 항등식이다.

$$\nabla_m \theta = c\,\nabla_m s - s\,\nabla_m c$$

$$\nabla^2_{mn} \theta = c\,\nabla^2_{mn} s - s\,\nabla^2_{mn} c - 2sc\left(\nabla_m s\,\nabla_n s^{\mathsf T} - \nabla_m c\,\nabla_n c^{\mathsf T}\right) + (s^2 - c^2)\left(\nabla_m s\,\nabla_n c^{\mathsf T} + \nabla_m c\,\nabla_n s^{\mathsf T}\right)$$

$$\nabla_m C = F^\prime\,\nabla_m\theta,\qquad \nabla^2_{mn} C = F^\prime\,\nabla^2_{mn}\theta + F^{\prime\prime}\,\nabla_m\theta\,\nabla_n\theta^{\mathsf T}$$

**원자의 정의** — 힌지 $(h, i, j, k)$

$$\mathbf{N}_1 = (\mathbf{p}_i - \mathbf{p}_h)\times(\mathbf{p}_j - \mathbf{p}_h),\qquad \mathbf{N}_2 = (\mathbf{p}_j - \mathbf{p}_k)\times(\mathbf{p}_i - \mathbf{p}_k),\qquad \tilde{e} = [\mathbf{p}_j - \mathbf{p}_i]_\times$$

$$n_1 = \mathbf{N}_1/|\mathbf{N}_1|,\qquad n_2 = \mathbf{N}_2/|\mathbf{N}_2|$$

$[\mathbf{a}]_\times$ 는 $\mathbf{a}\times$ 의 3 × 3 반대칭 행렬이다.

**원자의 1 차 미분** (3 × 3 블록)

$$\frac{\partial \mathbf{N}_1}{\partial p_h} = [\mathbf{p}_j - \mathbf{p}_i]_\times = \tilde{e},\qquad \frac{\partial \mathbf{N}_1}{\partial p_i} = [\mathbf{p}_h - \mathbf{p}_j]_\times,\qquad \frac{\partial \mathbf{N}_1}{\partial p_j} = [\mathbf{p}_i - \mathbf{p}_h]_\times$$

$$\frac{\partial \mathbf{N}_2}{\partial p_k} = [\mathbf{p}_i - \mathbf{p}_j]_\times = -\tilde{e},\qquad \frac{\partial \mathbf{N}_2}{\partial p_j} = [\mathbf{p}_k - \mathbf{p}_i]_\times,\qquad \frac{\partial \mathbf{N}_2}{\partial p_i} = [\mathbf{p}_j - \mathbf{p}_k]_\times$$

$$\frac{\partial n}{\partial p_m} = \frac{1}{|\mathbf{N}|}\left(\mathbf{I} - n\,n^{\mathsf T}\right)\frac{\partial \mathbf{N}}{\partial p_m}$$

$$\frac{\partial \tilde{e}}{\partial p_j}[\boldsymbol\delta] = [\boldsymbol\delta]_\times,\qquad \frac{\partial \tilde{e}}{\partial p_i}[\boldsymbol\delta] = -[\boldsymbol\delta]_\times$$

**원자의 2 차 미분** — 방향 $\mathbf{a}$(첫 꼭짓점) · $\mathbf{b}$(둘째 꼭짓점)에 대한 값. $\mathbf{N} = (P - H)\times(Q - H)$ 로 쓰면($\mathbf{N}_1$: $H{=}h, P{=}i, Q{=}j$; $\mathbf{N}_2$: $H{=}k, P{=}j, Q{=}i$)

$$D^2\mathbf{N}[\mathbf{a}, \mathbf{b}] = \sigma\;\mathbf{a}\times\mathbf{b},\qquad \sigma_{PQ} = \sigma_{HP} = \sigma_{QH} = +1,\quad \sigma_{QP} = \sigma_{PH} = \sigma_{HQ} = -1,\quad \sigma_{PP} = \sigma_{QQ} = \sigma_{HH} = 0$$

$$D^2 n[\mathbf{a}, \mathbf{b}] = \frac{\mathbf{P}\,D^2\mathbf{N}[\mathbf{a}, \mathbf{b}]}{|\mathbf{N}|} - \frac{(n\cdot \mathbf{N}_b)\,\mathbf{P}\mathbf{N}_a + (n\cdot \mathbf{N}_a)\,\mathbf{P}\mathbf{N}_b + (\mathbf{N}_a\cdot\mathbf{P}\mathbf{N}_b)\,n}{|\mathbf{N}|^2},\qquad \mathbf{P} = \mathbf{I} - n\,n^{\mathsf T},\quad \mathbf{N}_a = D\mathbf{N}[\mathbf{a}]$$

$\tilde{e}$ 는 좌표의 1 차식이라 2 차 미분이 0 이다.

**단위 변** $\hat{\mathbf{e}} = \mathbf{e}/|\mathbf{e}|$, $\tilde{\hat{e}} = [\hat{\mathbf{e}}]_\times$ — 단위 법선과 같은 꼴이다($\mathbf{e}_a = D\mathbf{e}[\mathbf{a}]$: 꼭짓점 $j$ 이면 $\mathbf{a}$, $i$ 이면 $-\mathbf{a}$).

$$D\hat{\mathbf{e}}[\mathbf{a}] = \frac{\mathbf{P}_e\,\mathbf{e}_a}{|\mathbf{e}|},\qquad D^2\hat{\mathbf{e}}[\mathbf{a}, \mathbf{b}] = -\frac{(\hat{\mathbf{e}}\cdot\mathbf{e}_b)\,\mathbf{P}_e\mathbf{e}_a + (\hat{\mathbf{e}}\cdot\mathbf{e}_a)\,\mathbf{P}_e\mathbf{e}_b + (\mathbf{e}_a\cdot\mathbf{P}_e\mathbf{e}_b)\,\hat{\mathbf{e}}}{|\mathbf{e}|^2},\qquad \mathbf{P}_e = \mathbf{I} - \hat{\mathbf{e}}\hat{\mathbf{e}}^{\mathsf T}$$



검사(명세의 성분 미분과): {'samples': 3.0, 'max_rel_C': 1.5624440492200542e-16, 'max_rel_grad': 3.2066383784035913e-16, 'max_rel_hess': 1.5721708232970782e-15, 'max_abs_s2_c2_minus_1': 4.440892098500626e-16}
