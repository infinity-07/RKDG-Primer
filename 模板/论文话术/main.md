# 用词用语库
We solve this scalar advection equation up to time t = 20. 

The numerical results are displayed in Figure 4. 

with no spurious velocity and pressure oscillations across a material interface

Certainly, the results of WENO-NE6 are less diffusive than those of other sixth-order WENO schemes in the discontinuities, showing excellent agreement with the exact solution.

# 方法
This section presents some numerical experiments to demonstrate the performance of the proposed WENO scheme. The sixth-order numerical flux $\hat{f}$ interpolating cell-average values on the stencil $S_6$ is obtained from the space

$\Gamma_6 = \text{span}\{1, x, e^{\lambda x}, e^{-\lambda x}, \cos(\lambda x), \sin(\lambda x)\}$

with $\lambda > 0$. The numerical results of WENO-NE6 are compared with those from the adaptive central-upwind WENO scheme (WENO-CU6) \cite{13} for several test problems. In addition, it is worth pointing out that the set of six algebraic polynomials, $\{x^n : n = 0, \ldots, 5\}$, is a special case of exponential polynomial bases. The WENO scheme based on these algebraic polynomials (with the same smoothness indicators in (5.3)) is termed as WENO-NP6. Thus, in order to address the potential of using exponential polynomials, the performance of WENO-NE6 is also compared with that of WENO-NP6.

---

# 算例

Next, we consider a sharp initial volume fraction given as follows

xxxxx

where εz = 10−4 and r0 = L/10 is the half length of the interface. The output time is T = 1s and all the other variables remain the same as in the previous smooth interface advection problem. This case is used to verify the ability of the present method in maintaining the oscillation-free property of velocity and pressure at an isolated phase interface. Since the solution contains sharp gradient regions, the nonlinear WENO limiter is applied in this test case.

In Fig. 2 we compare the numerical solutions of ρ obtained by the DG-WENO method with different polynomial approximation orders and mesh resolutions. We can observe that the initial sharp interface becomes smeared as time evolves. While with the increase of the polynomial approximation order and the number of mesh cells, the numerical solutions are getting closer to the initial solution. In order to verify the oscillation-free property of the proposed method, we present the numerical solutions and the oscillation evolutions of velocity and pressure in Fig. 3 and Fig. 4, respectively. It is obvious to see that uniform velocity and pressure fields with nearly zero oscillation, e.g. around O(10−10) and O(10−13) for velocity and pressure, respectively, can be ensured by the present DG-WENO method.

In the end, we test the feasibility of the positivity-preserving limiter. In order to initialize the lower and upper bounds of the volume fraction zk more close to 0 and 1, respectively, we decrease εz to 10−8 in (50). Under this case, since the high-order DG method violates the monotonicity, negative partial density appears in the vicinity of the phase interface which leads to a break down during the simulation. Thus, the DG-WENO-PP method is applied in this test case for improving robustness. Fig. 5(a) shows the mixture density obtained by the DG-WENO-PP method with 320 uniform cells. We can observe that the positivity-preserving limiter improves the computational robustness significantly and the present DG-WENO-PP method can capture the high and low mixture density regions accurately. In order to estimate the potential undershoots near the phase interface, we introduce the bound error of
