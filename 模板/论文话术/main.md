# 用词用语库

utilize classical

To avoid the denominator being zero, the small constant ε can be taken as 10−6.

In this paper, we denote by MDFLU the scheme (17)–(19). 

Combining (20) and (29) gives the high-order WENO scheme, which is termed as MDFLU-WENO5.
Together with (21) and (29), we obtain another high-order accurate WENO scheme, which is called DFLU-WENO5.
So, we obtain the fifth high-order accurate scheme with the well-balanced method and call the scheme DFLU-WENO5B.

We solve this scalar advection equation up to time t = 20. 

The computations are conducted on four different mesh resolutions with the number of cells taken as N = 502, 1002, 2002 and 4002.

The numerical results are displayed in Figure 4. 

The simulation was performed up to time  $t=2$.

We can clearly find that the AWENO-E scheme achieves the designed order of accuracy.

Oscillation evolution for the 1D sharp interface advection problem

Density and bound error for the 1D sharp interface advection problem 

Velocity and pressure for the 1D sharp interface advection problem

Density for the 1D sharp interface advection problem


with no spurious velocity and pressure oscillations across a material interface

Certainly, the results of WENO-NE6 are less diffusive than those of other sixth-order WENO schemes in the discontinuities, showing excellent agreement with the exact solution.

# 方法
This section presents some numerical experiments to demonstrate the performance of the proposed WENO scheme. The sixth-order numerical flux $\hat{f}$ interpolating cell-average values on the stencil $S_6$ is obtained from the space

$\Gamma_6 = \text{span}\{1, x, e^{\lambda x}, e^{-\lambda x}, \cos(\lambda x), \sin(\lambda x)\}$

with $\lambda > 0$. The numerical results of WENO-NE6 are compared with those from the adaptive central-upwind WENO scheme (WENO-CU6) \cite{13} for several test problems. In addition, it is worth pointing out that the set of six algebraic polynomials, $\{x^n : n = 0, \ldots, 5\}$, is a special case of exponential polynomial bases. The WENO scheme based on these algebraic polynomials (with the same smoothness indicators in (5.3)) is termed as WENO-NP6. Thus, in order to address the potential of using exponential polynomials, the performance of WENO-NE6 is also compared with that of WENO-NP6.

---

# 算例

## 精度
- In Fig. 7, we plot the computational accuracy which indicates that the present $\mathrm{DG}\left(P^{K}\right)$ method is able to achieve the optimal $(K+1)$-th convergence order for the two-phase vortex advection problem. 

- The initial and final mixture density fields obtained by the  $\operatorname{DG}\left(P^{2}\right)$  method on the finest mesh is shown in Fig. 8. It can be seen that the final result exhibits a good agreement with the initial solution.

- The convergence statistics from the fifth-order upwind linear scheme and the present hybrid TENO5-THINC are given in Table 2. It is clear that the desired fifth-order convergence is achieved without order degeneration with the present nonlinear shock-capturing scheme. Moreover, the absolute errors from TENO5-THINC are identical to those from the upwind linear scheme, indicating that neither the nonlinear adaptation of TENO nor the THINC reconstruction is activated in smooth
regions.
- A sequence of globally refined uniform grids is employed to investigate the convergence in terms of the L∞ norm. As shown in Fig. 4, the resolved profile from the present scheme agrees well with the analytical solution, and the desired fifth-order convergence is achieved without order degeneration.

- Example 4.2. In this example, we test the accuracy of the AWENO-E scheme for the two-dimensional system of Euler equations in (47) with the initial condition $\mathbf{q}(x, y, 0)=(\rho, u, v, p)=(1+0.5 \sin (4 \pi(x+y)), 1,-0.5,1), \quad(x, y) \in[-1,1] \times[-1,1],$ where periodic boundary conditions are used. The simulation was performed up to time  $t=2$. 

The approximation errors computed by the AWENO-E scheme (with LF and HLLC solvers) are smaller than those by the AWENO-JS and AWENO-Z schemes. In particular, the  $L_{1}$  errors obtained by the AWENO-E with HLLC solver are slightly smaller than those computed by AWENO-E with LF solver, especially on low-resolution grids. Further, the computational efficiencies of the tested AWENO methods are compared in terms of the CPU time versus  $L_{1}$  or $L_{\infty}$  errors using different grids. The graphs of CPU times against errors are displayed in Fig. 2. Each marker in these graphs depicts the pairs of CPU time and  L_{\infty} -approximation errors at  10 \times 2^{k}  grid points with  k=0, \ldots, 6 .

- Periodic boundary conditions are imposed at  x=0  and  x= 
The solution is advanced to  t=1  corresponding to one period in time, and a sequence of globally refined uniform grids is employed to investigate the  L_{\infty} -error convergence. The timestep
is decreased so that the time-integration error can be neglected. As shown in Fig. 10, the  L_{\infty} -error convergence histories of all proposed TENO schemes coincide with their corresponding background linear schemes and therefore the desirable order of accuracy is achieved. 

## 间断
Next, we consider a sharp initial volume fraction given as follows

xxxxx

where εz = 10−4 and r0 = L/10 is the half length of the interface. The output time is T = 1s and all the other variables remain the same as in the previous smooth interface advection problem. This case is used to verify the ability of the present method in maintaining the oscillation-free property of velocity and pressure at an isolated phase interface. Since the solution contains sharp gradient regions, the nonlinear WENO limiter is applied in this test case.

In Fig. 2 we compare the numerical solutions of ρ obtained by the DG-WENO method with different polynomial approximation orders and mesh resolutions. We can observe that the initial sharp interface becomes smeared as time evolves. While with the increase of the polynomial approximation order and the number of mesh cells, the numerical solutions are getting closer to the initial solution. In order to verify the oscillation-free property of the proposed method, we present the numerical solutions and the oscillation evolutions of velocity and pressure in Fig. 3 and Fig. 4, respectively. It is obvious to see that uniform velocity and pressure fields with nearly zero oscillation, e.g. around O(10−10) and O(10−13) for velocity and pressure, respectively, can be ensured by the present DG-WENO method.

In the end, we test the feasibility of the positivity-preserving limiter. In order to initialize the lower and upper bounds of the volume fraction zk more close to 0 and 1, respectively, we decrease εz to 10−8 in (50). Under this case, since the high-order DG method violates the monotonicity, negative partial density appears in the vicinity of the phase interface which leads to a break down during the simulation. Thus, the DG-WENO-PP method is applied in this test case for improving robustness. Fig. 5(a) shows the mixture density obtained by the DG-WENO-PP method with 320 uniform cells. We can observe that the positivity-preserving limiter improves the computational robustness significantly and the present DG-WENO-PP method can capture the high and low mixture density regions accurately. In order to estimate the potential undershoots near the phase interface, we introduce the bound error of
