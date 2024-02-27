!----------------------------------------------------------------------------------------------------------
! Exact (Godnov) Riemann solution  Ver 1.0 , 2010-05-20
! Copyright by Li Xinliang , LHD, Institute of Mechanics, CAS,  lixl@imech.ac.cn
! Riemann Problem :  1d Euler equation -1<x<1;  
!                    initial condition:  (r,u,p)= (rL,uL,pL) if x<0;  (r,u,p)=(rR,uR,pR)  else where
! Ref. :  http://www.cfluid.com/bbs/viewthread.php?tid=75637&extra=page%3D1     (see the PPT files: Lecture 2 and Lecture 9)
! Ref. :  傅德薰，马延文：《计算流体力学》 p29-34  
! 版权所有： 李新亮， 中国科学院力学研究所LHD实验室， lixl@imech.ac.cn
! 欢迎使用和传播本程序。 若在科研中使用本程序，请在论文（报告）中进行标注 （可引用作者的文献或致谢）。
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!  r_L, u_L, p_L : initial density, velocity and pressure in the  Lift of the initial discontinuity;  
!  r_R, u_R, p_R : initial density, velocity and pressure in the Right of the initial discontinuity;
!  r,u,p : dimension(nx), one dimensioanal distribution of density, velocity and pressure at time t;  
!  c : speed of sound, Z : speed of the shockwave ; Z_head, Z_tail : Speed of the rarefaction's head and rarafaction's tail  
!  x: diemsion (nx), x-coordinate;  t: time  ; gamma (=cp/cv=1.4)
!  subscript "_L" and "_R": variables in the Lift and the right of the perturbation waves ; 
!  subscript "_star" : variables in the central region  (between the lift wave and the right wave);
!  subscript "_star_L", "_star_R" : variables in the lift and right side of the centact discontinuity;   
!-------------------------------------------------------------------------------------------------------------
   program main
   implicit none 
   real*8 r_L,u_L,p_L,r_R,u_R,p_R,t,gamma
   integer nx,k
   real*8,allocatable:: r(:),u(:),p(:),x(:)
   gamma=1.4d0
   print*,  "---------Exact Riemann Solver, ver 1.0---------------------" 
   print*, "Copyright by Li Xinliang , Institute of Mechanics, CAS"
   print*, "please input nx (grid number) for plot,  e.g. 201"
   read(*,*) nx
   allocate(r(nx),u(nx),p(nx),x(nx))
   print*, "please input r_L, u_L, p_L, r_R, u_R, p_R"
   print*, "density,velocity and pressure in the lift and right side"
   print*, "example:  1, 0, 1, 0.125, 0, 0.1"
   read(*,*) r_L,u_L,p_L,r_R,u_R,p_R
   print*, "please input time t, example: t=0.14"
   read(*,*) t
   print*, "--------------------------------------"
   do k=1,nx
    x(k)=1.d0*(k-1.d0)/(nx-1.d0)-0.5d0
   enddo

   call Exact_Riemann_solver(r_L,u_L,p_L,r_R,u_R,p_R,r,u,p,x,nx,t,gamma)
   
   open(99,file="Riemann.dat")
   write(99,*) "variables=x,r,u,p"
   do k=1,nx
   write(99,"(4E20.12)") x(k),r(k),u(k),p(k)
   enddo
   close(99)
   print*, "OK, the flow data are writen to 'Riemann.dat' as a tecplot file"
   deallocate(r,u,p,x)
   end

!-------------------------------------------------------------------------------------------------------------

    subroutine Exact_Riemann_solver(r_L,u_L,p_L,r_R,u_R,p_R,r,u,p,x,nx,t,gamma)
    implicit none
    integer k, nx
    real*8:: r(nx),u(nx),p(nx),x(nx),t
    real*8:: p_L,r_L,u_L,p_R,r_R,u_R   ! p: pressure; r: density; u: velocity
    real*8:: p_star, u_star,r_star_L,r_star_R, c_star_L, c_star_R   ! variables in the central region between the lift and the right waves
    real*8:: gamma,f_Riemann, Fp_Riemann, du, Fp_L,Fp_R,Fp_0, c_L,c_R,A_L,A_R, Z_L,Z_R, Z_L_head,Z_L_tail,Z_R_head,Z_R_tail, ca

!   
    print*, "Exact Riemann Solver ......" 
    print*, "Copyright by Li Xinliang , LHD, Institute of Mechanics, CAS"

    c_L=sqrt(gamma*p_L/r_L) ; c_R=sqrt(gamma*p_R/r_R)    ! sound speed
    call get_pstar_Newton(p_star,p_L,r_L,u_L,p_R,r_R,u_R,gamma)    ! get the pressure in the central region
   
    u_star=(u_R+u_L+f_Riemann(p_star,p_R,r_R,gamma)-f_Riemann(p_star,p_L,r_L,gamma))/2.d0
    A_L=r_L*c_L*sqrt((gamma+1.d0)/(2.d0*gamma)*p_star/p_L+(gamma-1.d0)/(2.d0*gamma))    ! A1
    A_R=r_R*c_R*sqrt((gamma+1.d0)/(2.d0*gamma)*p_star/p_R+(gamma-1.d0)/(2.d0*gamma))    ! A2

!   To judge which case will be occur :  case 1 to case 5------------------------------------------------
    du=u_L-u_R 
    FP_L=FP_Riemann(p_L,p_L,r_L,p_R,r_R,gamma)
    FP_R=FP_Riemann(p_R,p_L,r_L,p_R,r_R,gamma)
    Fp_0=Fp_Riemann(0.d0,p_L,r_L,p_R,r_R,gamma)
 
  if( du .ge. max(Fp_L,FP_R) ) then
 !    case 1 : Lift & Right waves are all shockwaves; --------------------------------
      Z_L=u_L-A_L/r_L;  r_star_L=r_L*A_L/(A_L-r_L*(u_L-u_star))
      Z_R=u_R+A_R/r_R;  r_star_R=r_R*A_R/(A_R+r_R*(u_R-u_star))
      do k=1,nx
       if( x(k) .lt. Z_L*t) then
        r(k)=r_L ; u(k)=u_L ; p(k)=p_L
       else if ( x(k) .lt. u_star*t ) then
        r(k)=r_star_L ; u(k)=u_star ; p(k)=p_star
       else if (x(k) .lt. Z_R*t) then
        r(k)=r_star_R ; u(k)=u_star; p(k)=p_star
       else
        r(k)=r_R ; u(k)=u_R ; p(k)=p_R
       endif
      enddo

  else if( du .lt. Fp_0 ) then
 ! case 5: Lift & Right waves are all rarefactions, and the center is vacuum ;
     Z_L_head=u_L-c_L; Z_L_tail=u_L+2.d0/(gamma-1.d0)*c_L 
     Z_R_head=u_R+c_R; Z_R_tail=u_R-2.d0/(gamma-1.d0)*c_R 
      do k=1,nx
       if(x(k) .lt. Z_L_head*t) then
         r(k)=r_L ; u(k)=u_L ; p(k)=p_L
       else if (x(k) .lt. Z_L_tail*t) then
         ca=(gamma-1.d0)/(gamma+1.d0)*((u_L-x(k)/t)+2.d0/(gamma-1.d0)*c_L )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 32. !!
         u(k)=x(k)/t+ca ;    p(k)=p_L*(ca/c_L)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else if (x(k) .lt. Z_R_tail*t) then
         u(k)=0.d0; p(k)=0.d0; r(k)=0.d0 ;   ! Vacuum region
       else if (x(k) .lt. Z_R_head*t ) then
          ca=(gamma-1.d0)/(gamma+1.d0)*((x(k)/t-u_R)+2.d0/(gamma-1.d0)*c_R )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 33. !!
          u(k)=x(k)/t-ca ;    p(k)=p_R*(ca/c_R)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else
        r(k)=r_R ; u(k)=u_R ; p(k)=p_R
       endif
      enddo
  
  else if(du .lt. min(Fp_L,Fp_R)) then
   ! case 4 : Lift & Right waves are all rarefactions; 
     c_star_L=c_L+(gamma-1.d0)/2.d0*(u_L-u_star); 
     Z_L_head=u_L-c_L; Z_L_tail=u_star-c_star_L ; 
     r_star_L=gamma*p_star/(c_star_L*c_star_L) 
     c_star_R=c_R-(gamma-1.d0)/2.d0*(u_R-u_star); 
     Z_R_head=u_R+c_R; Z_R_tail=u_star+c_star_R ; 
     r_star_R=gamma*p_star/(c_star_R*c_star_R) 

      do k=1,nx
       if(x(k) .lt. Z_L_head*t) then
         r(k)=r_L ; u(k)=u_L ; p(k)=p_L
       else if (x(k) .lt. Z_L_tail*t) then
         ca=(gamma-1.d0)/(gamma+1.d0)*((u_L-x(k)/t)+2.d0/(gamma-1.d0)*c_L )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 32. !!
         u(k)=x(k)/t+ca ;    p(k)=p_L*(ca/c_L)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else if ( x(k) .lt. u_star*t ) then
        r(k)=r_star_L ; u(k)=u_star ; p(k)=p_star
       else if (x(k) .lt. Z_R_tail*t) then
        r(k)=r_star_R ; u(k)=u_star; p(k)=p_star
       else if (x(k) .lt. Z_R_head*t ) then
          ca=(gamma-1.d0)/(gamma+1.d0)*((x(k)/t-u_R)+2.d0/(gamma-1.d0)*c_R )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 33. !!
          u(k)=x(k)/t-ca ;    p(k)=p_R*(ca/c_R)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else
        r(k)=r_R ; u(k)=u_R ; p(k)=p_R
       endif
      enddo

  else  if(p_L .gt. p_R) then
 ! case 2: Lift wave is rarefaction, Right wave is shockwave;
     c_star_L=c_L+(gamma-1.d0)/2.d0*(u_L-u_star); Z_L_head=u_L-c_L; 
     Z_L_tail=u_star-c_star_L ; 
     r_star_L=gamma*p_star/(c_star_L*c_star_L)  ! lift rarefaction
     Z_R=u_R+A_R/r_R;  
     r_star_R=r_R*A_R/(A_R+r_R*(u_R-u_star))   !right shockwave

      do k=1,nx
       if(x(k) .lt. Z_L_head*t) then
         r(k)=r_L ; u(k)=u_L ; p(k)=p_L
       else if (x(k) .lt. Z_L_tail*t) then
         ca=(gamma-1.d0)/(gamma+1.d0)*((u_L-x(k)/t)+2.d0/(gamma-1.d0)*c_L )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 32. !!
         u(k)=x(k)/t+ca ;    p(k)=p_L*(ca/c_L)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else if ( x(k) .lt. u_star*t ) then
        r(k)=r_star_L ; u(k)=u_star ; p(k)=p_star
       else if (x(k) .lt. Z_R*t) then
        r(k)=r_star_R ; u(k)=u_star; p(k)=p_star
       else
        r(k)=r_R ; u(k)=u_R ; p(k)=p_R
       endif
      enddo
 else
! case 3:  Lift wave is shockwave, Right wave is rarefaction;
      Z_L=u_L-A_L/r_L;  r_star_L=r_L*A_L/(A_L-r_L*(u_L-u_star))
      c_star_R=c_R-(gamma-1.d0)/2.d0*(u_R-u_star); 
      Z_R_head=u_R+c_R; Z_R_tail=u_star+c_star_R ; 
      r_star_R=gamma*p_star/(c_star_R*c_star_R) 
      do k=1,nx
       if(x(k) .lt. Z_L*t) then
         r(k)=r_L ; u(k)=u_L ; p(k)=p_L
       else if ( x(k) .lt. u_star*t ) then
        r(k)=r_star_L ; u(k)=u_star ; p(k)=p_star
       else if (x(k) .lt. Z_R_tail*t) then
        r(k)=r_star_R ; u(k)=u_star; p(k)=p_star
       else if (x(k) .lt. Z_R_head*t ) then
          ca=(gamma-1.d0)/(gamma+1.d0)*((x(k)/t-u_R)+2.d0/(gamma-1.d0)*c_R )   ! sound speed in the rarefaction.  There is an error in Fu's Book page 33. !!
          u(k)=x(k)/t-ca ;    p(k)=p_R*(ca/c_R)**(2.d0*gamma/(gamma-1.d0)) ; r(k)=gamma*p(k)/(ca*ca)
       else
        r(k)=r_R ; u(k)=u_R ; p(k)=p_R
       endif
      enddo
  endif 
 end


!----------------------------------------------------------------------------------------------------
! get pressure in the central region by using Newton iteration
! solve the problem:  F(p_star)=u1-u2    
    subroutine get_pstar_Newton(p_star,p_L,r_L,u_L,p_R,r_R,u_R,gamma)
    implicit none
    real*8 p_star,p_L,r_L,u_L,p_R,r_R,u_R,gamma,delt_p,du,Fp_Riemann,dF_dp,p_star_new,Err,Err_level
    Err_level=1.d-20 ;   delt_p=1.d-4  ! get derivative dF/dp  
    du=u_L-u_R
    p_star=min(p_L,p_R)   ! initial data
10   dF_dp=(FP_Riemann(p_star+delt_p,p_L,r_L,p_R,r_R,gamma)-FP_Riemann(p_star,p_L,r_L,p_R,r_R,gamma))/delt_p  ! F'(p_star)
     p_star_new=p_star-(FP_Riemann(p_star,p_L,r_L,p_R,r_R,gamma)-du)/dF_dp 
     if (p_star_new .lt. 0) p_star_new=0
     Err=abs(p_star_new-p_star)
     p_star=p_star_new
     if(Err .gt. 1.d-20 ) goto 10 
     print*, "p_star=", p_star_new, "abs(p_star_new-p_star)=",Err
    end

!---------------------------------------------------------
    function FP_Riemann(p,p_L,r_L,p_R,r_R,gamma)
    implicit none
    real*8:: Fp_Riemann,p,p_L,r_L,p_R,r_R,gamma,f_Riemann
       Fp_Riemann=f_Riemann(p,p_L,r_L,gamma)+f_Riemann(p,p_R,r_R,gamma)
    end
   
!  function f(p) ---------------------------     
    function f_Riemann(p,p1,r1,gamma)
    implicit none
    real*8:: f_Riemann,p,p1,r1,c1,gamma
    c1=sqrt(gamma*p1/r1)
    if( p .lt. 0) then
      print*, "Error, pressure < 0 !!! "
    else if( p .lt. 1.d-40) then   !! if (p==0)
      f_Riemann=-2.d0*c1/(gamma-1.d0)
    else
        if(p .gt. p1) then
           f_Riemann=(p-p1)/(r1*c1*sqrt((gamma+1.d0)/(2.d0*gamma)*p/p1+(gamma-1.d0)/(2.d0*gamma)) )
        else
           f_Riemann=2.d0*c1/(gamma-1.d0)*((p/p1)**((gamma-1.d0)/(2.d0*gamma))-1.d0)
        endif
    endif 
    end
!--------------------------------------------------------