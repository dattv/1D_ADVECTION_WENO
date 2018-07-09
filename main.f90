! $HeadURL$
! $Id$

!> \mainpage
!! This code, which is based on one provided though <a
!! href="http://www.its.caltech.edu/~appelo/ae232/">Caltech's AE 232
!! course</a>, solves inviscid or viscid 1D scalar conservation law of the
!! form \f$ u_t + f(u)_x = 0 \f$ or \f$ u_t + f(u)_x = \nu u_{xx} \f$ using
!! WENO reconstruction and a Lax-Friedrichs flux.
!!
!! @see main() for more detailed information on the numerics.
!! @see \ref Makefile for details on how the code is built.


!> This program solves inviscid or viscid 1D scalar conservation laws
!! of the form \f$ u_t + f(u)_x = 0 \f$ or \f$ u_t + f(u)_x = \nu u_{xx} \f$
!! using WENO reconstruction and a Lax-Friedrichs flux on the domain
!! \f$x\in\left[0,1\right]\f$.  Time stepping is performed using a 3rd-order
!! TVD Runge-Kutta scheme.  Periodic boundary conditions are used, but that
!! is hidden in subroutines.  Different reconstruction and viscous term
!! behavior is controlled using the compile-time macro definitions of
!! <tt>WENOORDER</tt> and <tt>VISCOUSORDER</tt>.  All results are saved into
!! a single HDF5 file for easy import into analysis software.  Initial
!! conditions are currently hardcoded into this routine.
!!
!! @see flux() for the Lax-Friedrichs %flux details.
!! @see reconstruct3() or reconstruct5() for the WENO reconstruction details.
!! @see viscous2() or viscous4() for the viscous term details.
!! @see Section 2 of Liu, Osher, and Chan's 1994 JCP paper for a
!!      synopsis of the TVD Runge-Kutta scheme employed.
PROGRAM main
!  USE doublePrecision
    USE LIMITER_MODULE
    USE IR_Precision
    USE ADV2D_MODULE
  IMPLICIT NONE

! Problem parameters
  REAL(R8P), PARAMETER :: pi      = 4.d0*ATAN(1.d0)
  REAL(R8P), PARAMETER :: tend    = 0.20d0
  REAL(R8P), PARAMETER :: cfl     = 0.5d0
  REAL(R8P), PARAMETER :: nu      = 1.d0/500.d0
  REAL(R8P), PARAMETER :: rk(2,2) = RESHAPE( &              ! TVD RK3
             [0.75d0, (1.d0/3.d0), (0.25d0), (2.d0/3.d0)], SHAPE(rk))

! Working storage
  INTEGER :: i, n, nt, nsteps
  REAL(R8P) :: h, hi, t, dt
  REAL(R8P), DIMENSION(:),   ALLOCATABLE :: x, xh, u, frp, frm, fp, fm, f
  REAL(R8P), DIMENSION(:,:), ALLOCATABLE :: up
  REAL(R8P), DIMENSION(:), ALLOCATABLE :: UNM1
  REAL(R8P), DIMENSION(:), ALLOCATABLE :: UNEW
!  CHARACTER(LEN = 20) :: str
  
  REAL(R8P), DIMENSION(:),    ALLOCATABLE :: UL, UR

! Handles related to to HDF5 operations
!  INTEGER(HID_T) :: file_hid, filespace_hid, dataspace_hid, dataset_hid
!  INTEGER(HID_T) :: memspace_hid
  INTEGER        :: error, SCHEME 

! Sanity check incoming argument count
!  IF (command_argument_count() /= 2) THEN
!    CALL get_command_argument(0, str)
!    print *, "Usage: ", trim(str), " outputfile npoints"
!    CALL EXIT(1)
!  END IF

! Open the output file from arguments
!  CALL h5open_f(error)
!  CALL get_command_argument (1, str)
!  CALL h5fcreate_f(trim(str), H5F_ACC_TRUNC_F, file_hid, error)

! Determine the grid size from arguments and allocate storage
!  CALL get_command_argument (2, str)
!  READ (str, fmt = '(I10)') n

  !CALL ADV2D()
  !STOP
  n = 6000
  ALLOCATE ( x(0:n), xh(0:n), u(0:n), up(0:n,1:2), &
             frp(0:n), frm(0:n), f(0:n), fp(0:n), fm(0:n) )
  ALLOCATE(UNM1(0:N), UNEW(0:N))

  ALLOCATE ( UL(0:N), UR(0:N) )
! Setup grid x \in [0,1] containing n points
  h = 1.d0 / n
  hi = 1.d0 / h
  FORALL(i=0:n:1) x(i) = i*h
  xh = x + 0.5d0*h
!  CALL h5ltmake_dataset_double_f( &
!       file_hid, "x", 1, [INTEGER(HSIZE_T)::n], x, error)
!  CALL h5ltmake_dataset_double_f( &
!       file_hid, "xh", 1, [INTEGER(HSIZE_T)::n], xh, error)

! Determine the WENO reconstruction order from #defines
!#if WENOORDER == 3 
!#define RECONSTRUCT_FUNCTION reconstruct3
!#elif WENOORDER == 5
!#define RECONSTRUCT_FUNCTION reconstruct5
!#else
!  #error "WENOORDER not #defined or unknown"
!#endif
!  PRINT '(" WENO reconstruction order = ", I2)', WENOORDER
!  CALL h5ltmake_dataset_int_f( &
!       file_hid, "weno", 1, [INTEGER(HSIZE_T)::1], [WENOORDER], error)

! Determine viscous term presence and/or order from #defines
!#ifdef VISCOUSORDER
!  PRINT '(" Viscous term order = ", I2, " with viscosity = ", F16.12)', &
!        VISCOUSORDER, nu
!  CALL h5ltmake_dataset_int_f( &
!       file_hid, "viscous", 1, [INTEGER(HSIZE_T)::1], [VISCOUSORDER], error)
!  CALL h5ltmake_dataset_double_f( &
!       file_hid, "nu", 1, [INTEGER(HSIZE_T)::1], [nu], error)
!#if VISCOUSORDER == 2
!#define VISCOUS_FUNCTION viscous2
!#elif VISCOUSORDER == 4
!#define VISCOUS_FUNCTION viscous4
!#else
!  #error "VISCOUSORDER unknown"
!#endif
!#else
!#define VISCOUS_FUNCTION viscousnop
!  PRINT '(" No viscous term present ")'
!#endif

! Determine stable time step size
  dt = cfl * h /3. ! Convective stability
!#ifdef VISCOUSORDER
!! Constant from 1 + x + x**2/2 + x**3/6 ~ 1 along real axis per TVD RK3
!  dt = MIN(dt, 2.512_dp / (nu * pi**2 * n **2)) ! Diffusive stability
!#endif
  nsteps = INT(tend / dt)
  dt = tend / REAL(nsteps)
  PRINT '(" Number of time steps = ", I7, " with dt = ", F16.12)', &
        nsteps, dt
!  CALL h5ltmake_dataset_double_f( &
!       file_hid, "t", 1, [INTEGER(HSIZE_T)::nsteps+1], &
!       [(dt*i, i=0, nsteps, 1)], error)

! Initial condition 1: simple sine wave for debugging
! u = SIN(2._dp*pi*x)

! Initial condition 2: viscous analytic solution from Hopf-Cole
  t = 0
  u =   (                          198*pi*SIN(1 + 4*pi*x) ) &
      / ( 125*(100*EXP(4*pi*pi*t/125)+ 99*COS(1 + 4*pi*x)))
  do i = 0, n, 1
      if (x(i) <=0.5) u(i) = 1.D0
      if (x(i) >0.5) u(i) = 0.D0
      !if (x(i) <=0.1) u(i) = 0.
  end do
  
  write(*, *) "INITIAL: "
  open(unit = 12, file = "INITIAL.TEC")
    do i = 0, n, 1
        write(12, *) x(i), u(i)
    end do
  close(unit = 12)  

!! Create dataset to store initial condition and (space x time) solution
!  CALL h5screate_simple_f( &
!       2, [INTEGER(HSIZE_T)::n, nsteps + 1], dataspace_hid, error)
!  CALL h5dcreate_f( &
!       file_hid, "u", H5T_NATIVE_DOUBLE, dataspace_hid, dataset_hid, error)
!! Create dataspace describing how solution is stored in memory
!  CALL h5screate_simple_f(2, [INTEGER(HSIZE_T)::n, 1], memspace_hid, error)
!! Create dataspace describing how one solution time is stored on disk
!  CALL h5dget_space_f(dataset_hid, filespace_hid, error)

!! Write the initial condition
!  CALL h5sselect_hyperslab_f(filespace_hid, H5S_SELECT_SET_F, &
!       [INTEGER(HSIZE_T)::0, 0], [INTEGER(HSIZE_T)::n, 1], error)
!  CALL h5dwrite_f(dataset_hid, H5T_NATIVE_DOUBLE, u(0:n-1), &
!       [INTEGER(HSIZE_T)::n, 1], error, memspace_hid, filespace_hid)
    UNM1 = U
! Advance the solution in time using explicit TVD RK3
  SCHEME = 0
  timeloop: DO nt = 1, nsteps, 1
    t = REAL(nt-1) * dt
    
    
    IF (SCHEME == 0) THEN !@ SCHEME => 0 IS EXPLICIT SECOND ORDER TIME DIFFERENCE.
        
        !@ COMPUTE LEFT AND RIGHT STAGE ON EACH ELEMENT.
        CALL reconstruct3 (UL, U, n,  1)
        CALL reconstruct3 (UR, U, n,  -1)
        !CALL MUSCLI(N, U, UL,UR,3)
        !CALL MUSCLE(N, U, UL, UR, 3)
        
        !@ COMPUTE REIMANN PROBLEM WITH LEFT AND RIGH STAGE.
        DO I = 1, N
            CALL RIEMANN_FLUX(UL(I), UR(I-1), F(I))
            !CALL RIEMANN_FLUX(UL(I), UR(I-1), F(I))
        END DO
        
        !@ BOUNDARY CONDITION FOR NUMERICAL FLUX.
        F(0) = F(1);    F(N) = F(N-1)
        
        !@ COMPUTE RIGHT-HAND-SIDE.
        CALL rhside(fp, f, n, hi)
        
        !@ UPDATE NEW SOLUTION BASED ON SECOND-ORDER TIME DIFFERENCE.
        UNEW(:) = (2.D0*DT*FP - UNM1(:) + 4.D0*U(:))/3.D0
        
        !@ APPLY BOUNDARY CONDITION WITH NEW SOLUTION.
        CALL BOUNDARY_CONDITION(N, UNEW)
        
        !@ COPY SOLUTION
        UNM1 = U
        U = UNEW
        
    ELSE    !@ SECHEME => 1 EXPLICIT RUNGEKUTTA THREE STAGE.

!   Substep 1
    !CALL flux (fp, fm, u, n)
    !CALL reconstruct5 (frp, fp, n,  1)
    !CALL reconstruct5 (frm, fm, n, -1)
    !f = frp + frm
    CALL reconstruct5 (UL, U, n,  1)
    CALL reconstruct5 (UR, U, n,  -1)
    
    DO I = 1, N
        CALL RIEMANN_FLUX(UL(I), UR(I-1), F(I))
    END DO
    F(0) = F(1)
    F(N) = F(N-1)
    CALL rhside(fp, f, n, hi)
!    CALL VISCOUS_FUNCTION (fp, u, n, hi, nu)
    up(:,1) = u + dt * fp
    CALL BOUNDARY_CONDITION(N, U)
!
!   Substep 2
   !CALL flux (fp, fm, up(:,1), n)
   !CALL reconstruct5 (frp, fp, n,  1)
   !CALL reconstruct5 (frm, fm, n, -1)
   !f = frp + frm
   CALL reconstruct5 (UL, UP(:,1), n,  1)
   CALL reconstruct5 (UR, UP(:,1), n,  -1)
   
   DO I = 1, N
    CALL RIEMANN_FLUX(UL(I), UR(I-1), F(I))
   END DO  
   F(0) = F(1)
   F(N) = F(N-1)    
    CALL rhside (fp, f, n, hi)
!    CALL VISCOUS_FUNCTION (fp, up(:,1), n, hi, nu)
    up(:,2) = rk(1,1) * u + rk(1,2) * (up(:,1) + dt * fp)
    CALL BOUNDARY_CONDITION(N, U)
!
!   Substep 3
    !CALL flux (fp, fm, up(:,2), n)
    !CALL reconstruct5 (frp, fp, n,  1)
    !CALL reconstruct5 (frm, fm, n, -1)
    !f = frp + frm
    CALL reconstruct5 (UL, UP(:,2), n,  1)
    CALL reconstruct5 (UR, UP(:,2), n,  -1)
    
    DO I = 1, N
        CALL RIEMANN_FLUX(UL(I), UR(I-1), F(I))
    END DO  
    F(0) = F(1)
    F(N) = F(N-1)    
    CALL rhside(fp, f, n, hi)
!    CALL VISCOUS_FUNCTION (fp, up(:,2), n, hi, nu)
    u = rk(2,1) * u + rk(2,2) * (up(:,2) + dt * fp)
    CALL BOUNDARY_CONDITION(N, U)
    
    END IF
    !DO I = 0, N
    !    IF (U(I) >=1 )U(I) = 1.
    !    IF (U(I) <=0)U(I) = 0.
    !END DO
    
!
!   Write the solution at the current time step
!    CALL h5sselect_hyperslab_f(filespace_hid, H5S_SELECT_SET_F, &
!         [INTEGER(HSIZE_T)::0, nt], [INTEGER(HSIZE_T)::n, 1], error)
!    CALL h5dwrite_f(dataset_hid, H5T_NATIVE_DOUBLE, u(0:n-1), &
!         [INTEGER(HSIZE_T)::n, 1], error, memspace_hid, filespace_hid)
    write(*, *) "INTERMEDIATE: "
    open(unit = 10, file = "INTERMEDIATE.TEC")
    do i = 0, n, 1
        write(10, *) x(i), u(i)
    end do
    
    close(unit = 10)
  END DO timeloop

  ! Tear down resources
  DEALLOCATE (x, xh, u, up, frp, frm, f, fp, fm)
!  CALL h5sclose_f(memspace_hid, error)
!  CALL h5dclose_f(dataset_hid, error)
!  CALL h5sclose_f(dataspace_hid, error)
!  CALL h5fclose_f(file_hid, error)
!  CALL h5close_f(error)

END PROGRAM main
   
! $HeadURL$
! $Id$

!> Compute the Lax-Friedrichs flux
!! \f$\hat{f}^{\mbox{LF}}\left(u^{-},u^{+}\right)\f$ given values of
!! \f$u\f$.  See section 3.1 of Shu's 2009 SIAM Review paper
!! or section 2 of Liu, Sher, and Chan's 1994 JCP paper for more details.
!!
!! @param fp
!! @param fm
!! @param u  \f$u\f$
!! @param n  Grid size
SUBROUTINE flux (fp, fm, u, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: n
  REAL(8), INTENT(IN)  :: u(0:n)
  REAL(8), INTENT(OUT) :: fp(0:n), fm(0:n)
  REAL(8)              :: alpha

! Euler equation where f(u) := u**2/2
!  alpha = MAXVAL(ABS(u))
!  fp    = 0.5_dp * (0.5d0*u**2 + alpha*u)
!  fm    = 0.5_dp * (0.5d0*u**2 - alpha*u)

! The scalar advection equation where f(u) := u
 alpha = 1.d0
 fp    = 0.5d0 * (u + alpha*u)
 fm    = 0.5d0 * (u - alpha*u)

END SUBROUTINE flux 
    
SUBROUTINE RIEMANN_FLUX(UL, UR, F)
IMPLICIT NONE
REAL(8), INTENT(IN) :: UL, UR
REAL(8), INTENT(OUT)    :: F
INTEGER :: I

    F = 0.5*(UL + UR) - 0.5D0*(UR - UL)
RETURN 
END SUBROUTINE RIEMANN_FLUX
    
! $HeadURL$
! $Id$

!> Compute a 3rd-order WENO reconstruction.
!! Given function values \f$u(x_i)\f$ for \f$i\in\left\{0,\ldots,n\right\}\f$
!! compute the reconstruction at \f$u_{r}\left(x_{i+1/2}\right)\f$ following
!! section 3.5 of Liu, Osher, and Chan's 1994 JCP paper.
!!
!! @param ur    Reconstruction \f$u_{r}\left(x_{i+1/2}\right)\f$
!! @param u     Function values \f$u(x_i)\f$
!! @param n     Grid size
!! @param bias  If strictly positive, bias stencil to the left.
!!              Otherwise, bias stencil to the right.
SUBROUTINE reconstruct3 (ur, u, n, bias)
! Equation numbers in the implementation refer to Liu, Osher, and Chan's paper.
    USE IR_Precision
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, bias
  REAL(R8P), INTENT(IN) :: u(0:n)
  REAL(R8P), INTENT(OUT) :: ur(0:n)
  REAL(R8P), PARAMETER :: eps = 1.d-14  ! guarantee nonzero denominator
  REAL(R8P) :: beta(1:2)               ! smoothness indicators
  REAL(R8P) :: w(1:2)                  ! nonlinear weights
  REAL(R8P) :: wt(1:2), wtsumi         ! temporary nonlinear weights
  REAL(R8P) :: urloc(1:2)              ! the two local reconstructions
  REAL(R8P) :: a(1:2,1:2)              ! weights in reconstruction
  INTEGER :: i
  REAL(R8P) :: v(-1:n+2)               ! add on periodic BCs
  REAL(R8P) :: v0, vp, vm              ! local values

  a(1,1) = -1.d0 / 2.d0
  a(1,2) =  3.d0 / 2.d0
  a(2,1) =  1.d0 / 2.d0
  a(2,2) =  1.d0 / 2.d0

! Add on periodic boundary conditions
! this is wasteful but results in a single loop so the code is easier to read
  v(0:n) = u(0:n)
  v(-1)  = u(n-1)
  v(n+1:n+2) = u(1:2)

 IF (bias > 0) THEN ! Bias to the left, case 1 in section 3.5
   DO i = 0, n, 1
     v0 = v(i)
     vp = v(i+1)
     vm = v(i-1)
! The reconstructed values at x(i+1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
     urloc(1) = a(1,1) * vm + a(1,2) * v0
     urloc(2) = a(2,1) * v0 + a(2,2) * vp
! Smoothness indicators from p206 just above equation 3.16
     beta(1) = (v0 - vm)**2
     beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
     wt(1) = 0.5d0 / ((eps + beta(1))**2)
     wt(2) = 1.0d0 / ((eps + beta(2))**2)
     wtsumi = 1.d0 / (wt(1) + wt(2))
     w(1) = wt(1) * wtsumi
     w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
     ur(i) = w(1) * urloc(1) + w(2) * urloc(2)
   END DO
  ELSE ! biased to the right, case 2 in section 3.5
    DO i = 1, n+1, 1
      v0 = v(i)
      vp = v(i+1)
      vm = v(i-1)
! The reconstructed values at x(i-1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
      urloc(1) = a(2,1) * vm + a(2,2) * v0
      urloc(2) = a(1,2) * v0 + a(1,1) * vp
! Smoothness indicators from p206 just above equation 3.16
      beta(1) = (v0 - vm)**2
      beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
      wt(1) = 1.0d0 / ((eps + beta(1))**2)
      wt(2) = 0.5d0 / ((eps + beta(2))**2)
      wtsumi = 1.d0 / (wt(1) + wt(2))
      w(1) = wt(1) * wtsumi
      w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
      ur(i-1) = w(1) * urloc(1) + w(2) * urloc(2)
    END DO
  END IF
END SUBROUTINE reconstruct3   
    
    
! $HeadURL$
! $Id$

!> Compute a 5th-order WENO reconstruction.
!! Given function values \f$u(x_i)\f$ for \f$i\in\left\{0,\ldots,n\right\}\f$
!! compute the reconstruction at \f$u_{r}\left(x_{i+1/2}\right)\f$ following
!! Shu's 2009 SIAM Review paper.
!!
!! @param ur    Reconstruction \f$u_{r}\left(x_{i+1/2}\right)\f$
!! @param u     Function values \f$u(x_i)\f$
!! @param n     Grid size
!! @param bias  If strictly positive, bias stencil to the left.
!!              Otherwise, bias stencil to the right.
SUBROUTINE reconstruct5 (ur, u, n, bias)
! Equation numbers in the implementation refer to Shu's paper.
    USE IR_Precision
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, bias
  REAL(R8P), INTENT(IN) :: u(0:n)
  REAL(R8P), INTENT(OUT) :: ur(0:n)
  REAL(R8P), PARAMETER :: eps = 1.d-14  ! guarantee nonzero denominator
  REAL(R8P) :: beta(1:3)                ! smoothness indicators
  REAL(R8P) :: w(1:3)                   ! nonlinear weights
  REAL(R8P) :: wt(1:3), wtsumi          ! temporary nonlinear weights
  REAL(R8P) :: gam(1:3)                 ! linear weights
  REAL(R8P) :: urloc(1:3)               ! the three local reconstructions
  REAL(R8P) :: a(1:3,1:3)               ! weights in reconstruction
  REAL(R8P) :: b(1:2)                   ! constants for beta computation
  INTEGER(I4P) :: i
  REAL(8) :: v(-2:n+3)               ! add on periodic BCs
  REAL(8) :: v0, vp, vpp, vm, vmm    ! local values

  a(1,1) =  1._R8P / 3._R8P
  a(1,2) = -7._R8P / 6._R8P
  a(1,3) = 11._R8P / 6._R8P
  a(2,1) = -1._R8P / 6._R8P
  a(2,2) =  5._R8P / 6._R8P
  a(2,3) =  1._R8P / 3._R8P
  a(3,1) =  1._R8P / 3._R8P
  a(3,2) =  5._R8P / 6._R8P
  a(3,3) = -1._R8P / 6._R8P

  b(1) = 13._R8P / 12._R8P
  b(2) =  1._R8P /  4._R8P
! just below (2.15)
  gam(1) = 1._R8P / 10._R8P
  gam(2) = 3._R8P /  5._R8P
  gam(3) = 3._R8P / 10._R8P

! add on periodic boundary condition
! this is wasteful but results in a single loop so the code is easier to read
  v(0:n)     = u(0:n)
  v(-2:-1)   = u(n-2:n-1)
  v(n+1:n+3) = u(1:3)

  IF (bias > 0) THEN ! Bias to the left
    DO i = 0, n, 1
      v0  = v(i)
      vp  = v(i+1)
      vpp = v(i+2)
      vm  = v(i-1)
      vmm = v(i-2)
! The three reconstructed values at x(i+1/2)
! Formulas (2.11), (2.12), (2.13)
      urloc(1) = a(1,1) * vmm + a(1,2) * vm + a(1,3) * v0
      urloc(2) = a(2,1) * vm + a(2,2) * v0 + a(2,3) * vp
      urloc(3) = a(3,1) * v0 + a(3,2) * vp + a(3,3) * vpp
! Smoothness indicators, formula (2.17)
      beta(1) = b(1) * (vmm - 2._R8P * vm + v0)**2 + b(2) * (vmm - 4._R8P * vm + 3._R8P * v0)**2
      beta(2) = b(1) * (vm - 2._R8P * v0 + vp)**2 + b(2) * (vm - vp)**2
      beta(3) = b(1) * (v0 - 2._R8P * vp + vpp)**2 + b(2) * (3._R8P * v0 - 4._R8P * vp + vpp)**2
! Compute nonlinear weights (2.10)
      wt(1) = gam(1) / ((eps + beta(1))**2)
      wt(2) = gam(2) / ((eps + beta(2))**2)
      wt(3) = gam(3) / ((eps + beta(3))**2)
      wtsumi = 1._R8P / (wt(1) + wt(2) + wt(3))
      w(1) = wt(1) * wtsumi
      w(2) = wt(2) * wtsumi
      w(3) = wt(3) * wtsumi
! Finally reconstruct, formula (2.16)
      ur(i) = w(1) * urloc(1) + w(2) * urloc(2) + w(3) * urloc(3)
    END DO
  ELSE ! biased to the right
    DO i = 1, n+1, 1
       v0 = v(i  )
       vp = v(i+1)
      vpp = v(i+2)
       vm = v(i-1)
      vmm = v(i-2)
! The three reconstructed values at x(i-1/2)
! Slightly different formulas than (2.11), (2.12), (2.13)
      urloc(1) = a(2,1) * vmm + a(2,2) * vm + a(2,3) * v0
      urloc(2) = a(3,1) * vm  + a(3,2) * v0 + a(3,3) * vp
      urloc(3) = a(1,3) * v0  + a(1,2) * vp + a(1,1) * vpp
! Smoothness indicators, formula (2.17)
      beta(1) = b(1) * (vmm - 2._R8P * vm + v0 )**2 + b(2) *(vmm - 4._R8P * vm + 3._R8P * v0)**2
      beta(2) = b(1) * ( vm - 2._R8P * v0 + vp )**2 + b(2) *( vm - vp)**2
      beta(3) = b(1) * ( v0 - 2._R8P * vp + vpp)**2 + b(2) *(3._R8P * v0 - 4._R8P * vp + vpp)**2
! Compute nonlinear weights (2.10)
      wt(1)  = gam(3) / ((eps + beta(1))**2)
      wt(2)  = gam(2) / ((eps + beta(2))**2)
      wt(3)  = gam(1) / ((eps + beta(3))**2)
      wtsumi = 1._R8P / (wt(1) + wt(2) + wt(3))
      w(1)   = wt(1) * wtsumi
      w(2)   = wt(2) * wtsumi
      w(3)   = wt(3) * wtsumi
! Finally reconstruct! Formula (2.16)
      ur(i-1) = w(1) * urloc(1) + w(2) * urloc(2) + w(3) * urloc(3)
    END DO
  END IF
END SUBROUTINE reconstruct5        
    
    
! $HeadURL$
! $Id$

!> Evaluate the right hand side of the inviscid equation's time evolution
!! \f$
!! \partial_{t} \bar{u}_j = -\frac{1}{h}\left[
!!   f\left(u\left(x_{j+1/2},t\right)\right)
!!   -
!!   f\left(u\left(x_{j+1/2},t\right)\right)
!! \right]\f$
!! assuming periodic boundary conditions.  See section 2 of Liu,
!! Osher, and Chan's 1994 JCP paper for more details.
!!
!! @param fout The evaluated right hand side
!! @param fin  The input data \f$f\left(u\left(x_{j+1/2}\right)\right)\f$
!!             for \f$j\in\left\{0,\dots,n\right\}\f$.  Usually
!!             this will be an approximation found through reconstruction.
!! @param n    Grid size
!! @param hi   \f$\frac{1}{h}\f$
    SUBROUTINE rhside (fout, fin, n, hi)
    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(IN) :: fin(0:n), hi
      REAL(8), INTENT(OUT) :: fout(0:n)
    
      fout(1:n) = -hi * (fin(1:n) - fin(0:n-1))
      fout(0)   = fout(n)
    END SUBROUTINE rhside  
    
    SUBROUTINE BOUNDARY_CONDITION(N,U)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL(8), DIMENSION(0:N), INTENT(INOUT)  :: U
    
    U(0) = 1.
    RETURN 
    END SUBROUTINE 

    
	SUBROUTINE MUSCLI(IMAX, AL, QRB,QLB,LIMITER)
	IMPLICIT NONE   
    REAL(8), DIMENSION(1:IMAX), INTENT(IN)  :: AL
	REAL(8), INTENT(INOUT) :: QRB(1:IMAX),QLB(1:IMAX)
	REAL(8) BT,KAPPA,MINMOD,ALBADA
	REAL(8) DELPRM,DELM,DERP
	REAL(8) DELB,DELBB,MUSCL_LMTER
	INTEGER I,J,K,LIMITER,II
    INTEGER, INTENT(IN) :: IMAX
    REAL(8) :: WWI = 1.0D0
	KAPPA=1.D0/3.D0
	
	DO 100 I = 2, IMAX, 1
	    IF (I.EQ.2) THEN
	        DELM = 0.D0
	    ELSE
	        DELM = AL(I-1)-AL(I-2)
	    ENDIF

        IF (I.EQ.IMAX)THEN
            DERP = 0.D0
        ELSE
	        DERP = AL(I+1)-AL(I)
	    ENDIF
	    DELPRM = AL(I)-AL(I-1)

! CAL. QRIGHT

            DELB  = MUSCL_LMTER(DELPRM,DERP,LIMITER)
	        DELBB = MUSCL_LMTER(DERP,DELPRM,LIMITER)
	
	    QRB(I) = AL(I) - 0.25D0*WWI*((1.D0-KAPPA)*DELB+(1.D0+KAPPA)*DELBB)
	
! CAL. Q_LEFT

            DELB    =MUSCL_LMTER(DELPRM,DELM,LIMITER)
	        DELBB   =MUSCL_LMTER(DELM,DELPRM,LIMITER)

	    QLB(I) = AL(I-1)+0.25D0*WWI*((1.D0-KAPPA)*DELB+(1.D0+KAPPA)*DELBB)

    IF (I.EQ.2) QLB(I) = (AL(I-1)+AL(I))/2.D0

	IF (I.EQ.IMAX)  QRB(I) = (AL(I-1)+AL(I))/2.D0

	!IF (QLB(I).LE.0.D0) QLB(I) = 0.D0  ! CONSTRAINS
	!IF (QLB(I).GT.1.D0) QLB(I) = 1.D0  ! CONSTRAINS
	!IF (QRB(I).LE.0.D0) QRB(I) = 0.D0  ! CONSTRAINS
	!IF (QRB(I).GT.1.D0) QRB(I) = 1.D0  ! CONSTRAINS


100	CONTINUE
    !!$OMP END PARALLEL DO

      RETURN
    END	
    
    FUNCTION MUSCL_LMTER(X,Y,LMT)
	IMPLICIT NONE
	INTEGER LMT
	REAL*8 X,Y,SIGX,SIGY,BT,MUSCL_LMTER,R
	IF(LMT.EQ.1)THEN 
	    BT=4.D0
	    IF (X.GT.0.D0)THEN
	        SIGX=1.0D0
	    ELSE
	        SIGX=-1.0D0
	    ENDIF
	    IF (Y.GT.0.D0)THEN
	        SIGY=1.0D0
	    ELSE
	        SIGY=-1.0D0
	    ENDIF
	    MUSCL_LMTER=SIGX*DMAX1(0.D0,DMIN1(X*SIGY,SIGX*BT*Y))!MINMOD	    
	ELSE
	    IF(DABS(Y).LT.1.0D-14)THEN
            R =1.0D0
        ELSE
            R =X/Y
        ENDIF
        IF(LMT.EQ.2)THEN
            MUSCL_LMTER = (R+R*R)/(1.0D0+R*R) !ALBADA1
            MUSCL_LMTER =Y*MUSCL_LMTER
        ELSEIF(LMT.EQ.3)THEN
            MUSCL_LMTER = 2.D0*R/(1.0D0+R*R)  !ALBADA2
            MUSCL_LMTER =Y*MUSCL_LMTER
        ELSEIF(LMT.EQ.4)THEN
	        MUSCL_LMTER =DMAX1(0.D0,DMIN1(2.D0*R,1.D0),DMIN1(R,2.D0)) !SUPERBEE
            MUSCL_LMTER =Y*MUSCL_LMTER
        ELSEIF(LMT.EQ.5)THEN
	        MUSCL_LMTER =(DABS(R) +R) /(1.D0+DABS(R)) !VAN LEER
            MUSCL_LMTER =Y*MUSCL_LMTER
        ENDIF
	ENDIF
	RETURN
    END    
