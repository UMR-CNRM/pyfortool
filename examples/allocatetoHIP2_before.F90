!#PYFT transfo: --allocatetoHIP
!=== COPYRIGHT AND LICENSE STATEMENT ===
!
!  This file is part of the TensorProductMultigrid code.
!  
!  (c) The copyright relating to this work is owned jointly by the
!  Crown, Met Office and NERC [2014]. However, it has been created
!  with the help of the GungHo Consortium, whose members are identified
!  at https://puma.nerc.ac.uk/trac/GungHo/wiki .
!  
!  Main Developer: Eike Mueller
!  
!  TensorProductMultigrid is free software: you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public License as
!  published by the Free Software Foundation, either version 3 of the
!  License, or (at your option) any later version.
!  
!  TensorProductMultigrid is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!  
!  You should have received a copy of the GNU Lesser General Public License
!  along with TensorProductMultigrid (see files COPYING and COPYING.LESSER).
!  If not, see <http://www.gnu.org/licenses/>.
!
!=== COPYRIGHT AND LICENSE STATEMENT ===


!==================================================================
!
!  Discretisation module of the model problem
!
!
!  -omega2 * (d^2/dx^2 + d^2/dy^2 + lambda2 * d^2/dz^2 ) u
!                                                + delta u = RHS
!  [Cartesian]
!
!  or
!
!  -omega2 * (laplace_{2d} + lambda2/r^2 d/dr (r^2 d/dr)) u
!                                                + delta u = RHS
!  [Spherical]
!
!  We use a cell centered finite volume discretisation with
!  The equation is discretised either in a unit cube or on 1/6th
!  of a cubed sphere grid.
!
!  The vertical grid spacing is not necessarily uniform and can
!  be chosen by specifying the vertical grid in a vector.
!
!  The following boundary conditions are used:
!
!  * Dirichlet in the horizontal
!  * Neumann in the vertical
!
!  For delta = 0 the operator reduces to the Poisson operator.
!
!    Eike Mueller, University of Bath, Feb 2012
!
!==================================================================

module discretisation

  use parameters_mg
  use messages
  use datatypes
  use communication
#ifndef MNH
  use mpi
  use mpi, only NMNH_COMM_WORLD => MPI_COMM_WORLD
  use mpi, only MNH_STATUSES_IGNORE => MPI_STATUSES_IGNORE
#else
  use modd_mpif
  USE MODD_VAR_ll,      ONLY: NMNH_COMM_WORLD, MNH_STATUSES_IGNORE
#endif 

  implicit none

private

  type model_parameters
    real(kind=rl) :: omega2   ! omega^2
    real(kind=rl) :: lambda2  ! lambda^2
    real(kind=rl) :: delta    ! delta
  end type model_parameters

! --- Stencil ---
!

! Grid traversal direction in SOR
  integer, parameter :: DIRECTION_FORWARD = 1
  integer, parameter :: DIRECTION_BACKWARD = 2

! Ordering in SOR
  ! Lexicographic ordering
  integer, parameter :: ORDERING_LEX = 1
  ! Red-black ordering
  integer, parameter :: ORDERING_RB = 2

  type smoother_parameters
    ! smoother
    integer :: smoother
    ! relaxation parameter
    real(kind=rl) :: rho
    ! ordering of degrees of freedom
    integer :: ordering
  end type smoother_parameters

  ! Allowed smoothers
  integer, parameter :: SMOOTHER_LINE_SOR = 3
  integer, parameter :: SMOOTHER_LINE_SSOR = 4
  integer, parameter :: SMOOTHER_LINE_JAC = 6

  ! Number of levels
  integer :: nlev

  ! Grid parameters
  type(grid_parameters) :: grid_param

  ! Model parameters
  type(model_parameters) :: model_param

  ! Smoother parameters
  type(smoother_parameters) :: smoother_param

  ! Arrays for measuring the residual reduction
  real(kind=rl), allocatable :: log_resreduction(:)
  integer, allocatable :: nsmooth_total(:)

  ! Data structure for storing the vertical discretisation
  type vertical_coefficients
    real(kind=rl), pointer , contiguous :: a(:)
    real(kind=rl), pointer , contiguous :: b(:)
    real(kind=rl), pointer , contiguous :: c(:)
    real(kind=rl), pointer , contiguous :: d(:)
  end type vertical_coefficients

  ! Stoarge for vertical coefficients
  type(vertical_coefficients) :: vert_coeff

  type Temp_jacobi
     real(kind=rl), pointer , contiguous :: r(:)
     real(kind=rl), pointer , contiguous :: c(:), utmp(:)
     real(kind=rl), pointer , contiguous :: u0(:,:,:) 
     real(kind=rl), pointer , contiguous :: ut0(:,:,:)
     type(scalar3d) , pointer            :: Sr,Sc,Sut0,Sutmp
  end type Temp_jacobi

  integer , parameter :: max_lev = 128
  type (Temp_jacobi) , save , dimension(max_lev) :: Tjacobi

  real(kind=rl), pointer , contiguous , dimension(:) :: zt1d_discretisation
  INTEGER                                            :: zt1d_discretisation_size = 10000000
  INTEGER , parameter                                :: zt1d_discretisation_size_factor = 4 * 8
  INTEGER                                            :: nt1d_discretisation_top = 0 , nt1d_discretisation_max = 0
  
  INTEGER , ALLOCATABLE, DIMENSION (:)               :: nt1d_discretisation_pointer , nt1d_discretisation_pointer_ksize
  INTEGER , parameter                                :: zt1d_discretisation_pointer_size = 128
  INTEGER                                            :: nt1d_discretisation_pointer_top = 0 , nt1d_discretisation_pointer_max = 0
  
  logical                                            :: gfirst_call_zt1d_discretisation = .true.
  

public::discretisation_initialise_mnh
public::discretisation_initialise
public::discretisation_finalise
public::smooth_mnh
public::smooth
public::line_SOR
public::line_SSOR
public::line_jacobi_mnh
public::line_jacobi
public::calculate_residual_mnh
public::calculate_residual
public::apply_mnh
public::apply
public::model_parameters
public::smoother_parameters
public::volume_of_element
public::SMOOTHER_LINE_SOR
public::SMOOTHER_LINE_SSOR
public::SMOOTHER_LINE_JAC
public::DIRECTION_FORWARD
public::DIRECTION_BACKWARD
public::ORDERING_LEX
public::ORDERING_RB
public::zt1d_discretisation_init
public::zt1d_discretisation_allocate3d

contains

!==================================================================
! Initialise module
!==================================================================

  subroutine zt1d_discretisation_init(kiu,kju,kku)
    implicit none

    integer , optional , intent(in) :: kiu,kju,kku ! size of 1 3D MNH array

    if (gfirst_call_zt1d_discretisation) then
       !
       gfirst_call_zt1d_discretisation = .false.
       !
       if (present(kiu)) then
          zt1d_discretisation_size = kiu*kju*kku * zt1d_discretisation_size_factor
       end if
       allocate(zt1d_discretisation(zt1d_discretisation_size))
       zt1d_discretisation = 0.0
       !$acc enter data copyin(zt1d_discretisation)
       !
       allocate(nt1d_discretisation_pointer(zt1d_discretisation_pointer_size))
       nt1d_discretisation_pointer = 0.0
       allocate(nt1d_discretisation_pointer_ksize(zt1d_discretisation_pointer_size))
       nt1d_discretisation_pointer_ksize = 0.0
       !
    end if
    
  end subroutine zt1d_discretisation_init

  function zt1d_discretisation_allocate3d(ptab3d,kib,kie,kjb,kje,kkb,kke) result(kindex)
    implicit none
    real(kind=rl), pointer , contiguous , dimension(:,:,:) :: ptab3d
    !real(kind=rl), pointer , dimension(:,:,:) :: ptab3d
    integer , intent(in)                                   :: kib,kie,kjb,kje,kkb,kke
    integer                                                :: kindex

    ! local var

    integer :: ksize, kbeg, kend

    if ( nt1d_discretisation_pointer_top == zt1d_discretisation_pointer_size ) then
       print*,"ERROR zt1d_discretisation_allocate3d:: zt1d_discretisation_pointer_size to small=",zt1d_discretisation_pointer_size
       print*,"ERROR zt1d_discretisation_allocate3d:: Augment it "
       stop 
    end if

    ksize = (kie-kib+1)*(kje-kjb+1)*(kke-kkb+1)

    kbeg = nt1d_discretisation_top + 1
    kend  = nt1d_discretisation_top + ksize

    if ( kend > zt1d_discretisation_size ) then
       print*,"ERROR zt1d_discretisation_allocate3d:: zt1d_discretisation_size to small=",zt1d_discretisation_size
       print*,"ERROR zt1d_discretisation_allocate3d:: Augment it "
       stop       
    end if

    ptab3d(kib:kie,kjb:kje,kkb:kke) => zt1d_discretisation(kbeg:)

    nt1d_discretisation_pointer_top = nt1d_discretisation_pointer_top + 1
    nt1d_discretisation_top = kend
    
    nt1d_discretisation_pointer(nt1d_discretisation_pointer_top)       = kend
    nt1d_discretisation_pointer_ksize(nt1d_discretisation_pointer_top) = ksize

    if (  nt1d_discretisation_pointer_top > nt1d_discretisation_pointer_max ) then
       nt1d_discretisation_pointer_max =  nt1d_discretisation_pointer_top
!!$       print*,"zt1d_discretisation_allocate3d:: nt1d_discretisation_pointer_max=",nt1d_discretisation_pointer_max
!!$       call flush(6)
    endif
    
    if (  nt1d_discretisation_top > nt1d_discretisation_max ) then
       nt1d_discretisation_max =  nt1d_discretisation_top
!!$       print*,"zt1d_discretisation_allocate3d:: nt1d_discretisation_max=",nt1d_discretisation_max
!!$       call flush(6)
    endif

    kindex = nt1d_discretisation_pointer_top
  end function zt1d_discretisation_allocate3d

  
  
 subroutine discretisation_initialise_mnh(grid_param_in, &
                                       model_param_in, &
                                       smoother_param_in, &
                                       nlev_in, &
                                       PA_K,PB_K,PC_K,PD_K)
    implicit none
    type(grid_parameters), intent(in)   :: grid_param_in
    type(model_parameters), intent(in)   :: model_param_in
    type(smoother_parameters), intent(in)   :: smoother_param_in
    integer, intent(in) :: nlev_in
#ifndef MNH    
    real(kind=rl) , optional , intent (in) :: PA_K(:),PB_K(:),PC_K(:),PD_K(:)
#else    
    real(kind=MNH_REAL) , optional , intent (in) :: PA_K(:),PB_K(:),PC_K(:),PD_K(:)
#endif

    ! local var
    integer :: k

    grid_param = grid_param_in
    model_param = model_param_in
    smoother_param = smoother_param_in
    nlev = nlev_in
    allocate(log_resreduction(nlev))
    allocate(nsmooth_total(nlev))
    log_resreduction(:) = 0.0_rl
    nsmooth_total(:) = 0
    allocate(r_grid(grid_param%nz+1))
    if (grid_param%graded) then
      do k=1,grid_param%nz+1
        r_grid(k) = grid_param%H*(1.0_rl*(k-1.0_rl)/grid_param%nz)**2
      end do
    else
      do k=1,grid_param%nz+1
        r_grid(k) = grid_param%H*(1.0_rl*(k-1.0_rl)/grid_param%nz)
      end do
    end if
    ! Allocate arrays for vertical discretisation matrices
    ! and calculate matrix entries
    allocate(vert_coeff%a(grid_param%nz))
    allocate(vert_coeff%b(grid_param%nz))
    allocate(vert_coeff%c(grid_param%nz))
    allocate(vert_coeff%d(grid_param%nz))
    call construct_vertical_coeff_mnh(PA_K,PB_K,PC_K,PD_K)
    !$acc enter data copyin(vert_coeff%a,vert_coeff%b,vert_coeff%c,vert_coeff%d)
  end subroutine discretisation_initialise_mnh

  subroutine discretisation_initialise(grid_param_in, &
                                       model_param_in, &
                                       smoother_param_in, &
                                       nlev_in)
    implicit none
    type(grid_parameters), intent(in)   :: grid_param_in
    type(model_parameters), intent(in)   :: model_param_in
    type(smoother_parameters), intent(in)   :: smoother_param_in
    integer, intent(in) :: nlev_in
    integer :: k
    grid_param = grid_param_in
    model_param = model_param_in
    smoother_param = smoother_param_in
    nlev = nlev_in
    allocate(log_resreduction(nlev))
    allocate(nsmooth_total(nlev))
    log_resreduction(:) = 0.0_rl
    nsmooth_total(:) = 0
    allocate(r_grid(grid_param%nz+1))
    if (grid_param%graded) then
      do k=1,grid_param%nz+1
        r_grid(k) = grid_param%H*(1.0_rl*(k-1.0_rl)/grid_param%nz)**2
      end do
    else
      do k=1,grid_param%nz+1
        r_grid(k) = grid_param%H*(1.0_rl*(k-1.0_rl)/grid_param%nz)
      end do
    end if
#ifdef CARTESIANGEOMETRY
#else
    r_grid(:) = 1.0_rl + r_grid(:)
#endif
    ! Allocate arrays for vertical discretisation matrices
    ! and calculate matrix entries
    allocate(vert_coeff%a(grid_param%nz))
    allocate(vert_coeff%b(grid_param%nz))
    allocate(vert_coeff%c(grid_param%nz))
    allocate(vert_coeff%d(grid_param%nz))
    call construct_vertical_coeff()
  end subroutine discretisation_initialise

!==================================================================
! Finalise module
!==================================================================
  subroutine discretisation_finalise()
    implicit none
    integer :: level
    real(kind=rl) :: rho_avg
#ifdef MEASURESMOOTHINGRATE
    if (i_am_master_mpi) then
      write(STDOUT,'("Average smoothing rates:")')
      do level=nlev,1,-1
        if (nsmooth_total(level) > 0) then
          rho_avg = exp(log_resreduction(level)/nsmooth_total(level))
        else
          rho_avg = 1.0_rl
        end if
        write(STDOUT,'("rho_{avg}(",I3,") = ",E10.4," ( ",I5," x )")') &
          level, rho_avg, nsmooth_total(level)
      end do
    end if
#endif
    deallocate(log_resreduction)
    deallocate(nsmooth_total)
    deallocate(r_grid)
    ! Deallocate storage for vertical discretisation matrices
    deallocate(vert_coeff%a)
    deallocate(vert_coeff%b)
    deallocate(vert_coeff%c)
    deallocate(vert_coeff%d)
  end subroutine discretisation_finalise

!==================================================================
! Construct alpha_{i',j'} and |T_{ij}| needed for the
! horizontal stencil
! ( alpha_{i+1,j},
!   alpha_{i-1,j},
!   alpha_{i,j+1},
!   alpha_{i,j-1},
!   alpha_{ij})
! (ix,iy) are LOCAL indices of the grid boxes, which are
! converted to global indices
!==================================================================
  subroutine construct_alpha_T_mnh(grid_param,ix,iy,alpha_T,Tij)
    implicit none
    type(grid_parameters), intent(in) :: grid_param
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    real(kind=rl), intent(inout), dimension(5) :: alpha_T
    real(kind=rl), intent(out) :: Tij

    !local var
    real(kind=rl)              :: h, rho_i, sigma_j
    real(kind=rl)              :: xcoef 
    logical                    :: l2nd       

    h = grid_param%L/grid_param%n
    ! Cartesian coefficients
    Tij = h**2
    ! optimisation for newman MNH case = all coef constant
    alpha_T(1:4) = 1.0_rl
    alpha_T(5) = 4.0_rl
    return
    xcoef = 0.5_rl ! 0.0
    l2nd = .false. ! .true. ! .false.
    alpha_T(1) = 1.0
    alpha_T(2) = 1.0
    if (ix == grid_param%n) then
      alpha_T(1) = xcoef * 2.0_rl
      if (l2nd) alpha_T(2) = 2.0_rl
    end if
    if (ix == 1) then
      alpha_T(2) = xcoef * 2.0_rl
      if (l2nd) alpha_T(1) = 2.0_rl
    end if
    alpha_T(3) = 1.0
    alpha_T(4) = 1.0
    if (iy == grid_param%n) then
      alpha_T(3) = xcoef * 2.0_rl
      if (l2nd) alpha_T(4) = 2.0
    end if
    if (iy == 1) then
      alpha_T(4) = xcoef * 2.0_rl
      if (l2nd) alpha_T(3) = 2.0
    end if

    alpha_T(5) = alpha_T(1) + alpha_T(2) + alpha_T(3) + alpha_T(4)
  end subroutine construct_alpha_T_mnh
! constant coef for MNH
  subroutine construct_alpha_T_cst_mnh(grid_param,alpha_T,Tij)
    implicit none
    type(grid_parameters), intent(in) :: grid_param
    real(kind=rl), intent(inout), dimension(5) :: alpha_T
    real(kind=rl), intent(out) :: Tij

    !local var
    real(kind=rl)              :: h

    h = grid_param%L/grid_param%n
    ! Cartesian coefficients
    Tij = h**2
    ! optimisation for newman MNH case = all coef constant
    alpha_T(1:4) = 1.0_rl
    alpha_T(5) = 4.0_rl
 
  end subroutine construct_alpha_T_cst_mnh
!==================================================================
  subroutine construct_alpha_T(grid_param,ix,iy,alpha_T,Tij)
    implicit none
    type(grid_parameters), intent(in) :: grid_param
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    real(kind=rl), intent(inout), dimension(5) :: alpha_T
    real(kind=rl), intent(out) :: Tij
    real(kind=rl)              :: h, rho_i, sigma_j
#ifdef CARTESIANGEOMETRY
    h = grid_param%L/grid_param%n
    ! Cartesian coefficients
    Tij = h**2
    if (ix == grid_param%n) then
      alpha_T(1) = 2.0_rl  
    else
      alpha_T(1) = 1.0_rl
    end if
    if (ix == 1) then
      alpha_T(2) = 2.0_rl  
    else
      alpha_T(2) = 1.0_rl
    end if
    if (iy == grid_param%n) then
      alpha_T(3) = 2.0_rl  
    else
      alpha_T(3) = 1.0_rl
    end if
    if (iy == 1) then
      alpha_T(4) = 2.0_rl  
    else
      alpha_T(4) = 1.0_rl
    end if
#else
    ! Coefficients in cubed sphere geometry
    ! (rho_i,sigma_j) \in [-1,1] x [-1,1] are the coordinates of the
    ! cubed sphere segment
    h = 2.0_rl/grid_param%n
    Tij = volume_of_element(ix,iy,grid_param)
    rho_i = 2.0_rl*(1.0_rl*ix-0.5_rl)/grid_param%n-1.0_rl
    sigma_j = 2.0_rl*(1.0_rl*iy-0.5_rl)/grid_param%n-1.0_rl
    ! alpha_{i+1,j}
    if (ix == grid_param%n) then
      alpha_T(1) = 2.0_rl*SQRT((1.0_rl+(rho_i+0.25_rl*h)**2)/(1.0_rl+sigma_j**2))
    else
      alpha_T(1) = SQRT((1.0_rl+(rho_i+0.5_rl*h)**2)/(1.0_rl+sigma_j**2))
    end if
    ! alpha_{i-1,j}
    if (ix == 1) then
      alpha_T(2) = 2.0_rl*SQRT((1.0_rl+(rho_i-0.25_rl*h)**2)/(1.0_rl+sigma_j**2))
    else
      alpha_T(2) = SQRT((1.0_rl+(rho_i-0.5_rl*h)**2)/(1.0_rl+sigma_j**2))
    end if
    ! alpha_{i,j+1}
    if (iy == grid_param%n) then
      alpha_T(3) = 2.0_rl*SQRT((1.0_rl+(sigma_j+0.25_rl*h)**2)/(1.0_rl+rho_i**2))
    else
      alpha_T(3) = SQRT((1.0_rl+(sigma_j+0.5_rl*h)**2)/(1.0_rl+rho_i**2))
    end if
    ! alpha_{i,j-1}
    if (iy == 1) then
      alpha_T(4) = 2.0_rl*SQRT((1.0_rl+(sigma_j-0.25_rl*h)**2)/(1.0_rl+rho_i**2))
    else
      alpha_T(4) = SQRT((1.0_rl+(sigma_j-0.5_rl*h)**2)/(1.0_rl+rho_i**2))
    end if
#endif
    alpha_T(5) = alpha_T(1) + alpha_T(2) + alpha_T(3) + alpha_T(4)
  end subroutine construct_alpha_T
!==================================================================
! Construct coefficients of tridiagonal matrix A_T
! describing the coupling in the vertical direction and the
! diagonal matrix diag(d)
!==================================================================
subroutine construct_vertical_coeff_mnh(PA_K,PB_K,PC_K,PD_K)
  implicit none
#ifndef MNH    
  real(kind=rl) , optional , intent (in) :: PA_K(:),PB_K(:),PC_K(:),PD_K(:)
#else    
  real(kind=MNH_REAL) , optional , intent (in) :: PA_K(:),PB_K(:),PC_K(:),PD_K(:)
#endif
  !local var
  real(kind=rl) :: a_k_tmp, b_k_tmp, c_k_tmp, d_k_tmp
  real(kind=rl) :: omega2, lambda2, delta, vol_r, surface_k, surface_kp1
  integer :: k
  
  IF (.NOT. PRESENT(PA_K)) THEN
  omega2 = model_param%omega2
  lambda2 = model_param%lambda2
  delta = model_param%delta
  do k = 1, grid_param%nz

    vol_r = r_grid(k+1)-r_grid(k)
    surface_k = 1.0_rl
    surface_kp1 = 1.0_rl

    ! Diagonal element
    a_k_tmp = delta*vol_r
    ! off diagonal elements
    ! Boundary conditions
    ! Top
    if (k == grid_param%nz) then
      if (grid_param%vertbc == VERTBC_DIRICHLET) then
        b_k_tmp = - 2.0_rl * omega2*lambda2 &
                           * surface_kp1/(r_grid(k+1)-r_grid(k))
      else
        b_k_tmp = 0.0_rl
      end if
    else
      b_k_tmp = - 2.0_rl*omega2*lambda2 &
              * surface_kp1/(r_grid(k+2)-r_grid(k))
    end if
    ! Bottom
    if (k == 1) then
      if (grid_param%vertbc == VERTBC_DIRICHLET) then
        c_k_tmp = - 2.0_rl * omega2*lambda2 &
                           * surface_k/(r_grid(k+1)-r_grid(k))
      else
        c_k_tmp = 0.0_rl
      end if
    else
      c_k_tmp = - 2.0_rl * omega2 * lambda2 &
              * surface_k/(r_grid(k+1)-r_grid(k-1))
    end if
    ! Diagonal matrix d_k
    d_k_tmp = - omega2 * (r_grid(k+1)-r_grid(k))
    vert_coeff%a(k) = a_k_tmp/d_k_tmp
    vert_coeff%b(k) = b_k_tmp/d_k_tmp
    vert_coeff%c(k) = c_k_tmp/d_k_tmp
    vert_coeff%d(k) = d_k_tmp
  end do
  ELSE
  do k = 1, grid_param%nz
    vert_coeff%a(k) = PA_K(k)
    vert_coeff%b(k) = PB_K(k)
    vert_coeff%c(k) = PC_K(k)
    vert_coeff%d(k) = PD_K(k)
  end do
  ENDIF
end subroutine construct_vertical_coeff_mnh

subroutine construct_vertical_coeff()
  implicit none
  real(kind=rl) :: a_k_tmp, b_k_tmp, c_k_tmp, d_k_tmp
  real(kind=rl) :: omega2, lambda2, delta, vol_r, surface_k, surface_kp1
  integer :: k
  omega2 = model_param%omega2
  lambda2 = model_param%lambda2
  delta = model_param%delta
  do k = 1, grid_param%nz
#ifdef  CARTESIANGEOMETRY
    vol_r = r_grid(k+1)-r_grid(k)
    surface_k = 1.0_rl
    surface_kp1 = 1.0_rl
#else
    vol_r = (r_grid(k+1)**3 - r_grid(k)**3)/3.0_rl
    surface_k = r_grid(k)**2
    surface_kp1 = r_grid(k+1)**2
#endif
    ! Diagonal element
    a_k_tmp = delta*vol_r
    ! off diagonal elements
    ! Boundary conditions
    ! Top
    if (k == grid_param%nz) then
      if (grid_param%vertbc == VERTBC_DIRICHLET) then
        b_k_tmp = - 2.0_rl * omega2*lambda2 &
                           * surface_kp1/(r_grid(k+1)-r_grid(k))
      else
        b_k_tmp = 0.0_rl
      end if
    else
      b_k_tmp = - 2.0_rl*omega2*lambda2 &
              * surface_kp1/(r_grid(k+2)-r_grid(k))
    end if
    ! Bottom
    if (k == 1) then
      if (grid_param%vertbc == VERTBC_DIRICHLET) then
        c_k_tmp = - 2.0_rl * omega2*lambda2 &
                           * surface_k/(r_grid(k+1)-r_grid(k))
      else
        c_k_tmp = 0.0_rl
      end if
    else
      c_k_tmp = - 2.0_rl * omega2 * lambda2 &
              * surface_k/(r_grid(k+1)-r_grid(k-1))
    end if
    ! Diagonal matrix d_k
    d_k_tmp = - omega2 * (r_grid(k+1)-r_grid(k))
    vert_coeff%a(k) = a_k_tmp/d_k_tmp
    vert_coeff%b(k) = b_k_tmp/d_k_tmp
    vert_coeff%c(k) = c_k_tmp/d_k_tmp
    vert_coeff%d(k) = d_k_tmp
  end do
end subroutine construct_vertical_coeff

!==================================================================
! Calculate local residual r = b - A.u
!==================================================================
  subroutine calculate_residual_mnh(level,m,b,u,r)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    type(scalar3d), intent(inout)      :: r
    integer :: ix,iy,iz
    integer ::  iib,iie,ijb,ije,ikb,ike

    real(kind=rl) , dimension(:,:,:) , pointer , contiguous :: zr_st , zb_st

    integer :: nib,nie,njb,nje,nzb,nze

    ! r <- A.u
    !call boundary_mnh(u)
    call apply_mnh(u,r)
    ! r <- b - r = b - A.u
    if (LUseO) then
       do ix=u%icompx_min,u%icompx_max
          do iy=u%icompy_min,u%icompy_max
             do iz=1,u%grid_param%nz
                r%s(iz,iy,ix) = b%s(iz,iy,ix) - r%s(iz,iy,ix)
             end do
          end do
       end do
    endif
    if (LUseT) then
!!$       do iz=1,u%grid_param%nz
!!$          do iy=u%icompy_min,u%icompy_max
!!$             do ix=u%icompx_min,u%icompx_max
!!$                r%st(ix,iy,iz) = b%st(ix,iy,iz) - r%st(ix,iy,iz)
!!$             end do
!!$          end do
!!$       end do
       !-----------------
       iib=u%icompx_min
       iie=u%icompx_max
       ijb=u%icompy_min
       ije=u%icompy_max
       ikb=1
       ike=u%grid_param%nz

       zr_st => r%st
       zb_st => b%st

       nib = Lbound(zr_st,1) ; nie = Ubound(zr_st,1)
       njb = Lbound(zr_st,2) ; nje = Ubound(zr_st,2)
       nzb = Lbound(zr_st,3) ; nze = Ubound(zr_st,3)

       call calculate_residual_mnh_dim(zr_st,zb_st)

    endif

  contains

    subroutine calculate_residual_mnh_dim(pzr_st,pzb_st)
      implicit none

      real(kind=rl) :: pzr_st(nib:nie,njb:nje,nzb:nze), &
              pzb_st(nib:nie,njb:nje,nzb:nze)

      !$acc kernels present(pzr_st,pzb_st)
        pzr_st(iib:iie,ijb:ije,ikb:ike) = pzb_st(iib:iie,ijb:ije,ikb:ike) - pzr_st(iib:iie,ijb:ije,ikb:ike)
      !$acc end kernels
      
    end subroutine calculate_residual_mnh_dim

  end subroutine calculate_residual_mnh

  subroutine calculate_residual(level,m,b,u,r)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    type(scalar3d), intent(inout)      :: r
    integer :: ix,iy,iz


    ! r <- A.u
    call apply(u,r)
    ! r <- b - r = b - A.u
    do ix=u%icompx_min,u%icompx_max
      do iy=u%icompy_min,u%icompy_max
        do iz=1,u%grid_param%nz
          r%s(iz,iy,ix) = b%s(iz,iy,ix) - r%s(iz,iy,ix)
        end do
      end do
    end do
  end subroutine calculate_residual

!==================================================================
! Apply operator v = A.u
!==================================================================
  subroutine apply_mnh(u,v)
    implicit none
    type(scalar3d), intent(inout)      :: u
    type(scalar3d), intent(inout)      :: v

    ! local var
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) :: Tij
    real(kind=rl) :: a_k, b_k, c_k, d_k
    integer :: ix,iy,iz
    real(kind=rl) :: tmp
    integer :: iib,iie,ijb,ije

    real(kind=rl), dimension(:,:,:) , pointer , contiguous :: zv_st , zu_st  
    real(kind=rl), dimension(:)     , pointer , contiguous :: za_k, zb_k, zc_k, zd_k
    integer :: ii,ij
    integer :: ize

    integer :: nib,nie,njb,nje,nzb,nze

    call boundary_mnh(u)

    if (LUseO) then     
       do ix=u%icompx_min,u%icompx_max
          do iy=u%icompy_min,u%icompy_max
             ! Construct horizontal part of stencil
             call construct_alpha_T_mnh(u%grid_param,  &
                  ix+u%ix_min-1, &
                  iy+u%iy_min-1, &
                  alpha_T,Tij)
             do iz=1,u%grid_param%nz
                a_k = vert_coeff%a(iz)
                b_k = vert_coeff%b(iz)
                c_k = vert_coeff%c(iz)
                d_k = vert_coeff%d(iz)
                tmp = ((a_k-b_k-c_k)*Tij ) * u%s(iz,iy,ix)
                if (iz < grid_param%nz) then
                   tmp = tmp + b_k*Tij * u%s(iz+1,iy,ix)
                end if
                if (iz > 1) then
                   tmp = tmp + c_k*Tij * u%s(iz-1,iy,ix)
                end if
                if ((iz > 1) .and. (iz < grid_param%nz)) then
                   tmp = tmp - alpha_T(5) * u%s(iz,  iy  ,ix  ) &
                        + alpha_T(1) * u%s(iz,  iy  ,ix+1) &
                        + alpha_T(2) * u%s(iz,  iy  ,ix-1) &
                        + alpha_T(3) * u%s(iz,  iy+1,ix  ) &
                        + alpha_T(4) * u%s(iz,  iy-1,ix  )
                end if
                v%s(iz,iy,ix) = d_k*tmp 
             end do
          end do
       end do
    endif
    if (LUseT) then 
       call construct_alpha_T_cst_mnh(u%grid_param,alpha_T,Tij)    
       !-----------------------------------------------------------
       iib=u%icompx_min
       iie=u%icompx_max
       ijb=u%icompy_min
       ije=u%icompy_max
       ize=u%grid_param%nz
       !
       zv_st => v%st
       zu_st => u%st
       zb_k => vert_coeff%b
       zc_k => vert_coeff%c
       zd_k => vert_coeff%d

       nib = Lbound(zv_st,1) ; nie = Ubound(zv_st,1)
       njb = Lbound(zv_st,2) ; nje = Ubound(zv_st,2)
       nzb = Lbound(zv_st,3) ; nze = Ubound(zv_st,3)

       call apply_mnh_dim(zv_st,zu_st,zb_k,zc_k,zd_k)
    
    endif

  contains

    subroutine apply_mnh_dim(pzv_st,pzu_st,pzb_k,pzc_k,pzd_k)

      implicit none

      real(kind=rl) :: pzv_st(nib:nie,njb:nje,nzb:nze), &
              pzu_st(nib:nie,njb:nje,nzb:nze), &
              pzb_k(ize),pzc_k(ize),pzd_k(ize)

      !$acc kernels present_cr(pzb_k,pzv_st,pzu_st)
      iz=1
      !$mnh_do_concurrent( ii=iib:iie , ij=ijb:ije )
         pzv_st(ii,ij,iz) = pzd_k(iz)* ( (-pzb_k(iz)-pzc_k(iz))*Tij * pzu_st(ii,ij,iz  )  &
                                          +pzb_k(iz)           *Tij * pzu_st(ii,ij,iz+1)  )
      !$mnh_end_do()
      
      !
      !$mnh_do_concurrent( ii=iib:iie , ij=ijb:ije , iz=2:ize-1 )
         pzv_st(ii,ij,iz) = pzd_k(iz)* ( ((-pzb_k(iz)-pzc_k(iz))*Tij - 4.0_rl ) * pzu_st(ii,ij,iz)   &
                                           +pzb_k(iz)           *Tij            * pzu_st(ii,ij,iz+1) &
                                                     +pzc_k(iz) *Tij            * pzu_st(ii,ij,iz-1) &
                                     +                                            pzu_st(ii+1,ij,iz) &
                                     +                                            pzu_st(ii-1,ij,iz) &
                                     +                                            pzu_st(ii,ij+1,iz) &
                                     +                                            pzu_st(ii,ij-1,iz) &
                                       )
       !$mnh_end_do()
       !
       iz=ize
       !$mnh_do_concurrent( ii=iib:iie , ij=ijb:ije )       
             pzv_st(ii,ij,iz) = pzd_k(iz)*  (  (-pzb_k(iz)-pzc_k(iz))*Tij  * pzu_st(ii,ij,iz)    &
                                                          +pzc_k(iz) *Tij  * pzu_st(ii,ij,iz-1)  )
       !$mnh_end_do()
       !$acc end kernels     
    
     end subroutine apply_mnh_dim

  end subroutine apply_mnh

  subroutine apply(u,v)
    implicit none
    type(scalar3d), intent(in)         :: u
    type(scalar3d), intent(inout)      :: v
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) :: Tij
    real(kind=rl) :: a_k, b_k, c_k, d_k
    integer :: ix,iy,iz
    real(kind=rl) :: tmp

    do ix=u%icompx_min,u%icompx_max
      do iy=u%icompy_min,u%icompy_max
        ! Construct horizontal part of stencil
        call construct_alpha_T(u%grid_param,  &
                               ix+u%ix_min-1, &
                               iy+u%iy_min-1, &
                               alpha_T,Tij)
        do iz=1,u%grid_param%nz
          a_k = vert_coeff%a(iz)
          b_k = vert_coeff%b(iz)
          c_k = vert_coeff%c(iz)
          d_k = vert_coeff%d(iz)
          tmp = ((a_k-b_k-c_k)*Tij - alpha_T(5)) * u%s(iz,iy,ix)
          if (iz < grid_param%nz) then
            tmp = tmp + b_k*Tij * u%s(iz+1,iy,ix)
          end if
          if (iz > 1) then
            tmp = tmp + c_k*Tij * u%s(iz-1,iy,ix)
          end if
          tmp = tmp + alpha_T(1) * u%s(iz,  iy  ,ix+1) &
                    + alpha_T(2) * u%s(iz,  iy  ,ix-1) &
                    + alpha_T(3) * u%s(iz,  iy+1,ix  ) &
                    + alpha_T(4) * u%s(iz,  iy-1,ix  )
          v%s(iz,iy,ix) = d_k*tmp
        end do
      end do
    end do
  end subroutine apply
!==================================================================
!==================================================================
!
!     S M O O T H E R S
!
!==================================================================
!==================================================================

!==================================================================
! Perform nsmooth smoother iterations
!==================================================================
  subroutine smooth_mnh(level,m,nsmooth,direction,b,u)
    implicit none
    integer, intent(in) :: level
    integer, intent(in) :: m
    integer, intent(in) :: nsmooth       ! Number of smoothing steps
    integer, intent(in) :: direction     ! Direction
    type(scalar3d), intent(inout) :: b   ! RHS
    type(scalar3d), intent(inout) :: u   ! solution vector
    integer :: i
    real(kind=rl) :: log_res_initial, log_res_final
    type(scalar3d) :: r
    integer :: halo_size
    integer :: nlocal, nz

    integer :: nlocalx,nlocaly

#ifdef MEASURESMOOTHINGRATE
    r%ix_min = u%ix_min
    r%ix_max = u%ix_max
    r%iy_min = u%iy_min
    r%iy_max = u%iy_max
    r%icompx_min = u%icompx_min
    r%icompx_max = u%icompx_max
    r%icompy_min = u%icompy_min
    r%icompy_max = u%icompy_max
    r%halo_size = u%halo_size
    r%isactive = u%isactive
    r%grid_param = u%grid_param
    nlocal = r%ix_max-r%ix_min+1
    nlocalx = r%ix_max-r%ix_min+1
    nlocaly = r%iy_max-r%iy_min+1
    halo_size = r%halo_size
    nz = r%grid_param%nz

    if (LUseO) then
       allocate(r%s(0:nz+1,                       &
                    1-halo_size:nlocal+halo_size, &
                    1-halo_size:nlocal+halo_size))
    end if

    if (LUseT) then
       allocate(r%st(1-halo_size:nlocalx+halo_size, &
                     1-halo_size:nlocaly+halo_size, &
                     0:nz+1))
       !$acc enter data create (r%st)
    end if

    call calculate_residual(level,m,b,u,r)
    log_res_initial = log(l2norm(r))
#endif
    ! Carry out nsmooth iterations of the smoother
    if (smoother_param%smoother == SMOOTHER_LINE_SOR) then
      do i=1,nsmooth
        call line_SOR_mnh(level,m,direction,b,u)
      end do
    else if (smoother_param%smoother == SMOOTHER_LINE_SSOR) then
      do i=1,nsmooth
        call line_SSOR_mnh(level,m,direction,b,u)
      end do
    else if (smoother_param%smoother == SMOOTHER_LINE_JAC) then
      do i=1,nsmooth
        call line_jacobi_mnh(level,m,b,u)
      end do
    end if
#ifdef MEASURESMOOTHINGRATE
    call calculate_residual_mnh(level,m,b,u,r)
    log_res_final = log(l2norm(r))
    log_resreduction(level) = log_resreduction(level) &
                            + (log_res_final - log_res_initial)
    nsmooth_total(level) = nsmooth_total(level) + nsmooth
    if (LUseO) deallocate(r%s)
    if (LUseT) then
       !$acc exit data delete(r%st)
       deallocate(r%st)       
    end if
#endif
  end subroutine smooth_mnh
!==================================================================
! Perform nsmooth smoother iterations
!==================================================================
  subroutine smooth(level,m,nsmooth,direction,b,u)
    implicit none
    integer, intent(in) :: level
    integer, intent(in) :: m
    integer, intent(in) :: nsmooth       ! Number of smoothing steps
    integer, intent(in) :: direction     ! Direction
    type(scalar3d), intent(inout) :: b   ! RHS
    type(scalar3d), intent(inout) :: u   ! solution vector
    integer :: i
    real(kind=rl) :: log_res_initial, log_res_final
    type(scalar3d) :: r
    integer :: halo_size
    integer :: nlocal, nz

#ifdef MEASURESMOOTHINGRATE
    r%ix_min = u%ix_min
    r%ix_max = u%ix_max
    r%iy_min = u%iy_min
    r%iy_max = u%iy_max
    r%icompx_min = u%icompx_min
    r%icompx_max = u%icompx_max
    r%icompy_min = u%icompy_min
    r%icompy_max = u%icompy_max
    r%halo_size = u%halo_size
    r%isactive = u%isactive
    r%grid_param = u%grid_param
    nlocal = r%ix_max-r%ix_min+1
    halo_size = r%halo_size
    nz = r%grid_param%nz
    if (LUseO) then
    allocate(r%s(0:nz+1,                       &
                 1-halo_size:nlocal+halo_size, &
                 1-halo_size:nlocal+halo_size))
    call calculate_residual(level,m,b,u,r)
    endif
    if (LUseT) then
    allocate(r%st(1-halo_size:nlocal+halo_size,   &
                  1-halo_size:nlocal+halo_size,   &
                  0:nz+1) )
    !$acc enter data create (r%st)
    endif
    log_res_initial = log(l2norm(r))
#endif
    ! Carry out nsmooth iterations of the smoother
    if (smoother_param%smoother == SMOOTHER_LINE_SOR) then
      do i=1,nsmooth
        call line_SOR(level,m,direction,b,u)
      end do
    else if (smoother_param%smoother == SMOOTHER_LINE_SSOR) then
      do i=1,nsmooth
        call line_SSOR(level,m,direction,b,u)
      end do
    else if (smoother_param%smoother == SMOOTHER_LINE_JAC) then
      do i=1,nsmooth
        call line_jacobi(level,m,b,u)
      end do
    end if
#ifdef MEASURESMOOTHINGRATE
    call calculate_residual(level,m,b,u,r)
    log_res_final = log(l2norm(r))
    log_resreduction(level) = log_resreduction(level) &
                            + (log_res_final - log_res_initial)
    nsmooth_total(level) = nsmooth_total(level) + nsmooth
    if (LUseO) deallocate(r%s)
    if (LUseT) then
       !$acc exit data delete(r%st)
       deallocate(r%st)
    end if
#endif
  end subroutine smooth
!==================================================================
! SOR line smoother mnh
!==================================================================
  subroutine line_SOR_mnh(level,m,direction,b,u)

    implicit none

    integer, intent(in)                :: level
    integer, intent(in)                :: m
    integer, intent(in)                :: direction
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u

    !local Var
    real(kind=rl), allocatable :: r(:)
    integer :: nz, nlocal
    real(kind=rl), allocatable :: c(:), utmp(:)
    integer :: ixmin(5), ixmax(5), dix
    integer :: iymin(5), iymax(5), diy
    integer :: color
    integer :: nsweeps, isweep
    integer :: ordering
    real(kind=rl) :: rho
    integer, dimension(4) :: send_requests, recv_requests
    integer, dimension(4) :: send_requestsT, recv_requestsT
    integer :: tmp, ierr
    integer :: iblock
    logical :: overlap_comms

    integer :: nlocalx,nlocaly

    call boundary_mnh(u)

    ordering = smoother_param%ordering
    rho = smoother_param%rho

    nz = u%grid_param%nz

    ! Create residual vector
    allocate(r(nz))
    ! Allocate memory for auxiliary vectors for Thomas algorithm
    allocate(c(nz))
    allocate(utmp(nz))
    nlocal = u%ix_max-u%ix_min+1
    nlocalx = u%ix_max-u%ix_min+1
    nlocaly = u%iy_max-u%iy_min+1
#ifdef OVERLAPCOMMS
    overlap_comms = (nlocal > 2)
#else
    overlap_comms = .false.
#endif
    ! Block 1 (N)
    ixmin(1) = 1
    ixmax(1) = nlocalx
    iymin(1) = 1
    iymax(1) = 1
    ! Block 2 (S)
    ixmin(2) = 1
    ixmax(2) = nlocalx
    iymin(2) = nlocaly
    iymax(2) = nlocaly
    ! Block 3 (W)
    ixmin(3) = 1
    ixmax(3) = 1
    iymin(3) = 2
    iymax(3) = nlocaly-1
    ! Block 4 (E)
    ixmin(4) = nlocalx
    ixmax(4) = nlocalx
    iymin(4) = 2
    iymax(4) = nlocaly-1
    ! Block 5 (INTERIOR)
    if (overlap_comms) then
      ixmin(5) = 2
      ixmax(5) = nlocalx-1
      iymin(5) = 2
      iymax(5) = nlocaly-1
    else
      ! If there are no interior cells, do not overlap
      ! communications and calculations, just loop over interior cells
      ixmin(5) = 1
      ixmax(5) = nlocalx
      iymin(5) = 1
      iymax(5) = nlocaly
    end if
    dix = +1
    diy = +1
    color = 1
    ! When iteration backwards over the grid, reverse the direction
    if (direction == DIRECTION_BACKWARD) then
      do iblock = 1, 5
        tmp = ixmax(iblock)
        ixmax(iblock) = ixmin(iblock)
        ixmin(iblock) = tmp
        tmp = iymax(iblock)
        iymax(iblock) = iymin(iblock)
        iymin(iblock) = tmp
      end do
      dix = -1
      diy = -1
      color = 0
    end if
    nsweeps = 1
    if (ordering == ORDERING_LEX) then
      nsweeps = 1
    else if (ordering == ORDERING_RB) then
      nsweeps = 2
    end if
    do isweep = 1, nsweeps
      if (overlap_comms) then
        ! Loop over cells next to boundary (iblock = 1,...,4)
        do iblock = 1, 4
          call loop_over_grid_mnh(iblock)
        end do
        ! Initiate halo exchange
        call ihaloswap_mnh(level,m,u,send_requests,recv_requests,send_requestsT,recv_requestsT)
      end if
      ! Loop over INTERIOR cells
      iblock = 5
      call loop_over_grid_mnh(iblock)
      if (overlap_comms) then
        if (m > 0) then
          if (LUseO) call mpi_waitall(4,recv_requests, MNH_STATUSES_IGNORE, ierr)
          if (LUseO) call mpi_waitall(4,send_requests, MNH_STATUSES_IGNORE, ierr)
          if (LUseT) call mpi_waitall(4,recv_requestsT, MNH_STATUSES_IGNORE, ierr)
          if (LUseT) call mpi_waitall(4,send_requestsT, MNH_STATUSES_IGNORE, ierr)
        end if
      else
        call haloswap_mnh(level,m,u)
      end if
      color = 1-color
    end do

    ! Free memory again
    deallocate(r)
    deallocate(c)
    deallocate(utmp)

    contains

    !------------------------------------------------------------------
    ! Loop over grid, for a given block
    !------------------------------------------------------------------
    subroutine loop_over_grid_mnh(iblock)
      implicit none
      integer, intent(in) :: iblock
      integer :: ix,iy,iz

      if (LUseO) then
         do ix=ixmin(iblock),ixmax(iblock),dix
            do iy=iymin(iblock),iymax(iblock),diy
               if (ordering == ORDERING_RB) then
                  if (mod((ix+u%ix_min)+(iy+u%iy_min),2) .ne. color) cycle
               end if
               call apply_tridiag_solve_mnh(ix,iy,r,c,b,         &
                    u%s(1:nz,iy  ,ix+1), &
                    u%s(1:nz,iy  ,ix-1), &
                    u%s(1:nz,iy+1,ix  ), &
                    u%s(1:nz,iy-1,ix  ), &
                    utmp)
               ! Add to field with overrelaxation-factor
               do iz=1,nz
                  u%s(iz,iy,ix) = (1.0_rl-rho)*u%s(iz,iy,ix) + rho*utmp(iz)
               end do
            end do
         end do
      end if
      if (LUseT) then
         do ix=ixmin(iblock),ixmax(iblock),dix
            do iy=iymin(iblock),iymax(iblock),diy
               if (ordering == ORDERING_RB) then
                  if (mod((ix+u%ix_min)+(iy+u%iy_min),2) .ne. color) cycle
               end if
               call apply_tridiag_solve_mnhT(ix,iy,r,c,b,         &
                    u%st(ix+1,iy  ,1:nz), &
                    u%st(ix-1,iy  ,1:nz), &
                    u%st(ix  ,iy+1,1:nz), &
                    u%st(ix  ,iy-1,1:nz), &
                    utmp)
               ! Add to field with overrelaxation-factor
               do iz=1,nz
                  u%st(ix,iy,iz) = (1.0_rl-rho)*u%st(ix,iy,iz) + rho*utmp(iz)
               end do
            end do
         end do
      end if   

    end subroutine loop_over_grid_mnh

  end subroutine line_SOR_mnh
!==================================================================
! SOR line smoother
!==================================================================
  subroutine line_SOR(level,m,direction,b,u)

    implicit none

    integer, intent(in)                :: level
    integer, intent(in)                :: m
    integer, intent(in)                :: direction
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    real(kind=rl), allocatable :: r(:)
    integer :: nz, nlocal
    real(kind=rl), allocatable :: c(:), utmp(:)
    integer :: ixmin(5), ixmax(5), dix
    integer :: iymin(5), iymax(5), diy
    integer :: color
    integer :: nsweeps, isweep
    integer :: ordering
    real(kind=rl) :: rho
    integer, dimension(4) :: send_requests, recv_requests
    integer :: tmp, ierr
    integer :: iblock
    logical :: overlap_comms

    ordering = smoother_param%ordering
    rho = smoother_param%rho

    nz = u%grid_param%nz

    ! Create residual vector
    allocate(r(nz))
    ! Allocate memory for auxiliary vectors for Thomas algorithm
    allocate(c(nz))
    allocate(utmp(nz))
    nlocal = u%ix_max-u%ix_min+1
#ifdef OVERLAPCOMMS
    overlap_comms = (nlocal > 2)
#else
    overlap_comms = .false.
#endif
    ! Block 1 (N)
    ixmin(1) = 1
    ixmax(1) = nlocal
    iymin(1) = 1
    iymax(1) = 1
    ! Block 2 (S)
    ixmin(2) = 1
    ixmax(2) = nlocal
    iymin(2) = nlocal
    iymax(2) = nlocal
    ! Block 3 (W)
    ixmin(3) = 1
    ixmax(3) = 1
    iymin(3) = 2
    iymax(3) = nlocal-1
    ! Block 4 (E)
    ixmin(4) = nlocal
    ixmax(4) = nlocal
    iymin(4) = 2
    iymax(4) = nlocal-1
    ! Block 5 (INTERIOR)
    if (overlap_comms) then
      ixmin(5) = 2
      ixmax(5) = nlocal-1
      iymin(5) = 2
      iymax(5) = nlocal-1
    else
      ! If there are no interior cells, do not overlap
      ! communications and calculations, just loop over interior cells
      ixmin(5) = 1
      ixmax(5) = nlocal
      iymin(5) = 1
      iymax(5) = nlocal
    end if
    dix = +1
    diy = +1
    color = 1
    ! When iteration backwards over the grid, reverse the direction
    if (direction == DIRECTION_BACKWARD) then
      do iblock = 1, 5
        tmp = ixmax(iblock)
        ixmax(iblock) = ixmin(iblock)
        ixmin(iblock) = tmp
        tmp = iymax(iblock)
        iymax(iblock) = iymin(iblock)
        iymin(iblock) = tmp
      end do
      dix = -1
      diy = -1
      color = 0
    end if
    nsweeps = 1
    if (ordering == ORDERING_LEX) then
      nsweeps = 1
    else if (ordering == ORDERING_RB) then
      nsweeps = 2
    end if
    do isweep = 1, nsweeps
      if (overlap_comms) then
        ! Loop over cells next to boundary (iblock = 1,...,4)
        do iblock = 1, 4
          call loop_over_grid(iblock)
        end do
        ! Initiate halo exchange
        call ihaloswap(level,m,u,send_requests,recv_requests)
      end if
      ! Loop over INTERIOR cells
      iblock = 5
      call loop_over_grid(iblock)
      if (overlap_comms) then
        if (m > 0) then
          call mpi_waitall(4,recv_requests, MNH_STATUSES_IGNORE, ierr)
        end if
      else
        call haloswap(level,m,u)
      end if
      color = 1-color
    end do

    ! Free memory again
    deallocate(r)
    deallocate(c)
    deallocate(utmp)

    contains

    !------------------------------------------------------------------
    ! Loop over grid, for a given block
    !------------------------------------------------------------------
    subroutine loop_over_grid(iblock)
      implicit none
      integer, intent(in) :: iblock
      integer :: ix,iy,iz
      do ix=ixmin(iblock),ixmax(iblock),dix
        do iy=iymin(iblock),iymax(iblock),diy
          if (ordering == ORDERING_RB) then
            if (mod((ix+u%ix_min)+(iy+u%iy_min),2) .ne. color) cycle
          end if
          call apply_tridiag_solve(ix,iy,r,c,b,         &
                                   u%s(1:nz,iy  ,ix+1), &
                                   u%s(1:nz,iy  ,ix-1), &
                                   u%s(1:nz,iy+1,ix  ), &
                                   u%s(1:nz,iy-1,ix  ), &
                                   utmp)
           ! Add to field with overrelaxation-factor
          do iz=1,nz
            u%s(iz,iy,ix) = (1.0_rl-rho)*u%s(iz,iy,ix) + rho*utmp(iz)
          end do
        end do
      end do
    end subroutine loop_over_grid

  end subroutine line_SOR

!==================================================================
! SSOR line smoother mnh
!==================================================================
  subroutine line_SSOR_mnh(level,m,direction,b,u)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    integer, intent(in)                :: direction
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    if (direction == DIRECTION_FORWARD) then
      call line_SOR_mnh(level,m,DIRECTION_FORWARD,b,u)
      call line_SOR_mnh(level,m,DIRECTION_BACKWARD,b,u)
    else
      call line_SOR_mnh(level,m,DIRECTION_BACKWARD,b,u)
      call line_SOR_mnh(level,m,DIRECTION_FORWARD,b,u)
    end if
  end subroutine line_SSOR_mnh

!==================================================================
! SSOR line smoother
!==================================================================
  subroutine line_SSOR(level,m,direction,b,u)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    integer, intent(in)                :: direction
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    if (direction == DIRECTION_FORWARD) then
      call line_SOR(level,m,DIRECTION_FORWARD,b,u)
      call line_SOR(level,m,DIRECTION_BACKWARD,b,u)
    else
      call line_SOR(level,m,DIRECTION_BACKWARD,b,u)
      call line_SOR(level,m,DIRECTION_FORWARD,b,u)
    end if
  end subroutine line_SSOR

!==================================================================
! Jacobi line smoother
!==================================================================
  subroutine line_Jacobi_mnh(level,m,b,u)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    integer :: ix,iy,iz, nz
    real(kind=rl), dimension(5) :: alpha_T
    integer :: nlocal, halo_size
    real(kind=rl) :: rho
    logical :: overlap_comms
    integer, dimension(4) :: send_requests, recv_requests
    integer, dimension(4) :: send_requestsT, recv_requestsT
    integer :: ixmin(5), ixmax(5)
    integer :: iymin(5), iymax(5)
    integer :: iblock, ierr

    integer :: iib,iie,ijb,ije
    
    real(kind=rl) , dimension(:,:,:) , pointer , contiguous :: zSut0_st , zu_st, zSr_st, zSutmp_st

    logical , save                      :: Ljacob_first_call = .true.
    logical , save , dimension(max_lev) :: Ljacob_first_call_level = .true.

    real(kind=rl), pointer , contiguous :: r(:)
    real(kind=rl), pointer , contiguous :: c(:), utmp(:)
    real(kind=rl), pointer , contiguous :: u0(:,:,:) 
    real(kind=rl), pointer , contiguous :: ut0(:,:,:)
    type(scalar3d) , pointer :: Sr,Sc,Sut0,Sutmp

    integer :: nib,nie,njb,nje,nkb,nke
    integer :: ii,ij,ik
    integer :: iindex
    integer :: n_lin

    integer :: nlocalx,nlocaly

    !
    !  init size , param
    !
    nz = u%grid_param%nz
    nlocal = u%ix_max -u%ix_min + 1
    nlocalx = u%ix_max -u%ix_min + 1
    nlocaly = u%iy_max -u%iy_min + 1
    halo_size = u%halo_size
    n_lin = (nlocalx+2*halo_size) * (nlocaly+2*halo_size) * (nz+2)

    nib = 1-halo_size
    nie = nlocalx+halo_size
    njb = 1-halo_size
    nje = nlocaly+halo_size
    nkb = 0
    nke = nz + 1
    
    ! Set optimal smoothing parameter on each level
    !rho = 2.0_rl/(2.0_rl+4.0_rl*model_param%omega2*u%grid_param%n**2/(1.0_rl+4.0_rl*model_param%omega2*u%grid_param%n**2))
    rho = smoother_param%rho
    
    ! Allocate data one for all by level
    if (Ljacob_first_call_level(level)) then
       Ljacob_first_call_level(level) = .false.
       if (LUseO) then
          allocate(Tjacobi(level)%u0(0:u%grid_param%nz+1,            &
               1-halo_size:nlocal+halo_size,   &
               1-halo_size:nlocal+halo_size) )
          !$acc enter data create(Tjacobi(level)%u0)
       end if
       if (LUseT) then
!!$          print*,"allocate(ut0),nib,nie,njb,nje,nkb,nke,level",nib,nie,njb,nje,nkb,nke, level
!!$          allocate(Tjacobi(level)%ut0(1-halo_size:nlocal+halo_size,   &
!!$               1-halo_size:nlocal+halo_size,   &
!!$               0:u%grid_param%nz+1) )
!!$          ut0 => Tjacobi(level)%ut0

!!$          allocate(ut0(nib:nie,njb:nje,nkb:nke) )
!!$          !$acc enter data create(ut0)          
          iindex = zt1d_discretisation_allocate3d(ut0,nib,nie,njb,nje,nkb,nke)
          
          Tjacobi(level)%ut0 => ut0
          
       end if
       ! Create residual vector       
       allocate(Tjacobi(level)%r(nz))
       r => Tjacobi(level)%r
       !$acc enter data create(r)
       ! Allocate memory for auxiliary vectors for Thomas algorithm      
       allocate(Tjacobi(level)%c(nz))
       c => Tjacobi(level)%c
       !$acc enter data create(c)
       allocate(Tjacobi(level)%utmp(nz))
       utmp => Tjacobi(level)%utmp
       !$acc enter data create(utmp)
       if (LUseT) then
          allocate(Tjacobi(level)%Sr)
!!$          allocate(Tjacobi(level)%Sr%st(1-halo_size:nlocal+halo_size,   &
!!$               1-halo_size:nlocal+halo_size,   &
!!$               0:u%grid_param%nz+1) )
!!$          zSr_st => Tjacobi(level)%Sr%st
!!$          !$acc enter data create(zSr_st)
          iindex = zt1d_discretisation_allocate3d(zSr_st,nib,nie,njb,nje,nkb,nke)
          Tjacobi(level)%Sr%st => zSr_st
          
          allocate(Tjacobi(level)%Sut0)
!!$          allocate(Tjacobi(level)%Sut0%st(1-halo_size:nlocal+halo_size,   &
!!$               1-halo_size:nlocal+halo_size,   &
!!$               0:u%grid_param%nz+1) )
!!$          zSut0_st => Tjacobi(level)%Sut0%st
!!$          !$acc enter data create(zSut0_st)
          iindex = zt1d_discretisation_allocate3d(zSut0_st,nib,nie,njb,nje,nkb,nke)
          Tjacobi(level)%Sut0%st => zSut0_st
          
          allocate(Tjacobi(level)%Sutmp)
!!$          allocate(Tjacobi(level)%Sutmp%st(1-halo_size:nlocal+halo_size,   &
!!$                            1-halo_size:nlocal+halo_size,   &
!!$                            0:u%grid_param%nz+1) )
!!$          zSutmp_st => Tjacobi(level)%Sutmp%st
!!$          !$acc enter data create(zSutmp_st)
          iindex = zt1d_discretisation_allocate3d(zSutmp_st,nib,nie,njb,nje,nkb,nke)
          Tjacobi(level)%Sutmp%st => zSutmp_st
       end if
    end if
    
#ifdef MG_DEBUG
    if (i_am_master_mpi)  write(STDOUT,*) ' ====== begin boundary_mnh_mnh ====== level=',level
    call print_scalaprod2(level , m  , u  , "Befor boundary_mnh u=" )
#endif
    call boundary_mnh(u)
#ifdef MG_DEBUG
    call print_scalaprod2(level , m  , u  , "After boundary_mnh u=" )
#endif
    
#ifdef OVERLAPCOMMS
    overlap_comms = (nlocal > 2)
#else
    overlap_comms = .false.
#endif

    ! Block 1 (N)
    ixmin(1) = 1
    ixmax(1) = nlocalx
    iymin(1) = 1
    iymax(1) = 1
    ! Block 2 (S)
    ixmin(2) = 1
    ixmax(2) = nlocalx
    iymin(2) = nlocaly
    iymax(2) = nlocaly
    ! Block 3 (W)
    ixmin(3) = 1
    ixmax(3) = 1
    iymin(3) = 2
    iymax(3) = nlocaly-1
    ! Block 4 (E)
    ixmin(4) = nlocalx
    ixmax(4) = nlocalx
    iymin(4) = 2
    iymax(4) = nlocaly-1
    ! Block 5 (INTERIOR)
    if (overlap_comms) then
      ixmin(5) = 2
      ixmax(5) = nlocalx-1
      iymin(5) = 2
      iymax(5) = nlocaly-1
    else
      ! If there are no interior cells, do not overlap
      ! communications and calculations, just loop over interior cells
      ixmin(5) = 1
      ixmax(5) = nlocalx
      iymin(5) = 1
      iymax(5) = nlocaly
    end if

    ! Temporary vector 
    if (LUseO) then
       u0 => Tjacobi(level)%u0
    end if
    if (LUseT) then
       ut0 => Tjacobi(level)%ut0 
    end if
    if (LUseO) u0(:,:,:) = u%s(:,:,:)
    if (LUseT) then
       zu_st =>  u%st
!!$       !$acc kernels present_cr(ut0,zu_st)
!!$       ! mnh_do_concurrent ( ii=nib:nie , ij=njb:nje , ik=nkb:nke )
!!$          ut0(ii,ij,ik) = zu_st(ii,ij,ik)
!!$       ! mnh_end_do()
!!$       !$acc end kernels
       call dcopy(n_lin,zu_st,1,ut0,1)   
    end if
    ! Create residual vector
    r => Tjacobi(level)%r
    ! Allocate memory for auxiliary vectors for Thomas algorithm
    c => Tjacobi(level)%c    
    utmp => Tjacobi(level)%utmp
    if (LUseT) then
       Sr => Tjacobi(level)%Sr
       Sut0 => Tjacobi(level)%Sut0
       
       zSut0_st => Sut0%st
       zu_st => u%st

!!$       !$acc kernels present(zSut0_st,zu_st)
!!$       ! mnh_do_concurrent( ii=nib:nie , ij=njb:nje , ik=nkb:nke )
!!$           zSut0_st(ii,ij,ik) = zu_st(ii,ij,ik)
!!$       ! mnh_end_do()
!!$       !$acc end kernels
       call dcopy(n_lin,zu_st,1,zSut0_st,1)    

       Sutmp => Tjacobi(level)%Sutmp
       
    endif

    ! Loop over grid
    if (overlap_comms) then
    ! Loop over cells next to boundary (iblock = 1,...,4)
      do iblock = 1, 4
        call loop_over_grid_jacobi_mnh(iblock)
      end do
      ! Initiate halo exchange
      call ihaloswap_mnh(level,m,u,send_requests,recv_requests,send_requestsT,recv_requestsT)
    end if
    ! Loop over INTERIOR cells
    iblock = 5
    call loop_over_grid_jacobi_mnh(iblock)
    if (overlap_comms) then
      if (m > 0) then
        if (LUseO) call mpi_waitall(4,recv_requests, MNH_STATUSES_IGNORE, ierr)
        if (LUseO) call mpi_waitall(4,send_requests, MNH_STATUSES_IGNORE, ierr)
        if (LUseT) call mpi_waitall(4,recv_requestsT, MNH_STATUSES_IGNORE, ierr)
        if (LUseT) call mpi_waitall(4,send_requestsT, MNH_STATUSES_IGNORE, ierr)
      end if
    else
      call haloswap_mnh(level,m,u)
    end if

    ! Free memory again
!!$    deallocate(r)
!!$    deallocate(c)
!!$    if (LUseO) deallocate(u0)
!!$    if (LUseT) deallocate(ut0)
!!$    deallocate(utmp)
!!$    if (LUseT) deallocate(Sr%st,Sut0%st,Sutmp%st)

  contains

  subroutine loop_over_grid_jacobi_mnh(iblock)
    implicit none
    integer, intent(in) :: iblock
    integer :: ix,iy,iz

    real(kind=rl) , dimension(:,:,:) , pointer , contiguous :: zu_st , zSutmp_st , zSut0_st
    
    if (LUseO) then
       do ix=ixmin(iblock),ixmax(iblock)
          do iy=iymin(iblock),iymax(iblock)
             call apply_tridiag_solve_mnh(ix,iy,r,c,b,        &
                  u0(1:nz,iy  ,ix+1), &
                  u0(1:nz,iy  ,ix-1), &
                  u0(1:nz,iy+1,ix  ), &
                  u0(1:nz,iy-1,ix  ), &
                  utmp)
             ! Add correction
             do iz=1,nz
                u%s(iz,iy,ix) = rho*utmp(iz) + (1.0_rl-rho)*u0(iz,iy,ix)
             end do
          end do
       end do
    end if

    if (LUseT) then
       
       iib=ixmin(iblock)
       iie=ixmax(iblock)
       ijb=iymin(iblock)
       ije=iymax(iblock)
       
       zu_st => u%st
       zSutmp_st => Sutmp%st
       zSut0_st => Sut0%st
       
       call apply_tridiag_solve_mnh_allT(iib,iie,ijb,ije,Sr,c,b,        &
            Sut0, &
            Sutmp,level )
       
       call loop_over_grid_jacobi_mnh_dim(&
            zu_st,zSutmp_st,zSut0_st,&
            iib,iie,ijb,ije)
       
    end if

  end subroutine loop_over_grid_jacobi_mnh

! contains

  subroutine loop_over_grid_jacobi_mnh_dim(&
            zu_st,zSutmp_st,zSut0_st,&
            iib,iie,ijb,ije)
    
    implicit none

    real(kind=rl) :: zu_st    (nib:nie,njb:nje,nkb:nke ),&
            zSutmp_st(nib:nie,njb:nje,nkb:nke ),&
            zSut0_st (nib:nie,njb:nje,nkb:nke )

    integer :: iib,iie,ijb,ije

    ! local var
    integer :: ix,iy,iz
    
    !$acc kernels present_cr(zsut0_st,zSutmp_st,zu_st)
    !$mnh_do_concurrent( ix=iib:iie  , iy=ijb:ije , iz=1:nz )
    zu_st(ix,iy,iz) = & 
         rho*zSutmp_st(ix,iy,iz) & 
         + (1.0_rl-rho)*zSut0_st(ix,iy,iz)
    !$mnh_end_do() ! concurrent
    !$acc end kernels
    
  end subroutine loop_over_grid_jacobi_mnh_dim
    
end subroutine line_Jacobi_mnh
!==================================================================
! Jacobi line smoother
!==================================================================
  subroutine line_Jacobi(level,m,b,u)
    implicit none
    integer, intent(in)                :: level
    integer, intent(in)                :: m
    type(scalar3d), intent(in)         :: b
    type(scalar3d), intent(inout)      :: u
    real(kind=rl), allocatable :: r(:)
    integer :: ix,iy,iz, nz
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl), allocatable :: c(:), utmp(:)
    real(kind=rl), allocatable :: u0(:,:,:)
    integer :: nlocal, halo_size
    real(kind=rl) :: rho
    logical :: overlap_comms
    integer, dimension(4) :: send_requests, recv_requests
    integer :: ixmin(5), ixmax(5)
    integer :: iymin(5), iymax(5)
    integer :: iblock, ierr

    ! Set optimal smoothing parameter on each level
    rho = 2.0_rl/(2.0_rl+4.0_rl*model_param%omega2*u%grid_param%n**2/(1.0_rl+4.0_rl*model_param%omega2*u%grid_param%n**2))

    nz = u%grid_param%nz
    nlocal = u%ix_max -u%ix_min + 1
    halo_size = u%halo_size

#ifdef OVERLAPCOMMS
    overlap_comms = (nlocal > 2)
#else
    overlap_comms = .false.
#endif

    ! Block 1 (N)
    ixmin(1) = 1
    ixmax(1) = nlocal
    iymin(1) = 1
    iymax(1) = 1
    ! Block 2 (S)
    ixmin(2) = 1
    ixmax(2) = nlocal
    iymin(2) = nlocal
    iymax(2) = nlocal
    ! Block 3 (W)
    ixmin(3) = 1
    ixmax(3) = 1
    iymin(3) = 2
    iymax(3) = nlocal-1
    ! Block 4 (E)
    ixmin(4) = nlocal
    ixmax(4) = nlocal
    iymin(4) = 2
    iymax(4) = nlocal-1
    ! Block 5 (INTERIOR)
    if (overlap_comms) then
      ixmin(5) = 2
      ixmax(5) = nlocal-1
      iymin(5) = 2
      iymax(5) = nlocal-1
    else
      ! If there are no interior cells, do not overlap
      ! communications and calculations, just loop over interior cells
      ixmin(5) = 1
      ixmax(5) = nlocal
      iymin(5) = 1
      iymax(5) = nlocal
    end if

    ! Temporary vector
    allocate(u0(0:u%grid_param%nz+1,            &
                1-halo_size:nlocal+halo_size,   &
                1-halo_size:nlocal+halo_size) )
    u0(:,:,:) = u%s(:,:,:)
    ! Create residual vector
    allocate(r(nz))
    ! Allocate memory for auxiliary vectors for Thomas algorithm
    allocate(c(nz))
    allocate(utmp(nz))

    ! Loop over grid
    if (overlap_comms) then
    ! Loop over cells next to boundary (iblock = 1,...,4)
      do iblock = 1, 4
        call loop_over_grid(iblock)
      end do
      ! Initiate halo exchange
      call ihaloswap(level,m,u,send_requests,recv_requests)
    end if
    ! Loop over INTERIOR cells
    iblock = 5
    call loop_over_grid(iblock)
    if (overlap_comms) then
      if (m > 0) then
        call mpi_waitall(4,recv_requests, MNH_STATUSES_IGNORE, ierr)
      end if
    else
      call haloswap(level,m,u)
    end if

    ! Free memory again
    deallocate(r)
    deallocate(c)
    deallocate(u0)
    deallocate(utmp)

  contains

  subroutine loop_over_grid(iblock)
    implicit none
    integer, intent(in) :: iblock
    integer :: ix,iy,iz
    do ix=ixmin(iblock),ixmax(iblock)
      do iy=iymin(iblock),iymax(iblock)
        call apply_tridiag_solve(ix,iy,r,c,b,        &
                                 u0(1:nz,iy  ,ix+1), &
                                 u0(1:nz,iy  ,ix-1), &
                                 u0(1:nz,iy+1,ix  ), &
                                 u0(1:nz,iy-1,ix  ), &
                                 utmp)
        ! Add correction
        do iz=1,nz
          u%s(iz,iy,ix) = rho*utmp(iz) + (1.0_rl-rho)*u0(iz,iy,ix)
        end do
      end do
    end do
  end subroutine loop_over_grid

  end subroutine line_Jacobi
!==================================================================
! At a given horizontal position (ix,iy) (local coordinates),
! calculate
!
! u_out = T(ix,iy)^{-1} (b_(ix,iy)
!       - sum_{ix',iy' != ix,iy} A_{(ix,iy),(ix',iy')}*u_in(ix',iy'))
!
!==================================================================
  subroutine apply_tridiag_solve_mnh(ix,iy,r,c,b, &
                                 u_in_1,      &
                                 u_in_2,      &
                                 u_in_3,      &
                                 u_in_4,      &
                                 u_out)

    implicit none
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    real(kind=rl), intent(inout), dimension(:) :: r
    real(kind=rl), intent(inout), dimension(:) :: c
    type(scalar3d), intent(in) :: b
    real(kind=rl), intent(in), dimension(:) :: u_in_1
    real(kind=rl), intent(in), dimension(:) :: u_in_2
    real(kind=rl), intent(in), dimension(:) :: u_in_3
    real(kind=rl), intent(in), dimension(:) :: u_in_4
    real(kind=rl), intent(inout), dimension(:) :: u_out
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) :: Tij
    real(kind=rl) :: alpha_div_Tij, tmp, b_k_tmp, c_k_tmp
    integer :: iz, nz

    real(kind=rl) :: xctop_boot 

    nz = b%grid_param%nz
    xctop_boot = 0.0

    call construct_alpha_T_mnh(b%grid_param,  &
                           ix+b%ix_min-1, &
                           iy+b%iy_min-1, &
                           alpha_T,Tij)
    ! Calculate r_i = b_i - A_{ij} u_i
    !alpha_T(5) = 4
    if (LUseO) then 
       iz=1 
       r(iz) = b%s(iz,iy,ix)
       do iz=2,nz-1
          r(iz) = b%s(iz,iy,ix) - vert_coeff%d(iz) * ( &
               alpha_T(1) * u_in_1(iz) + &
               alpha_T(2) * u_in_2(iz) + &
               alpha_T(3) * u_in_3(iz) + &
               alpha_T(4) * u_in_4(iz) )
       end do
       iz=nz
       r(iz) = b%s(iz,iy,ix)
       
       ! Thomas algorithm
       ! STEP 1: Create modified coefficients
       iz = 1
       alpha_div_Tij = alpha_T(5)/Tij
       tmp = (vert_coeff%a(iz)-vert_coeff%b(iz)-vert_coeff%c(iz)) &
            - xctop_boot*alpha_div_Tij
       c(iz) = vert_coeff%b(iz)/tmp
       u_out(iz) = r(iz) / (tmp*Tij*vert_coeff%d(iz))
       do iz=2,nz-1
          b_k_tmp = vert_coeff%b(iz)
          c_k_tmp = vert_coeff%c(iz)
          tmp = ((vert_coeff%a(iz)-b_k_tmp-c_k_tmp)-alpha_div_Tij) &
               - c(iz-1)*c_k_tmp
          c(iz) = b_k_tmp / tmp
          u_out(iz) = (r(iz) / (Tij*vert_coeff%d(iz)) - u_out(iz-1)*c_k_tmp) / tmp
       end do
       iz=nz
       b_k_tmp = vert_coeff%b(iz)
       c_k_tmp = vert_coeff%c(iz)
       tmp = ((vert_coeff%a(iz)-b_k_tmp-c_k_tmp)- xctop_boot*alpha_div_Tij) &
            - c(iz-1)*c_k_tmp
       c(iz) = b_k_tmp / tmp
       u_out(iz) = (r(iz) / (Tij*vert_coeff%d(iz)) - u_out(iz-1)*c_k_tmp) / tmp
       
       ! STEP 2: back substitution
       do iz=nz-1,1,-1
          u_out(iz) = u_out(iz) - c(iz) * u_out(iz+1)
       end do
    end if
    !

  end subroutine apply_tridiag_solve_mnh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tranpose version 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine apply_tridiag_solve_mnhT(ix,iy,r,c,b, &
                                 u_in_1,      &
                                 u_in_2,      &
                                 u_in_3,      &
                                 u_in_4,      &
                                 u_out)

    implicit none
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    real(kind=rl), intent(inout), dimension(:) :: r
    real(kind=rl), intent(inout), dimension(:) :: c
    type(scalar3d), intent(in) :: b
    real(kind=rl), intent(in), dimension(:) :: u_in_1
    real(kind=rl), intent(in), dimension(:) :: u_in_2
    real(kind=rl), intent(in), dimension(:) :: u_in_3
    real(kind=rl), intent(in), dimension(:) :: u_in_4
    real(kind=rl), intent(inout), dimension(:) :: u_out
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) :: Tij
    real(kind=rl) :: alpha_div_Tij, tmp, b_k_tmp, c_k_tmp
    integer :: iz, nz

    real(kind=rl) :: xctop_boot 

    nz = b%grid_param%nz
    xctop_boot = 0.0

    call construct_alpha_T_cst_mnh(b%grid_param,alpha_T,Tij)
    ! Calculate r_i = b_i - A_{ij} u_i
    if (LUseT ) then 
       iz=1 
       r(iz) = b%st(ix,iy,iz)
       do iz=2,nz-1
          r(iz) = b%st(ix,iy,iz) - vert_coeff%d(iz) * ( &
               alpha_T(1) * u_in_1(iz) + &
               alpha_T(2) * u_in_2(iz) + &
               alpha_T(3) * u_in_3(iz) + &
               alpha_T(4) * u_in_4(iz) )
       end do
       iz=nz
       r(iz) = b%st(ix,iy,iz)
       
       ! Thomas algorithm
       ! STEP 1: Create modified coefficients
       iz = 1
       alpha_div_Tij = alpha_T(5)/Tij
       tmp = (vert_coeff%a(iz)-vert_coeff%b(iz)-vert_coeff%c(iz)) &
            - xctop_boot*alpha_div_Tij
       c(iz) = vert_coeff%b(iz)/tmp
       u_out(iz) = r(iz) / (tmp*Tij*vert_coeff%d(iz))
       do iz=2,nz-1
          b_k_tmp = vert_coeff%b(iz)
          c_k_tmp = vert_coeff%c(iz)
          tmp = ((vert_coeff%a(iz)-b_k_tmp-c_k_tmp)-alpha_div_Tij) &
               - c(iz-1)*c_k_tmp
          c(iz) = b_k_tmp / tmp
          u_out(iz) = (r(iz) / (Tij*vert_coeff%d(iz)) - u_out(iz-1)*c_k_tmp) / tmp
       end do
       iz=nz
       b_k_tmp = vert_coeff%b(iz)
       c_k_tmp = vert_coeff%c(iz)
       tmp = ((vert_coeff%a(iz)-b_k_tmp-c_k_tmp)- xctop_boot*alpha_div_Tij) &
            - c(iz-1)*c_k_tmp
       c(iz) = b_k_tmp / tmp
       u_out(iz) = (r(iz) / (Tij*vert_coeff%d(iz)) - u_out(iz-1)*c_k_tmp) / tmp
       
       ! STEP 2: back substitution
       do iz=nz-1,1,-1
          u_out(iz) = u_out(iz) - c(iz) * u_out(iz+1)
       end do
    end if

  end subroutine apply_tridiag_solve_mnhT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tranpose version all xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine apply_tridiag_solve_mnh_allT(iib,iie,ijb,ije, &
                                 Sr,c,b,Su_in,Su_out,level)

    implicit none
    integer, intent(in) :: iib,iie,ijb,ije
    type(scalar3d), intent(inout) :: Sr
    real(kind=rl), intent(inout), dimension(:) :: c
    type(scalar3d), intent(in)    :: b
    type(scalar3d), intent(in)    :: Su_in
    type(scalar3d), intent(inout) :: Su_out
    integer, intent(in) :: level

    !local 
    !real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) ,pointer :: Tij
    real(kind=rl) ,pointer :: alpha_div_Tij ! b_k_tmp, c_k_tmp
    integer :: iz, nz

    real(kind=rl), dimension(:,:,:) , pointer, contiguous :: zSr_st , zb_st , zSu_in_st , zSu_out_st
    real(kind=rl), dimension(:)     , pointer, contiguous :: za_k, zb_k, zc_k, zd_k 
!!$    integer :: ii,ij

    integer , parameter :: max_lev = 128
    logical , save                                 :: Lfirst_call_tridiag_mnhallT = .true.
    logical , save , dimension(max_lev)            :: Lfirst_call_level_tridiag_mnhallT = .true.

    real(kind=rl), dimension(:), pointer, contiguous :: tmp_k,c_k
    
    type Temp_tridiag_mnh
       real(kind=rl), dimension(:), pointer, contiguous :: tmp_k,c_k
       real(kind=rl), pointer  :: Tij , alpha_div_Tij
    end type Temp_tridiag_mnh

    type(Temp_tridiag_mnh) , save , dimension(max_lev) :: Ttridiag_mnh

    real(kind=rl) :: Tijs

    integer :: nib,nie,njb,nje,nzb,nze
    
    if (LUseT ) then

      ! Calculate r_i = b_i - A_{ij} u_i
       
      nz = b%grid_param%nz
           
      zSr_st => Sr%st
      zb_st  => b%st
      zSu_in_st => Su_in%st
      zSu_out_st => Su_out%st
      za_k  => vert_coeff%a
      zb_k  => vert_coeff%b
      zc_k  => vert_coeff%c
      zd_k  => vert_coeff%d

      nib = Lbound(zb_st,1) ; nie = Ubound(zb_st,1)
      njb = Lbound(zb_st,2) ; nje = Ubound(zb_st,2)
      nzb = Lbound(zb_st,3) ; nze = Ubound(zb_st,3)

      if ( Lfirst_call_level_tridiag_mnhallT(level) ) then
         Lfirst_call_level_tridiag_mnhallT(level) = .false.

         allocate(Ttridiag_mnh(level)%Tij)
         allocate(Ttridiag_mnh(level)%alpha_div_Tij)

         Tij=>Ttridiag_mnh(level)%Tij
         alpha_div_Tij=>Ttridiag_mnh(level)%alpha_div_Tij
         !call construct_alpha_T_cst_mnh(b%grid_param,alpha_T,Tij)
         Tij = ( b%grid_param%L/b%grid_param%n ) ** 2
         alpha_div_Tij = 4.0_rl / Tij
         !$acc enter data copyin(Tij,alpha_div_Tij)
         !print*,"level=",level," Tij=",Tij," alpha_div_Tij=",alpha_div_Tij

         allocate(Ttridiag_mnh(level)%tmp_k(size(zb_k)))
         allocate(Ttridiag_mnh(level)%c_k(size(zb_k)))
         
         tmp_k => Ttridiag_mnh(level)%tmp_k
         c_k => Ttridiag_mnh(level)%c_k
         ! Thomas algorithm
         ! STEP 1: Create modified coefficients
         iz = 1
         tmp_k(iz) = (za_k(iz)-zb_k(iz)-zc_k(iz))
         c_k(iz) = zb_k(iz)/tmp_k(iz)
         !
         do iz=2,nz-1
            tmp_k(iz) = ((za_k(iz)-zb_k(iz)-zc_k(iz))-alpha_div_Tij) &
                 - c_k(iz-1)*zc_k(iz)
            c_k(iz) = zb_k(iz) / tmp_k(iz)
         end do
         !
         iz=nz
         tmp_k(iz) = ((za_k(iz)-zb_k(iz)-zc_k(iz))) &
              - c_k(iz-1)*zc_k(iz)
         c_k(iz) = zb_k(iz) / tmp_k(iz)

         !$acc enter data copyin(tmp_k,c_k)
         
      endif

      Tij=>Ttridiag_mnh(level)%Tij
      Tijs = Tij ! Bypass Cray Bug with allocatable/scalair pointer
      
      tmp_k => Ttridiag_mnh(level)%tmp_k
      c_k => Ttridiag_mnh(level)%c_k

!!$      print*,"apply_tridiag::level,Sr%ix_min/%ix_max/%icompx_min/%icompx_max",&
!!$           Sr%ix_min,Sr%ix_max,Sr%icompx_min,Sr%icompx_max
!!$      print*,"apply_tridiag::level,iib,iie,ijb,ije,nz",level,iib,iie,ijb,ije,nz
!!$      print*,"apply_tridiag::level,nib,nie,njb,nje,nzb,nze",level,&
!!$                                   nib,nie,njb,nje,nzb,nze
!!$      print*,"apply_tridiag::level, zSr_st=L/U/S",level,Lbound(zSr_st),Ubound(zSr_st),Shape(zSr_st)
!!$      print*,"apply_tridiag::lvel, zc_k  =L/U/S",level,Lbound(zc_k),Ubound(zc_k),Shape(zc_k)

      call  apply_tridiag_solve_mnh_allT_dim(&
           zSr_st,zSu_out_st,zb_st,zSu_in_st,&
           zc_k,zd_k,tmp_k,c_k &
      )
       
      end if

    contains
      
      subroutine  apply_tridiag_solve_mnh_allT_dim(&
           pzSr_st,pzSu_out_st,pzb_st,pzSu_in_st,&
           pzc_k,pzd_k,ptmp_k,pc_k &
        )

        implicit none

        real(kind=rl) :: pzSr_st    (nib:nie,njb:nje,nzb:nze) ,&
                pzSu_out_st(nib:nie,njb:nje,nzb:nze) ,&
                pzb_st     (nib:nie,njb:nje,nzb:nze) ,&
                pzSu_in_st (nib:nie,njb:nje,nzb:nze) ,&
                pzc_k(nz),pzd_k(nz),ptmp_k(nz),pc_k(nz)
        !
        ! local var
        !
        integer :: ii,ij,iz
                   
        ! acc kernels present(pzSr_st,pzSu_out_st,pzb_st,pzSu_in_st)
        !$acc parallel present(pzSr_st,pzSu_out_st,pzb_st,pzSu_in_st) &
        !$acc &        present(pzc_k,pzd_k,ptmp_k,pc_k)
        !$mnh_do_concurrent ( ii=iib:iie , ij=ijb:ije )   
        iz=1 
        pzSr_st(ii,ij,iz) = pzb_st(ii,ij,iz)
        !$acc loop seq
        do iz=2,nz-1
           pzSr_st(ii,ij,iz) = pzb_st(ii,ij,iz) - pzd_k(iz) * ( &
                pzSu_in_st(ii+1,ij,iz) + &
                pzSu_in_st(ii-1,ij,iz) + &
                pzSu_in_st(ii,ij+1,iz) + &
                pzSu_in_st(ii,ij-1,iz) )
        end do
        iz=nz
        pzSr_st(ii,ij,iz) = pzb_st(ii,ij,iz)
        !
        ! Thomas algorithm
        !
        iz = 1     
        pzSu_out_st(ii,ij,iz) = pzSr_st(ii,ij,iz) / (ptmp_k(iz)*Tijs*pzd_k(iz))
        !$acc loop seq
        do iz=2,nz-1
           pzSu_out_st(ii,ij,iz) = (pzSr_st(ii,ij,iz) / (Tijs*pzd_k(iz)) & 
                - pzSu_out_st(ii,ij,iz-1)*pzc_k(iz)) / ptmp_k(iz)
        end do
        iz=nz
        pzSu_out_st(ii,ij,iz) = (pzSr_st(ii,ij,iz) / (Tijs*pzd_k(iz)) &
             -  pzSu_out_st(ii,ij,iz-1)*pzc_k(iz)) / ptmp_k(iz)         
        !     
        ! STEP 2: back substitution
        !$acc loop seq
        do iz=nz-1,1,-1             
           pzSu_out_st(ii,ij,iz) = pzSu_out_st(ii,ij,iz) & 
                - pc_k(iz) * pzSu_out_st(ii,ij,iz+1)
        end do
        !$mnh_end_do()   
        !$acc end parallel
        ! acc end kernels
         
      end subroutine apply_tridiag_solve_mnh_allT_dim

  end subroutine apply_tridiag_solve_mnh_allT
  !==================================================================
! At a given horizontal position (ix,iy) (local coordinates),
! calculate
!
! u_out = T(ix,iy)^{-1} (b_(ix,iy)
!       - sum_{ix',iy' != ix,iy} A_{(ix,iy),(ix',iy')}*u_in(ix',iy'))
!
!==================================================================
  subroutine apply_tridiag_solve(ix,iy,r,c,b, &
                                 u_in_1,      &
                                 u_in_2,      &
                                 u_in_3,      &
                                 u_in_4,      &
                                 u_out)

    implicit none
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    real(kind=rl), intent(inout), dimension(:) :: r
    real(kind=rl), intent(inout), dimension(:) :: c
    type(scalar3d), intent(in) :: b
    real(kind=rl), intent(in), dimension(:) :: u_in_1
    real(kind=rl), intent(in), dimension(:) :: u_in_2
    real(kind=rl), intent(in), dimension(:) :: u_in_3
    real(kind=rl), intent(in), dimension(:) :: u_in_4
    real(kind=rl), intent(inout), dimension(:) :: u_out
    real(kind=rl), dimension(5) :: alpha_T
    real(kind=rl) :: Tij
    real(kind=rl) :: alpha_div_Tij, tmp, b_k_tmp, c_k_tmp
    integer :: iz, nz

    nz = b%grid_param%nz

    call construct_alpha_T(b%grid_param,  &
                           ix+b%ix_min-1, &
                           iy+b%iy_min-1, &
                           alpha_T,Tij)
    ! Calculate r_i = b_i - A_{ij} u_i
    !alpha_T(5) = 4.0
    do iz=1,nz
      r(iz) = b%s(iz,iy,ix) - vert_coeff%d(iz) * ( &
                alpha_T(1) * u_in_1(iz) + &
                alpha_T(2) * u_in_2(iz) + &
                alpha_T(3) * u_in_3(iz) + &
                alpha_T(4) * u_in_4(iz) )
    end do
    !r(1:nz) = b%s(1:nz,iy,ix)
    ! Thomas algorithm
    ! STEP 1: Create modified coefficients
    iz = 1
    alpha_div_Tij = alpha_T(5)/Tij
    tmp = (vert_coeff%a(iz)-vert_coeff%b(iz)-vert_coeff%c(iz)) &
             - alpha_div_Tij
    c(iz) = vert_coeff%b(iz)/tmp
    u_out(iz) = r(iz) / (tmp*Tij*vert_coeff%d(iz))
    do iz=2,nz
      b_k_tmp = vert_coeff%b(iz)
      c_k_tmp = vert_coeff%c(iz)
      tmp = ((vert_coeff%a(iz)-b_k_tmp-c_k_tmp)-alpha_div_Tij) &
          - c(iz-1)*c_k_tmp
      c(iz) = b_k_tmp / tmp
      u_out(iz) = (r(iz) / (Tij*vert_coeff%d(iz)) - u_out(iz-1)*c_k_tmp) / tmp
    end do
    ! STEP 2: back substitution
    do iz=nz-1,1,-1
      u_out(iz) = u_out(iz) - c(iz) * u_out(iz+1)
    end do
  end subroutine apply_tridiag_solve

end module discretisation
