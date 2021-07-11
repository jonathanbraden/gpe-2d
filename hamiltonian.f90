#include "macros.h"
!
! Currently solves decoupled nonlinear Schrodinger equations
!
module Hamiltonian
  use constants, only : dl, twopi
  use fftw3_cmplx
  use simulation
  implicit none

  integer, parameter :: n_Hamiltonian_terms = 2
  integer, parameter :: nPar = 1
  real(dl) :: mu1, mu2
  real(dl) :: lam = 1.2_dl
  real(dl) :: mu = 1._dl
  
contains

  subroutine set_model_params(par)
    real(dl), dimension(1:nPar), intent(in) :: par
    lam = par(1)
    mu = 0._dl  ! Compute mu using the RMS of the fields
  end subroutine set_model_params

  !>@brief
  !> Compute chemical potential using variance of underlying fields
  subroutine compute_chem_pot(this)
    type(Lattice), intent(in) :: this

    mu = sum(abs(this%psi)**2)/dble(size(this%psi))
    mu1 = sum(abs(this%psi(:,:,1)**2))/dble(this%nLat)**2
    mu2 = sum(abs(this%psi(:,:,2)**2))/dble(this%nLat)**2
  end subroutine compute_chem_pot
  
  subroutine write_model_header(fNum)
    integer, intent(in) :: fNum
    logical :: o

    inquire(opened=o,unit=fNum)
    if (o) then
       write(fNum,*) "# Model is : "
    endif
  end subroutine write_model_header

  subroutine Hamiltonian_Split(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call Hamiltonian_kinetic(this,dt)
    case (2)
       call Hamiltonian_nonlinear(this,dt)
    case default
       print*,"Undefined Hamiltonian term"
       stop
    end select
  end subroutine Hamiltonian_Split

  subroutine symp_o2_step(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_Hamiltonian_terms-1
       call Hamiltonian_Split(this,0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(this,w1*dt,n_Hamiltonian_terms)
    do i=n_Hamiltonian_terms-1,2,-1
       call Hamiltonian_Split(this,0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(this,0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2_step

  subroutine Hamiltonian_kinetic(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: n,nn, l, i,j,jj
    real(dl) :: norm, dk, rad2
    integer, dimension(1:this%nLat) :: n_x
    
    ! In here, it will probably be more convenient to replace the inner loop with precomputed kx2 values
    n = this%nLat; nn = n/2+1
    norm = 1._dl/dble(n)**2
    dk = this%dk

    do i=1,nn; n_x(i) = i-1; enddo; do i=nn+1,n; n_x(i) = n+1-i; enddo ! Precomputed for now.  Try using this to vectorize the inner loop below

#ifdef SPECTRAL
    do l=1,this%nFld
       this%tPair%realSpace(1:n,1:n) = this%psi(1:n,1:n,l)
       call forward_transform_2d_wtype_c(this%tPair)
       do j=1,n; if (j<=nn) then; jj = j-1; else; jj = n+1-j; endif
          do i=1,nn
             rad2 = dble(jj**2+(i-1)**2)*dk**2
             this%tPair%specSpace = norm*exp(-0.5_dl*iImag*rad2*dt)*this%tPair%specSpace
          enddo
          do i=nn+1,n
             rad2 = dble(jj**2+(n+1-i)**2)*dk**2
             this%tPair%specSpace = norm*exp(-0.5_dl*iImag*rad2*dt)*this%tPair%specSpace
          enddo
       enddo
       call backward_transform_2d_wtype_c(this%tPair)
       this%psi(1:n,1:n,l) = this%tPair%realSpace(1:n,1:n)
    enddo
#else
    call wrap_field(this)
    ! Write this
#endif
  end subroutine Hamiltonian_kinetic

  !>@brief
  !>
  !> Solves df = -i|f|^2f
  subroutine Hamiltonian_nonlinear(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: j

    do j=1,this%nFld
       ! This is probably the place to include the chem pot associated with |psi|^2
       this%psi(:,:,j) = exp(-iImag*(abs(this%psi(:,:,j))**2-mu)*dt)*this%psi(:,:,j)
    enddo
  end subroutine Hamiltonian_nonlinear
  
end module Hamiltonian
