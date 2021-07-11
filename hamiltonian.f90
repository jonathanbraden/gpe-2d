#include "macros.h"
module Hamiltonian
  use constants, only : dl, twopi
  use fftw3_cmplx
  use simulation
  implicit none

  integer, parameter :: n_Hamiltonian_terms = 2
  integer, parameter :: nPar = 1
  real(dl) :: lam = 1.2_dl
  real(dl) :: mu = 1._dl
  
contains

  subroutine set_model_params(par)
    real(dl), dimension(1:nPar), intent(in) :: par
    lam = par(1)
    mu = 0._dl
  end subroutine set_model_params

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
       call Hamiltonian_linear(this,dt)
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

  subroutine Hamiltonian_linear(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: n, j

    n = this%nLat
#ifdef SPECTRAL
    do j=1,this%nFld
       this%tPair%realSpace(1:n,1:n) = this%psi(1:n,1:n,j)
       call laplacian_2d_wtype_c(this%tPair,this%dk)
       this%psi(1:n,1:n,j) = this%psi(1:n,1:n,j) + dt*iImag*this%tPair%realSpace(1:n,1:n)
    enddo
#else
    call wrap_field(this)
    ! Write this
#endif
  end subroutine Hamiltonian_linear

  subroutine Hamiltonian_nonlinear(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    
  end subroutine Hamiltonian_nonlinear
  
end module Hamiltonian
