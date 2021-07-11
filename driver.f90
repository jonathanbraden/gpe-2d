#include "macros.h"
program gpe_2d
  use constants, only : dl, twopi
  use utils, only : newunit
  use Hamiltonian
  implicit none

  type Time_Params
     real(dl) :: dt, dt_out
     integer :: nstep, nout
  end type Time_Params

  type Spec_Params
     real(dl) :: dk, m2
     integer :: ncut
     integer :: type
  end type Spec_Params

  type(Lattice) :: mySim
  type(Time_Params) :: tPar
  
  integer :: u1, u2
  
  call create_lattice(mySim,4,1._dl,2)
  
  call initialize_rand(42,1315)
  open(unit=newunit(u1),file='field.dat',access='stream',status='replace')
  open(unit=newunit(u2),file='dfield.dat',access='stream',status='replace')

  mySim%psi(:,:,1) = 1./2.**0.5
  mySim%psi(:,:,2) = 1./2.**0.5*iImag
  
contains

  subroutine sample_ics(nSamp)
    integer, intent(in) :: nSamp
    integer :: i

    do i=1,nSamp

    enddo
  end subroutine sample_ics
  
  subroutine run_simulation(sim,tPar)
    type(Lattice), intent(inout) :: sim
    type(Time_Params), intent(in) :: tPar
    integer :: i

    do i=1,tPar%nstep
       !call step_lattice(sim,tPar%dt, )
       !call output_lattice_mean(sim)
    enddo
  end subroutine run_simulation
  
  subroutine initialize_rand(seed, seedFac)
    integer, intent(in) :: seed, seedFac
    integer :: nseed, i
    integer, allocatable, dimension(:) :: seeds

    call random_seed(size=nseed)
    print*,"Seed size is ",nseed
    allocate(seeds(1:nseed))
    seeds = seed + seedfac*(/ (i-1, i=1,nseed) /)
    call random_seed(put=seeds)
    deallocate(seeds)
  end subroutine initialize_rand

  subroutine add_white_noise(this)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nlat) :: amp, phase
    integer :: nx, nnx, i, j, jj, l
    real(dl) :: norm
    
    !! ADD: Boot random number generator if needed
    nx = this%nLat; nnx = nx/2+1
    norm = 1._dl / this%lSize

    do l=1,this%nFld
       do j=1,nx
          call random_number(amp(1:nx)); call random_number(phase(1:nx))
          this%tPair%specSpace(:,j) = norm*sqrt(-log(amp))*exp(iImag*twopi*phase)
       enddo
       this%tPair%specSpace(1,1) = 0._dl
       call fftw_execute_dft(this%tPair%planb,this%tPair%specSpace,this%tPair%realSpace)
       this%psi(:,:,l) = this%psi(:,:,l) + this%tPair%realSpace
    enddo
  end subroutine add_white_noise

  
end program gpe_2d
