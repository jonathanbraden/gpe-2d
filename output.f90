#include "macros.h"
module output
  use utils, only : newunit
  use simulation
  use Hamiltonian

contains

  subroutine output_lattice(this)
    type(Lattice), intent(in) :: this
  end subroutine output_lattice

  subroutine output_lattice_mean(this,new_sim)
    type(Lattice), intent(in) :: this
    logical, intent(in), optional :: new_sim
    
    integer, save :: u_mean
    logical :: o
    real(dl) :: vol
    integer :: n

    inquire(opened=o,file='log.out')
    if (.not.o) then
       open(unit=newunit(u_mean),file='log.out')
       call write_model_header(u_mean)
       call write_lattice_header(this,u_mean)
       write(u_mean,*) "# Time  <Psi_1>  <Psi_2>"
    endif

    if (present(new_sim).and.new_sim) write(u_mean,*)
    n = this%nlat
    vol = dble(this%nlat)*dble(this%nlat)

    write(u_mean,*) this%time, sum(this%psi(:,:,1))/vol, sum(this%psi(:,:,2))/vol
  end subroutine output_lattice_mean
    
end module output
