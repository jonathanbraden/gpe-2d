module fftw3_cmplx
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none
  include 'fftw3.f03'

  complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)

  type transformPair2D_c
     integer :: nx, ny, nnx, nny
     complex(C_DOUBLE_COMPLEX), pointer :: realSpace(:,:)
     complex(C_DOUBLE_COMPLEX), pointer :: specSpace(:,:)
     type(C_PTR) :: planf, planb
     type(C_PTR), private :: rPtr, sPtr
  end type transformPair2D_c

contains

  subroutine boot_openmp(nThread)
    integer, optional, intent(in) :: nThread
    integer :: errorOMP
#ifdef USEOMP
    errorOMP = fftw_init_threads()
    if (errorOMP == 0) then
       print*,"Error initializing OpenMP threading for FFTW"
       stop
    endif
    if (present(nThread)) then
       errorOMP = nThread
    else
       print*,"Defaulting to OMP_NUM_THREADS environment variable"
       errorOMP = omp_get_max_threads()
    endif
    call fftw_plan_with_nthreads(errorOMP)
    print*,"FFTW booted using ",errorOMP," threads"
#endif
  end subroutine boot_openmp

  subroutine initialize_transform_2d_c(this,n)
    type(transformPair2D_c), intent(out) :: this
    integer, intent(in), dimension(1:2) :: n

    this%nx = n(1); this%ny = n(2)
    this%nnx = n(1)/2+1; this%nny = n(2)/2+1  ! Check these

    call allocate_2d_array_c(n(1),n(2),this%realSpace,this%specSpace,this%rPtr,this%sPtr)
    this%planf = fftw_plan_dft_2d(n(2),n(1),this%realSpace,this%specSpace,FFTW_FORWARD,FFTW_MEASURE)
    this%planb = fftw_plan_dft_2d(n(2),n(1),this%specSpace,this%realSpace,FFTW_BACKWARD,FFTW_MEASURE) ! check ordering in this one
  end subroutine initialize_transform_2d_c

  subroutine destroy_transform_2d_c(this)
    type(transformPair2D_c), intent(inout) :: this

    call fftw_destroy_plan(this%planf); call fftw_destroy_plan(this%planb)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)
  end subroutine destroy_transform_2d_c

  subroutine allocate_2d_array_c(nx,ny,arr,Fk,fptr,fkptr)
    integer :: nx,ny
    complex(C_DOUBLE_COMPLEX), pointer :: arr(:,:), Fk(:,:)
    type(C_PTR) :: fptr, fkptr

    fptr = fftw_alloc_complex(int(nx*ny,C_SIZE_T))
    call c_f_pointer(fptr, arr, [nx,ny])
    fkptr = fftw_alloc_complex(int(nx*ny,C_SIZE_T))
    call c_f_pointer(fkptr,Fk, [nx,ny])
  end subroutine allocate_2d_array_c

  subroutine forward_transform_2d_wtype_c(this)
    type(transformPair2D_c), intent(inout) :: this
    call fftw_execute_dft(this%planf,this%realSpace,this%specSpace)
  end subroutine forward_transform_2d_wtype_c

  subroutine backward_transform_2d_wtype_c(this)
    type(transformPair2D_c), intent(inout) :: this
    call fftw_execute_dft(this%planb,this%specSpace,this%realSpace) ! Check ordering
  end subroutine backward_transform_2d_wtype_c
  
  subroutine laplacian_2d_wtype_c(this,dk)
    type(transformPair2D_c), intent(inout) :: this
    real(dl), intent(in) :: dk

    integer :: i,j,jj; real(dl) :: rad2
#ifdef VECTORIZE
    integer, dimension(1:this%nx) :: i_x

    do i=1,this%nnx
       i_x(i) = i-1
       i_x(this%nnx+i) = this%nx+1-i
    enddo
#endif
    call fftw_execute_dft(this%planf,this%realSpace,this%specSpace)
    rad2 = 0._dl
!$OMP PARALLEL DO FIRSTPRIVATE(dk) PRIVATE(i,j,jj,rad2)
    do j=1,this%ny; if (j>this%nny) then; jj = this%ny+1-j; else; jj=j-1; endif
#ifdef VECTORIZE
       this%specSpace(:,j) = -(i_x**2+jj**2)*dk**2*this%specSpace(:,j)
#else
       do i=1,this%nnx
          rad2 = dble((i-1)**2+jj**2)
          this%specSpace(i,j) = -rad2*dk**2*this%specSpace(i,j)
       enddo
       do i=this%nnx+1,this%nx
          rad2 = dble((this%nx+1-i)**2 + jj**2)
          this%specSpace(i,j) = -rad2*dk**2*this%specSpace(i,j)
       enddo
#endif
    enddo
!$OMP END PARALLEL DO
    call fftw_execute_dft(this%planb,this%specSpace,this%realSpace)
    this%realSpace = this%realSpace / dble(this%nx) / dble(this%ny)
  end subroutine laplacian_2d_wtype_c

  subroutine gradient_energy_parts_c(this)
    type(transformPair2D_c), intent(inout) :: this
  end subroutine gradient_energy_parts_c

  subroutine gradient_energy_c(this)
    type(transformPair2D_c), intent(inout) :: this
  end subroutine gradient_energy_c
  
end module fftw3_cmplx
