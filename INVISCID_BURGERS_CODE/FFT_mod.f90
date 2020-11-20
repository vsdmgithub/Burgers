module fft
	! HEADER FILES/MODULES INCLUSION
    ! ----------------------
    use,intrinsic::iso_c_binding ! Standard module which defines the equivalent of C types in fortran
    implicit none
    include 'fftw3.f03'  ! Fortran interface files for all of the C routines for FFTW operation
    ! C VARIABLES DECLARATION
    ! -----------------------
    integer(C_INT)::N0
    type(C_PTR)::plan_r2c,plan_c2r,cdata_r2c_in,cdata_r2c_out,cdata_c2r_in,cdata_c2r_out ! all fftw plans are of this datatype in FORTRAN
    complex(C_DOUBLE_COMPLEX),pointer::data_r2c_out(:)
    real(C_DOUBLE),pointer::data_r2c_in(:)
    complex(C_DOUBLE_COMPLEX),pointer::data_c2r_in(:)
    real(C_DOUBLE),pointer::data_c2r_out(:)
    contains
    subroutine fft_r2c(data_in,M,Mh,data_out)
    ! Call this with real input array 'data_in' and get spectral output array 'data_out'
        implicit none
        integer(kind=4),intent(in)::M,Mh
        double precision,dimension(M),intent(in)::data_in
        double complex,dimension(0:Mh),intent(out)::data_out
        integer(kind=4)::k
        double precision::M_d
        M_d=DBLE(M) ! Normalization factor in the r2c transform
        N0=M
        ! ALLOCATE ARRAYS - DYNAMIC
        ! NOTE:- The array dimensions are in reverse order for FORTRAN
        ! ---------------------------------------
        cdata_r2c_in=fftw_alloc_real(int(N0,C_SIZE_T))
        call c_f_pointer(cdata_r2c_in,data_r2c_in,[N0])
        cdata_r2c_out=fftw_alloc_complex(int((N0/2+1),C_SIZE_T))
        call c_f_pointer(cdata_r2c_out,data_r2c_out,[(N0/2+1)])
        ! PLAN FOR OUT-PLACE FORWARD DFT R2C
        ! -----------------------------------
        plan_r2c=fftw_plan_dft_r2c_1d(N0,data_r2c_in,data_r2c_out,FFTW_MEASURE)
        ! INITIALIZE INPUT DATA
        ! ---------------------
        data_r2c_in=data_in
         ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_r2c(plan_r2c,data_r2c_in,data_r2c_out)
        ! WRITE OUTPUT (in format of first Brillouin zone format)
        ! -----------
		data_out=data_r2c_out/M_d
        ! DESTROY PLANS
        ! -------------
        call fftw_destroy_plan(plan_r2c)
        call fftw_free(cdata_r2c_in)
        call fftw_free(cdata_r2c_out)
    end
    subroutine fft_c2r(data_in,Mh,M,data_out)
    ! Call this with complex input array 'in' and get real output array 'out'
        implicit none
        integer(kind=4),intent(in)::M,Mh
        double complex,dimension(0:Mh),intent(in)::data_in
        double precision,dimension(M),intent(out)::data_out
        integer(kind=4)::k
        double precision::M_d
        N0=M
        M_d=1.0D0 ! Normalization factor in the FFT
        ! ALLOCATE ARRAYS - DYNAMIC
        ! NOTE:- The array dimensions are in reverse order for FORTRAN
        ! ---------------------------------------
        cdata_c2r_in=fftw_alloc_complex(int((N0/2+1),C_SIZE_T))
        call c_f_pointer(cdata_c2r_in,data_c2r_in,[(N0/2+1)])
        cdata_c2r_out=fftw_alloc_real(int(N0,C_SIZE_T))
        call c_f_pointer(cdata_c2r_out,data_c2r_out,[N0])
        ! PLAN FOR OUT-PLACE FORWARD DFT R2C
        ! -----------------------------------
        plan_c2r=fftw_plan_dft_c2r_1d(N0,data_c2r_in,data_c2r_out,FFTW_MEASURE)    
        ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
        ! ---------------------
		data_c2r_in=data_in
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_c2r(plan_c2r,data_c2r_in,data_c2r_out)
        ! WRITE OUTPUT
        ! ------------
        data_out=data_c2r_out/M_d
        ! DESTROY PLANS
        ! -------------
        call fftw_destroy_plan(plan_c2r)
        call fftw_free(cdata_c2r_in)
        call fftw_free(cdata_c2r_out)
    end
end module fft
