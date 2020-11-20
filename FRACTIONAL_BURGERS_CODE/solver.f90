! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------------------------------------------------------------------------------------                                                                           
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
! #########################
! MODULE: solver
! LAST MODIFIED: 10 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER FOR  BURGERS EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE solver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and updates it by a step, using the subroutines
! 1. rk4_algorithm
! 2. time_derivative
! 3. convection_dissipation_terms 
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE initial_condition
	USE fft
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    IMPLICIT NONE
    ! _________________________
    ! SOLVER ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE::dv1, dv2, dv3, dv4
    DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE::vel_k_temp
    DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE::convection_k
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
	SUBROUTINE rk4_algorithm
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
    ! Alg: - Runga kutta 4th order
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)''

        ALLOCATE(dv1(0:Nh), dv2(0:Nh), dv3(0:Nh), dv4(0:Nh))
        ALLOCATE(vel_k_temp(0:Nh),  convection_k(0:Nh))
        
        vel_k_temp  =   truncator   *   vel_k
        CALL time_derivative(dv1) ! This call provides \vec{dv} for the existing \vec{v}

    	vel_k   =   vel_k_temp  +   hf * dv1
        CALL time_derivative(dv2)

        vel_k   =   vel_k_temp  +   hf * dv2
        CALL time_derivative(dv3)
        
		vel_k   =   vel_k_temp  +   dv3
        CALL time_derivative(dv4)

        ! Final increment for 'v(k)'
		vel_k   =   vel_k_temp  +   ( dv1 + two * dv2 + two * dv3 + dv4 ) / six

        CALL fft_c2r(vel_k,Nh,N,vel_x)
        ! FFT to get real velocity

        DEALLOCATE( dv1, dv2, dv3, dv4)
        DEALLOCATE( vel_k_temp,  convection_k)

    END
    
    SUBROUTINE time_derivative(dv)
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to get the time derivative matrix for matrix 'v(k)'
        ! This is the BURGERS EQUATION implemented for numerical computation
        ! spectral space.
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        DOUBLE COMPLEX,DIMENSION(0 : Nh),INTENT(OUT)::dv

        CALL convection_dissipation_terms
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   B  U  R  G  E  R  S           E   Q   N.
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        dv  =   convection_k 
        dv  =   dt * truncator * dv

	END

    SUBROUTINE convection_dissipation_terms
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to give convection term by the prescribed function f(x), usually it is the energy u^2/2 and
        ! the dissipation term.
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   C   O   N   V   E   C   T   I   O   N       T   E   R   M
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        CALL fft_c2r(vel_k, Nh, N,  vel_x)
        
        vel_x  =  vel_x  * vel_x
        vel_x  =  - hf * vel_x
        ! using vel_x array to store the convection term in real space
        
        CALL fft_r2c(vel_x, N,  Nh, convection_k)
        ! Fourier transform of real energy to spectral energy

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   D    I    S   S   I   P   A   T   I   O   N       T   E   R   M
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        convection_k  = i  *  k_axis  *  convection_k
        ! Derivative of energy in spectral space
        
        convection_k  = convection_k  -  viscosity  *  frac_laplacian_k  *  vel_k
        ! Adding the dissipation terms to the convection
        
    END
END MODULE solver
