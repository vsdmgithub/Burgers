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
! MODULE: initial_condition
! LAST MODIFIED: 10 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR BURGERS EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE initial_condition
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! In this module, various initial conditions are provided in spectral space and in real space.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE system_parameters
    USE fft
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    IMPLICIT  NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	INTEGER (KIND=4):: 
    ! ---------------------------------------------------------
!    DOUBLE PRECISION::
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS

    SUBROUTINE make_initial_condition(en0,vel0_k)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Initialize initial condition in spectral space
    ! INPUT : Energy 
    ! OUTPUT : Complex array "vel_k" 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION::A0
        DOUBLE PRECISION,DIMENSION(0:N-1)::vel0_x
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION,INTENT(IN)::en0
        DOUBLE COMPLEX,DIMENSION(0:Nh),INTENT(OUT)::vel0_k

        A0=one
        ! Normalization constant
        
        DO x_ind = 0, N-1

            vel0_x ( x_ind ) = A0 * DSIN( x_axis( x_ind ) )
            
        END DO

        CALL energy_system(vel0_x, A0)
        ! Returns energy in real space for the array.
        
        A0 = DSQRT( en0 / A0 )
        ! Normalization factor

        DO x_ind = 0, N-1

            vel0_x ( x_ind ) = A0 * DSIN( x_axis( x_ind ) )

        END DO

        CALL fft_r2c(vel0_x, N, Nh, vel0_k)
        ! FFT to get spectral velocity

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END

    SUBROUTINE make_initial_potential
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! To get potential and its periodic counterpart from its spectral velocity
    ! make its periodic counterpart, along with initial lagrangian map (which is identity)
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! ---------------------------------------------------------
        ! FINDING THE SPECTRAL POTENTIAL FROM SPECTRAL VELOCITY
        ! ---------------------------------------------------------
        DO k_ind = 1, Nh
            vel_k( k_ind )     =    i * vel_k ( k_ind ) / k_axis ( k_ind )
        END DO
            vel_k ( 0 )        =    zero
        
        CALL fft_c2r(vel_k,Nh,N,potential_fn)
        ! FFT TO GET REAL POTENTIAL
        
        ! ---------------------------------------------------------
        ! MAKING A COPY OF INITIAL CONDITION ON PEROIDIC POTENTIAL
        ! ---------------------------------------------------------
        DO x_ind    = -N + 1, -1
            potential_fn_per0( x_ind )  = potential_fn( x_ind + N )
        END DO
        
        DO x_ind    = 0, N - 1
            potential_fn_per0( x_ind )  = potential_fn( x_ind )
        END DO
            potential_fn_per0( N )      = potential_fn( 0 )
            potential_fn_per0(-N )      = potential_fn( 0 )


        DO x_ind    = -N , N
            a_map(x_ind)   =    x_ind
        END DO
        
        potential_fn_per    =   potential_fn_per0
        ! INITIAL CONDITION FOR POTENTIAL FUNCTION
    END
END MODULE initial_condition
