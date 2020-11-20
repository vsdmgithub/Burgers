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
! MODULE: system_parameters
! LAST MODIFIED: 10 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SYSTEM PARAMETERS FOR BURGERS EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE system_parameters
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the system parameters, corresponding to the different variant of burgers equation we are studying
! is defined here. Note global variables is common amongst these different variants.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE global_variables
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
	IMPLICIT  NONE
    ! _________________________
    ! SYSTEM VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::x_lef_ind,a_lef_ind
    INTEGER(KIND=4)::x_rig_ind,a_rig_ind
    INTEGER(KIND=4)::a_saddle
    INTEGER(KIND=4)::s_lab,m_lab,j_lab
    INTEGER(KIND=4)::x_swap,y_ind
    ! _________________________
    DOUBLE PRECISION::initial_en
    DOUBLE PRECISION::U_rms
    DOUBLE PRECISION::potential_saddle
    DOUBLE PRECISION::time_grid
    ! _________________________
    CHARACTER(LEN=30)::name_sys
    ! _________________________
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vel_x
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::potential_fn_per0
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::potential_fn_per
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::potential_fn,sub_array
    INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE::a_map,a_map_net
    INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE::a_map_global
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_time 
    DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE::vel_k
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
	SUBROUTINE init_system_parameters
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the system related parameters.
    ! PREREQUISITE: SUBROUTINE 'init_global_arrays' (in global_variables module)
    !  has to be called before initializing  these variables.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        name_sys    =   'inviscid_'

        initial_en  =   one

        U_rms      =  DSQRT( two *  initial_en )

        time_grid  = 0.1D0 *  dx / U_rms

        m_lab      = N_log_2

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          
	END

    SUBROUTINE energy_system(vel0_x,en0)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Find the energy in real space
    ! INPUT : Velocity array 
    ! OUTPUT : Energy of it 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION,INTENT(OUT)::en0
        DOUBLE PRECISION,DIMENSION(0:N-1),INTENT(IN)::vel0_x

        en0  =   SUM( vel0_x ** two)
        en0  =   hf * en0 / N_db
        
     END
    
    SUBROUTINE gradient_potential 
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! To get velocity from its potential.
    ! Returns central first order difference, error of (\delta x)^2
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! -------------------------
        ! CENTRAL DIFFERENCE METHOD
        ! -------------------------
        DO x_ind=1,N-2
            vel_x( x_ind )    =    potential_fn( x_ind + 1 ) - potential_fn( x_ind - 1 )
        END DO
        
        vel_x( 0 )      =   potential_fn( 1 )- potential_fn( N - 1 )
        vel_x( N - 1 )  =   potential_fn( 0 )- potential_fn( N - 2 )
        vel_x           =   - hf * vel_x / dx
        ! because u=-d\psi/dx

    END
    
    SUBROUTINE net_lagrangian_map
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Call this to find the Lagrangian map at time t relative to initial condition
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER(KIND=4)::t_lag

        a_map_net   =   a_map_global( t_step, : )

        ! ---------------------------------------------------------
        ! THIS IS COMPOSITION OF STEPWISE LAGRANGIAN MAPS 
        ! ---------------------------------------------------------
        DO t_lag = 1, t_step 
            DO x_ind = -N, N
                a_map_net( x_ind )  =   a_map_global( t_step - t_lag, a_map_net( x_ind ))
            END DO
        END DO
        
    END 

END MODULE system_parameters
