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
    INTEGER(KIND=4)::t_step_purging
    INTEGER(KIND=4)::k_P
    ! _________________________
    DOUBLE PRECISION::purging_alpha,purging_beta
    DOUBLE PRECISION::time_purging
    DOUBLE PRECISION::initial_en
    DOUBLE PRECISION::U_rms
    DOUBLE PRECISION::time_grid
    ! _________________________
    CHARACTER(LEN=30)::name_sys
    ! _________________________
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::purger
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vel_x
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

        ALLOCATE( purger (0 : Nh) )
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        name_sys        =   'purged_'

        purging_alpha   =   0.8D0

        purging_beta    =   0.6D0

        k_P             =   k_G - FLOOR( DBLE(k_G) ** purging_beta )

        time_purging    =   one / ( DBLE(k_G) ** purging_alpha )

        initial_en      =   one

        U_rms           =  DSQRT( two *  initial_en )

        time_grid       = 0.1D0 *  dx / U_rms

        DO k_ind = 0, Nh
            IF ( k_axis(k_ind) .LE. k_P ) THEN
               
                purger (k_ind)  =   one

            END IF
        END DO

        CALL time_to_step_convert(time_purging, t_step_purging)
        ! Convert purging time to no of time steps to purge
        
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

END MODULE system_parameters
