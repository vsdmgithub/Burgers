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
! MODULE: global_variables
! LAST MODIFIED: 10 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! GLOBAL VARIABLES FOR BURGERS EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE global_variables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temporary (IF necessary) are declared within the subroutines.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE constants
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
	IMPLICIT  NONE
    ! _________________________
    ! SPACE VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::N,Nh,N_log_2
	INTEGER (KIND=4)::x_ind,k_ind    
    ! ---------------------------------------------------------
    DOUBLE PRECISION::N_db
    DOUBLE PRECISION::length,dx
    ! _________________________
    ! TIME (SIMULATION) VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::t_step,t_step_total
	INTEGER (KIND=4)::t_step_save
	INTEGER (KIND=4)::t_step_debug
	INTEGER (KIND=4)::save_no,save_total    
    ! ---------------------------------------------------------
    DOUBLE PRECISION::time_total,time_now,time_save
    DOUBLE PRECISION::dt
    ! _________________________
    ! OUTPUT AND INPUT FILES AND DIR
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER (KIND=4)::state_sim
    ! _________________________
    DOUBLE PRECISION::rd_no
    ! _________________________
    CHARACTER(LEN=60)::input_file
    CHARACTER(LEN=80)::name_sim
    CHARACTER(LEN=60)::N_char
    CHARACTER(LEN=60)::path_dir
    CHARACTER(LEN=60)::name_dir
    ! _________________________
    ! GLOBAL ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::x_axis,k_axis,t_axis
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
    SUBROUTINE read_input
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Read simulation parameters from a file 'input_file'
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE

        input_file='simulation_parameters.dat'
        ! This file contains all major Input parameters to be fed
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        OPEN(unit=1001,file=TRIM(ADJUSTL(input_file)))
        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_i4,ADVANCE='yes'),N_log_2 	        ! Resolution of real domain in log_2 scale.
        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_d8p4,ADVANCE='yes'),dt                  ! Time step for marching forward in time
        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_d8p4,ADVANCE='yes'),time_total       ! Total time to simulate
        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_i4,ADVANCE='yes'),save_total             ! No of saves for the velocity data. 
        CLOSE(1001)
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END

	SUBROUTINE init_global_variables
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the global variables that are used through out the code.
    ! PREREQUISITE: SUBROUTINE 'read_input' has to be called before initializing  these variables.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE

        state_sim=0
        ! This being the first variable to start the simulation. At last, it will be set to '1'

        CALL get_simulation_name(name_sim)
        ! Creating dated and timed name for the simulation
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! S  P  A  C  E       A  N  D         T  I  M  E  
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        length      =       two_pi
        ! Length of the periodic box (natural units)

        N       =   2 ** ( N_log_2 )
        ! Resolution

        N_db    =   DBLE( N )
        ! Double formatted 'N'

        Nh     =    N / 2
        ! Size of the fourier transform
      
        dx      =  length / N_db
        ! Grid Gap

        dt      = DBLE( FLOOR(  100000000.0D0 * dt ) ) / 100000000.D0
        ! Making sure that it doesnot have any error, so that time_steps comes nicely as
        ! 3.0000 rather than 2.998
        
        CALL time_to_step_convert(time_total,t_step_total)
        ! Converts time to steps
        
    	t_step_save   =     t_step_total / save_total
        ! Determines how many time steps after the save has to be made.

        CALL step_to_time_convert(t_step_save,time_save)
        ! Converts steps to time
        
    	t_step_debug =      t_step_total / 10
        ! No of times, the data will be checked for 'NaN' during the simul

        WRITE (N_char, f_i8),N
        ! converting resolution value to CHARACTER

        path_dir    =   '../burgers_data/'
        ! path of the main directory relative to this file.

        name_dir    =   'N'//TRIM(ADJUSTL(N_char))//'/'
        ! name of the main directory
        
	END

	SUBROUTINE init_global_arrays
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the global arrays that are declared/allotted here.
    ! PREREQUISITE: SUBROUTINE 'init_global_variables' has to be called before initializing  these arrays.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        
        IMPLICIT  NONE

        ALLOCATE ( x_axis(0 : N-1), k_axis(0 : Nh) , t_axis(0 : t_step_total) )
        
        DO x_ind = 0, N-1

            x_axis ( x_ind )    =  x_ind * dx

        END DO

        DO k_ind = 0, Nh

            k_axis ( k_ind )    =  k_ind * one

        END DO

        DO t_step = 0 , t_step_total

            t_axis ( t_step )   =    dt *  t_step
              
        END DO

    END

    SUBROUTINE step_to_time_convert(step,time)
	! CALL this to convert time step into actual time of simulation

    	IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(IN)::step
		DOUBLE PRECISION,INTENT(OUT)::time
		time    =   DBLE(step) * dt

    END

    SUBROUTINE time_to_step_convert(time,step)
	! CALL this to convert time step into actual time of simulation

    	IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(OUT)::step
		DOUBLE PRECISION,INTENT(IN)::time
        step    =   FLOOR(time / dt)

    END

    SUBROUTINE get_simulation_name( sim_char )
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This provides a string  of format
    !     'run_d121120_t104022' for run dated 12/11/2020 timed 10:40:22 
    !  
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=8)::year_char,month_char,date_char
        CHARACTER(LEN=8)::hour_char,min_char,sec_char
        INTEGER,DIMENSION(8)::values
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=*),INTENT(OUT)::sim_char
 
        CALL DATE_AND_TIME(VALUES=values)
        ! Gets the date and time as integer array
        
        values(1)   =   MOD(values(1),2000)
        ! writing only last two digits of year
        
        WRITE(year_char,f_i4), values(1)
        WRITE(month_char,f_i2), values(2)
        WRITE(date_char,f_i2), values(3)
        WRITE(hour_char,f_i2), values(5)
        WRITE(min_char,f_i2), values(6)
        WRITE(sec_char,f_i2), values(7)
        ! Self-explained

        sim_char    =  'run_d'//TRIM(ADJUSTL(date_char))//TRIM(ADJUSTL(month_char))//&
        TRIM(ADJUSTL(year_char))//'_t'//TRIM(ADJUSTL(hour_char))//&
        TRIM(ADJUSTL(min_char))//TRIM(ADJUSTL(sec_char))

    END
END MODULE global_variables
