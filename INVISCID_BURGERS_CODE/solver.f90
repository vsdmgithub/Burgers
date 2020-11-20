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
    SUBROUTINE fast_legendre_transform
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Using potential method, from the saddle point evaluation of the convolution in
    ! integral solution of Cole-Hopf transformed burgers. This method uses Legendre transformation
    ! which numerically costs NxN operations. Using a regularity in the transformation we can reduce
    ! it to N Log(N) operations, called Fast Legendre Transformation
    ! ------------------------------------------------------------------------------------------
    ! Since we have a periodic function, it is sufficient to take two periods of interval to do it. we
    ! from x=-N to x=N. Finally translating x=-N/2,0 to x=N/2,N replacing all the data.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    IMPLICIT NONE
    
    ! -------------------
    ! FOR X=0 INDEX
    ! -------------------
    x_lef_ind   =   -N
    x_rig_ind   =    N
    x_ind       =    0

    CALL find_saddle_point
    
    a_map ( x_ind )             =   a_map ( x_lef_ind ) + a_saddle - 1
    potential_fn_per( x_ind )   =   potential_saddle
    
    DO s_lab=0,m_lab-1
    ! S RUNS THROUGH LOG_2(N) LOOPS, SPLITTING THE ARRAY INTO 2^S SUBARRAYS EVERYTIME.

        DO j_lab=1,2**s_lab
        ! THE LOOP THAT RUNS THROUGH 2^S SUBARRAYS. EACH J IS IN A SUBARRAY OF SIZE N/2^S

            ! -------------------
            ! PART FROM X=0 TO X=N
            ! -------------------
            x_lef_ind       =   (j_lab-1) * ( 2 ** (m_lab-s_lab) )
            ! DETERMINNG THE LEFT MOST ELEMENT OF THE SUBARRAY
            x_rig_ind       =   j_lab * ( 2 ** (m_lab-s_lab) )
            ! DETERMINING THE RIGHT MOST ELEMENT OF THE SUBARRAY
            x_ind           =   ( x_lef_ind + x_rig_ind ) / 2
            ! THE ELEMENT OF THE SUBARRAY, FOR WHICH THE LEGENDRE TRANSFORMATION HAS TO BE CARRIED             

            CALL find_saddle_point                
            ! FINDS THE INDEX OF THE SADDLE POINT BETWEEN LEFT AND RIGHT INDICES FOR THE SOLUTION TO POTENTIAL

            a_map ( x_ind )             = a_map ( x_lef_ind ) + a_saddle - 1
            ! THE LOCATION WHERE MINIMUM POTENTIAL IS OBTAINED, ALSO THE INVERSE LAGRANGIAN IMAGE OF 'X'
            potential_fn_per ( x_ind )  = potential_saddle
            ! ASSIGNING TO X, THE MINMIUM EFFECTIVE POTENTIAL DETERMINED BY THE SADDLE POINT APPX TO INTEGRAL

            ! -------------------
            ! PART FROM X=-N TO X=0
            ! -------------------
            x_swap          =     x_lef_ind
            ! TEMPORARY VARIABLE TO SWAP LEFT AND RIGHT INDICES
            x_lef_ind       =   - x_rig_ind
            ! DETERMINING THE LEFT MOST ELEMENT OF THE SUBARRAY
            x_rig_ind       =   - x_swap
            ! DETERMINNG THE RIGHT MOST ELEMENT OF THE SUBARRAY
            x_ind           =   ( x_lef_ind + x_rig_ind ) / 2
            ! THE ELEMENT OF THE SUBARRAY, FOR WHICH THE LEGENDRE TRANSFORMATION HAS TO BE CARRIED             

            CALL find_saddle_point
            ! FINDS THE INDEX OF THE SADDLE POINT BETWEEN LEFT AND RIGHT INDICES FOR THE SOLUTION TO POTENTIAL
            
            a_map ( x_ind )             = a_map ( x_lef_ind ) + a_saddle - 1

            ! THE LOCATION WHERE MINIMUM POTENTIAL IS OBTAINED, ALSO THE INVERSE LAGRANGIAN IMAGE OF 'X'
            potential_fn_per ( x_ind )  = potential_saddle
            ! ASSIGNING TO X, THE MINMIUM EFFECTIVE POTENTIAL DETERMINED BY THE SADDLE POINT APPX TO INTEGRAL
            
        END DO
        
    END DO
    
    !-----------------------
    ! RETAINING PERIODICITY
    !-----------------------
    potential_fn_per( N )       = potential_fn_per( 0 )
    potential_fn_per(-N )       = potential_fn_per( 0 )
  
    !----------------------------------------------------
    ! EXTRACTING DATA FROM 0-N FOR PSI FROM BOTH SIDES 0
    !----------------------------------------------------
    DO  x_ind   = -Nh , -1
        potential_fn ( x_ind + N )      =   potential_fn_per ( x_ind )
        potential_fn_per ( x_ind + N )  =   potential_fn_per ( x_ind )
    END DO
    DO x_ind    = 0 , Nh-1
        potential_fn ( x_ind )          =   potential_fn_per ( x_ind )
        potential_fn_per ( x_ind - N )  =   potential_fn_per ( x_ind )
    END DO

    !---------------------------------------------------
    ! UPDATING THE POTENTIAL FOR NEXT STEP
    !---------------------------------------------------
!    potential_fn_per0        =   potential_fn_per

    !---------------------------------------------------
    ! SAVING THE STEPWISE LAGRANGIAN MAP TO A GLOBAL MAP
    !---------------------------------------------------
    a_map_global (t_step,:)  =   a_map

    END
    
    SUBROUTINE find_saddle_point 
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! The potential solution at time 't' is known as the integral. In the viscosity going to zero
    ! the limit can be taken as saddle point approximation.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    a_lef_ind    =     a_map ( x_lef_ind )
    a_rig_ind    =     a_map ( x_rig_ind )
    ! FINDING THE INVERSE LAGRANGIAN MAP OF THE SUBARRAY

    ALLOCATE(sub_array (1 : a_rig_ind - a_lef_ind + 1) )
    ! ALLOTING THE SIZE OF SUBARRAY TO THE FIND ITS SADDLE POINT

    DO y_ind = a_lef_ind, a_rig_ind
        sub_array( y_ind - a_lef_ind + 1 ) =  potential_fn_per0( y_ind )-(( DBLE(y_ind-x_ind) * dx )**two) / (two*(time_now + dt) )
    END DO
    
    a_saddle           = MAXLOC( sub_array , DIM = 1 )
    ! FINDS THE MAX LOCATION IN THE SUBARRAY

    potential_saddle   = sub_array( a_saddle )
    ! THIS IS THE ANSWER TO THE POTENTIAL AT THE GIVEN POINT INSIDE THE SUBARRAY

    DEALLOCATE ( sub_array )
    
    END

  END MODULE solver
