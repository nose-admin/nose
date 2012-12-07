module resources_montecarlo

    use resources
    use sci_redfield

    implicit none

 	logical	:: 	origin = .false., g_functions = .true., depository = .false., variant360 = .false., 	&
 										modified_unraveling2 = .false., fixed_seed = .false., 			&
 										vynechani_G_ifu = .false., use_exciton_basis = .false., 		&
 										load_evops = .false.
 	real(dp)	::  debug_gamma = 1e-10 ! nonzero for debug only

 	character(len=256) :: methodMC



contains

    !
    ! Initialization of the MC specific variables
    !
    subroutine init_resources_montecarlo()


    end subroutine init_resources_montecarlo


    !
    ! Cleaning the QME specific variables
    !
    subroutine clean_resources_montecarlo()



    end subroutine clean_resources_montecarlo


end module resources_montecarlo
