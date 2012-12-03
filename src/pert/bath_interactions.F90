#include "util_allocation.h"

!
! This modules prepares all variables associated with bath interaction
! using inputs from wizards
!
module bath_interactions

    use goft
    use excitons
    use rates
    !use input_wizard

    implicit none

contains


    !
    !
    !
    subroutine prepare_bath_interactions(err)
        integer, intent(out) :: err
         ! local
        real(dp), dimension(:), allocatable :: ll
        integer :: i

        err = 0

		call print_log_message("Bath interactions ...", 8)

        !
        ! Calculate g(t) of the exciton states
        !
        call goft_prepare_excitonic_goft()

        !
        ! Calculate relaxation rates
        !
        call init_rates(err)

		call print_log_message("... bath interactions done.", 8)

    end subroutine prepare_bath_interactions

end module bath_interactions
