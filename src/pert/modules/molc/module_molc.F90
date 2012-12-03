!
! Module MOLC
!

module module_molc

    use excband
    
    implicit none 

    public  :: do_molc_work

contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_molc_work(err)
        integer, intent(out) :: err

        write(*,*) "MOLC working"
        call excband_structure()

    end subroutine do_molc_work

end module module_molc

