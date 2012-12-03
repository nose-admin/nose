!
! List of available modules and 
!   
module modules
    
    use module_tdpt3
    use module_qme
    use module_montecarlo
    use module_molc
    use module_test
    
    implicit none
    
    private              :: used_module
    character(len=64)    :: used_module
    
contains
    
    !
    ! Calls appropriate module
    !
    subroutine do_module_work(modname,err)
        character(len=*), intent(in) :: modname
        integer, intent(out)         :: err
               
        err = 0
        
        if (modname == "TDPT-3") then
            
            used_module = modname
            call do_tdpt3_work(err)
          
        else if (modname == "QME") then
            used_module = modname
            call do_qme_work(err)
            
        else if (modname == "MC") then
            used_module = modname
            call do_montecarlo_work(err)

        else if (modname == "TEST") then
            used_module = modname
            call do_test_work(err)

        else if (modname == "MOLC") then
            used_module = modname    
            call do_molc_work(err)
        else 
        
            print *, "Error: unknown module name ", modname
            err = -1
          
        end if
        
    end subroutine do_module_work

    !
    ! Collects module data from parallel calculation
    !
    subroutine collect_module_data(modname,err)
        character(len=*), intent(in) :: modname
        integer, intent(out)         :: err
               
        err = 0
        
        if (modname == "TDPT-3") then
            
            used_module = modname
            call collect_tdpt3_data(err)
          
        else if (modname == "QME") then
            
            used_module = modname
            call collect_qme_data(err)

        else if (modname == "MC") then

            used_module = modname
            call collect_montecarlo_data(err)

        else if (modname == "TEST") then

            used_module = modname
            call collect_test_data(err)

        else if (modname == "MOLC") then
            used_module = modname    
            !call do_molc_work(err)
            
        else 
        
            print *, "Error: unknown module name ", modname
            err = -1
          
        end if
                
    end subroutine collect_module_data

    !
    ! Cleans all module variables
    !
    subroutine clean_module(err)
        integer, intent(out)         :: err
               
        err = 0
        if (used_module == "TDPT-3") then
            
            call clean_tdpt3
            
        else if (used_module == "QME") then
            
            call clean_qme

        else if (used_module == "MC") then
            
            !call clean_montecarlo

        else if (used_module == "TEST") then
!            call clean_test
            
        else if (modname == "MOLC") then
!            call clean_molc

        else
        
            err = -1
                
        end if
    end subroutine clean_module
    
end module modules
