module main_comm

    ! libnose
    use std_types
    use std_parallel
    use util_timing
    use util_char_graphics

    ! local
    use modules
    use prepare

    implicit none

    !
    character(len=256) :: part_name

contains

!
! OUTPUT ORGANIZATION
!
    subroutine warning_parallel()
        if (logging) then
            call print_log_message(" ",5)
            call print_log_message("            ***********  Starting parallel processing ***********",5)
        end if
    end subroutine warning_parallel

    subroutine warning_output()
        if (logging) then
            call print_log_message("              Following output is from the master process only",5)
            call print_log_message(" ",5)
        end if
    end subroutine warning_output

    subroutine start_section(pname)
        character(len=*), intent(in) :: pname
        call print_log_message(pname,5)
        part_name = pname
        call flush()
    end subroutine start_section

    subroutine start_section_parallel(pname)
        character(len=*), intent(in) :: pname
        if (logging) then
            call print_log_message(pname,5)
            part_name = pname
            call flush()
        end if
    end subroutine start_section_parallel

    subroutine next_section(spec)
        character(len=*), intent(in) :: spec
        call print_log_message(trim(part_name)//"- section finished",5)
        call next_code_position()
        call print_elapsed_time(elapsed_time(spec))
        call print_devider(2)
        call flush()
    end subroutine next_section

    subroutine next_section_parallel(spec)
        character(len=*), intent(in) :: spec
        character(len=256) :: char_buff
        if (logging) then
            write(char_buff,'(2a)') trim(part_name), ": section finished"
            call print_log_message(trim(char_buff),5)
            call next_code_position()
            call print_elapsed_time(elapsed_time(spec))
            call print_devider(2)
            call flush()
        end if
    end subroutine next_section_parallel

    subroutine finish_sections(spec)
        character(len=*), intent(in) :: spec
        call print_devider(1)
        call print_log_message("Execution summary", 5)
        call next_code_position()
        call print_total_elapsed_time(elapsed_time(spec,.false.))
        call print_devider(2)
        call flush()
    end subroutine finish_sections

    subroutine finish_sections_parallel(spec)
        character(len=*), intent(in) :: spec
        if (logging) then
            call print_devider(1)
            call print_log_message("Execution summary", 5)
            call next_code_position()
            call print_total_elapsed_time(elapsed_time(spec,.false.))
            call print_devider(2)
            call flush()
        end if

    end subroutine finish_sections_parallel

    subroutine print_header()
        integer :: n
        n = 1
        call print_log_message(" ",n)
        call print_log_message(" ",n)
        call print_log_message(&
"                           NOSE (ver. 0.5.33) ",n)
        call print_log_message(" ",n)
        call print_log_message(&
"                 Institute of Physics of Charles University,",n)
        call print_log_message(&
"                    Faculty of Mathematics and Physics,",n)
        call print_log_message(&
"                       Charles University in Prague",n)
!        call print_log_message(&
!"                                   & ",n)
!        call print_log_message(&
!"                       Department of Chemistry",n)
!        call print_log_message(&
!"                  University of California, Berkeley",n)
!        call print_log_message(&
!"                                   & ",n)
!        call print_log_message(&
!"              Physical Biosciences Division, LBNL, Berkeley",n)
        call print_log_message(" ",n)
        call print_log_message(&
"                        last change 10-15-2011 ",n)
        call print_log_message(" ",n)
    end subroutine print_header

    subroutine print_header_parallel()
        if (logging) then
            call print_header()
        end if
    end subroutine print_header_parallel

    subroutine print_final_parallel()
        if (logging) then
            call print_log_message(" ",5)
            call print_log_message(" This is the end, my friend ... ",5)
            call print_log_message(" ... bye!",5)
            call print_log_message(" ",5)
        end if
    end subroutine print_final_parallel

end module main_comm
