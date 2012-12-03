!
!  NOSE Input Stream
!  Implementation of a simple protocol of the NOSE communication
!  with the user layer and between processors
!
!  Tomas Mancal, tomas.mancal@mff.cuni.cz
!
!  last changed: 2007-01-16
!
module nis ! NOSE INPUT STREAM

    use std_types
    use std_io
    use std_parallel
    use numer_random

    implicit none

    !
    ! Public
    !

    ! this is the random generator seed
    integer, dimension(:), allocatable :: rseed

    ! constants
    integer, public,                 parameter :: NIS_CONST_LENGTH = 3
    character(len=NIS_CONST_LENGTH), parameter :: NIS_MODE_STD     = "STD"
    character(len=NIS_CONST_LENGTH), parameter :: NIS_MODE_MPI     = "MPI"
    character(len=NIS_CONST_LENGTH), parameter :: NIS_MODE_FILE    = "FLE"
    character(len=*), parameter :: NIS_DEFAULT_FILENAME = "nis_logfile"

    integer                                    :: NIS_IO_UNIT
    logical                                    :: output_file_set

    character(len=64) :: output_file
    character(len=64), private :: output_dir

    !
    ! interfaces
    !

    interface nis_add
        module procedure add_new_data_i11, add_new_data_r11, add_new_data_rmn, add_new_data_imn, add_new_data_string
    end interface

    interface nis_read_int
        module procedure nis_read_int_nm, nis_read_int_11
    end interface

    interface nis_read_real
        module procedure nis_read_real_nm, nis_read_real_11
    end interface

    interface nis_send_int
        module procedure nis_send_int_nm, nis_send_int_11
    end interface

    interface nis_send_real
        module procedure nis_send_real_nm, nis_send_real_11
    end interface

    interface nis_get_real
        module procedure get_real_m, get_real_v, get_real_s
    end interface

    interface nis_get_integer
        module procedure get_integer_m, get_integer_v, get_integer_s
    end interface

    !
    ! Private
    !


    ! nis operating mode
    character(len=NIS_CONST_LENGTH), private   :: nis_mode

    !
    ! input stream private types
    !

    type, private :: real_node
        integer                               :: index
        real(dp), dimension(:,:), pointer     :: data
        type(real_node), pointer              :: next
    end type real_node

    type, private :: complex_node
        integer                               :: index
        complex(dpc), dimension(:,:), pointer :: data
        type(complex_node), pointer           :: next
    end type complex_node

    type, private :: integer_node
        integer                               :: index
        integer, dimension(:,:), pointer      :: data
        type(integer_node), pointer           :: next
    end type integer_node

    type, private :: string_node
        integer                               :: index
        character, dimension(:), pointer      :: data
        type(string_node), pointer            :: next
    end type string_node

    type, private :: char_node
        integer                               :: index
        character                             :: data
        type(char_node), pointer              :: next
    end type char_node

    type(char_node), pointer, private      :: stream_types
    type(char_node), pointer, private      :: current_type
    type(integer_node), pointer, private   :: stream_integers
    type(integer_node), pointer, private   :: current_integer
    type(real_node), pointer, private      :: stream_reals
    type(real_node), pointer, private      :: current_real
    type(string_node), pointer, private    :: stream_strings
    type(string_node), pointer, private    :: current_string
    type(char_node), pointer, private      :: p_current
    type(integer_node), pointer, private   :: i_current
    type(real_node), pointer, private      :: r_current
    type(string_node), pointer, private    :: c_current

    private :: add_new_data_i11, add_new_data_r11
    private :: add_new_data_rmn, add_new_data_imn
    private :: add_new_data_string
    private :: nis_read_int_nm, nis_read_int_11
    private :: nis_read_real_nm, nis_read_real_11
    private :: nis_send_int_nm, nis_send_int_11
    private :: nis_send_real_nm, nis_send_real_11
    private :: get_real_m, get_real_v, get_real_s
    private :: get_integer_m, get_integer_s, get_integer_v

    character(len=256), private :: cbuff

contains

    !
    ! Initialization
    !
    subroutine init_nis

        ! types
        allocate(stream_types)
        stream_types%index = 0
        stream_types%data  = '0'
        nullify(stream_types%next)
        current_type => stream_types

        ! integers
        allocate(stream_integers)
        stream_integers%index = 0
        nullify(stream_integers%data)
        nullify(stream_integers%next)
        current_integer => stream_integers

        ! reals
        allocate(stream_reals)
        stream_reals%index = 0
        nullify(stream_reals%data)
        nullify(stream_reals%next)
        current_real => stream_reals

        ! string
        allocate(stream_strings)
        stream_strings%index = 0
        nullify(stream_strings%data)
        nullify(stream_strings%next)
        current_string => stream_strings

        allocate(p_current,i_current,r_current,c_current)
        p_current => stream_types
        i_current => stream_integers
        r_current => stream_reals
        c_current => stream_strings

        call nis_set_mode(NIS_MODE_STD)

        NIS_IO_UNIT = 50 + parallel_id
        output_file_set = .false.

    end subroutine init_nis

    !
    !
    !
    subroutine nis_open_file(un,fname)
        integer, intent(in) :: un
        character(len=*), intent(in) :: fname
        NIS_IO_UNIT = un
        output_file = fname
        output_file_set = .true.

    end subroutine nis_open_file

    !
    !
    !
    subroutine nis_close_file()
        integer :: ierr
        if (nis_mode==NIS_MODE_FILE) then
        	close(unit=NIS_IO_UNIT,iostat=ierr)
        end if
    end subroutine nis_close_file

    !
    ! Rewinds the position of the stream pointers
    !
    subroutine nis_rewind()
        p_current => stream_types
        i_current => stream_integers
        r_current => stream_reals
        c_current => stream_strings
    end subroutine nis_rewind

    !
    ! Returns .true. if the nis has at least one more piece of data
    !
    function nis_has_next() result (hn)
        logical :: hn

        hn = .false.
        if (associated(p_current%next)) then
            hn = .true.
        end if

    end function nis_has_next

    !
    ! Returns current nis data type
    !
    function nis_get_type() result (tp)
        character :: tp
        tp = p_current%data
    end function nis_get_type

    !
    ! Returns current nis data size
    !
    function nis_get_size() result (sz)
        integer, dimension(2) :: sz
        if (p_current%data == "i") then
            sz(1) = size(i_current%data,1)
            sz(2) = size(i_current%data,2)
        else if (p_current%data == "r") then
            sz(1) = size(r_current%data,1)
            sz(2) = size(r_current%data,2)
        else if (p_current%data == "c") then
            sz(1) = size(c_current%data,1)
            sz(2) = 1 !size(c_current%data,2)
        end if
    end function nis_get_size

    !
    ! Positions all pointers to the next piece of data
    !
    subroutine nis_next(err)
        integer, intent(out) :: err
        err = 0
        if (.not.associated(p_current%next)) then
            err = -1
            return
        end if
        p_current => p_current%next
        if (p_current%data == "i") then
            if (.not.associated(i_current%next)) then
                err = -2
                return
            end if
            i_current => i_current%next
        else if (p_current%data == "r") then
            if (.not.associated(r_current%next)) then
                err = -3
                return
            end if
            r_current => r_current%next
        else if (p_current%data == "c") then
            if (.not.associated(c_current%next)) then
                err = -4
                return
            end if
            c_current => c_current%next
        end if
    end subroutine nis_next


    !
    ! Add integer scalar
    !
    subroutine add_new_data_i11(datum)
        integer, intent(in) :: datum
        type(char_node), pointer     :: new_type
        type(integer_node), pointer  :: new_integer

        !
        ! Type
        !

        ! create new member of the chain
        allocate(new_type)
        new_type%index = current_type%index + 1
        new_type%data  = 'i'
        nullify(new_type%next)

        ! point current to the new_type
        current_type%next => new_type

        ! new current is the new_type
        current_type => new_type

        !
        ! Data
        !

        ! create new member of the chain
        allocate(new_integer)
        new_integer%index = current_integer%index + 1
        allocate(new_integer%data(1,1))
        new_integer%data = datum
        nullify(new_integer%next)

        current_integer%next => new_integer

        current_integer => new_integer

    end subroutine add_new_data_i11

    !
    ! Add real scalar
    !
    subroutine add_new_data_r11(datum)
        real, intent(in)             :: datum
        type(char_node), pointer     :: new_type
        type(real_node), pointer     :: new_real

        !
        ! Type
        !

        ! create new member of the chain
        allocate(new_type)
        new_type%index = current_type%index + 1
        new_type%data  = 'r'
        nullify(new_type%next)

        ! point current to the new_type
        current_type%next => new_type

        ! new current is the new_type
        current_type => new_type

        !
        ! Data
        !

        ! create new member of the chain
        allocate(new_real)
        new_real%index = current_real%index + 1
        allocate(new_real%data(1,1))
        new_real%data = datum
        nullify(new_real%next)

        current_real%next => new_real

        current_real => new_real

    end subroutine add_new_data_r11

    !
    ! Add real matrix
    !
    subroutine add_new_data_rmn(datum)
        real, dimension(:,:) , intent(in)  :: datum
        type(char_node), pointer     :: new_type
        type(real_node), pointer     :: new_real

        !
        ! Type
        !

        ! create new member of the chain
        allocate(new_type)
        new_type%index = current_type%index + 1
        new_type%data  = 'r'
        nullify(new_type%next)

        ! point current to the new_type
        current_type%next => new_type

        ! new current is the new_type
        current_type => new_type

        !
        ! Data
        !

        ! create new member of the chain
        allocate(new_real)
        new_real%index = current_real%index + 1
        allocate(new_real%data(size(datum,1),size(datum,2)))
        new_real%data = datum
        nullify(new_real%next)

        current_real%next => new_real

        current_real => new_real

    end subroutine add_new_data_rmn

    !
    ! Add integer matrix
    !
    subroutine add_new_data_imn(datum)
        integer, dimension(:,:) , intent(in)  :: datum
        type(char_node), pointer     :: new_type
        type(integer_node), pointer     :: new_integer

        !
        ! Type
        !

        ! create new member of the chain
        allocate(new_type)
        new_type%index = current_type%index + 1
        new_type%data  = 'i'
        nullify(new_type%next)

        ! point current to the new_type
        current_type%next => new_type

        ! new current is the new_type
        current_type => new_type

        !
        ! Data
        !

        ! create new member of the chain
        allocate(new_integer)
        new_integer%index = current_integer%index + 1
        allocate(new_integer%data(size(datum,1),size(datum,2)))
        new_integer%data = datum
        nullify(new_integer%next)

        current_integer%next => new_integer

        current_integer => new_integer

    end subroutine add_new_data_imn

    !
    ! Add character matrix
    !
    subroutine add_new_data_string(datum)
        character, dimension(:) , intent(in)  :: datum
        type(char_node), pointer       :: new_type
        type(string_node), pointer     :: new_string

        !
        ! Type
        !

        ! create new member of the chain
        allocate(new_type)
        new_type%index = current_type%index + 1
        new_type%data  = 'c'
        nullify(new_type%next)

        ! point current to the new_type
        current_type%next => new_type

        ! new current is the new_type
        current_type => new_type

        !
        ! Data
        !

        ! create new member of the chain
        allocate(new_string)
        new_string%index = current_string%index + 1
        allocate(new_string%data(size(datum,1))) !,size(datum,2)))
        new_string%data = datum
        nullify(new_string%next)

        current_string%next => new_string

        current_string => new_string

    end subroutine add_new_data_string

    !
    ! Returns current real
    !
    function get_real_m(sz) result (rr)
        integer, intent(in), dimension(2) :: sz
        real, dimension(sz(1),sz(2)) :: rr

        rr = r_current%data

    end function get_real_m

    function get_real_v(sz) result (rr)
        integer, intent(in) :: sz
        real, dimension(sz) :: rr

        rr = r_current%data(1:sz,1)

    end function get_real_v

    function get_real_s() result (rr)
        real(dp):: rr

        rr = r_current%data(1,1)

    end function get_real_s

    !
    ! Returns current integer
    !
    function get_integer_m(sz) result (rr)
        integer, intent(in), dimension(:)    :: sz
        integer(i4b), dimension(sz(1),sz(2)) :: rr

        rr = i_current%data

    end function get_integer_m

    !
    ! Returns current integer
    !
    function get_integer_s() result (rr)
        integer :: rr

        rr = i_current%data(1,1)

    end function get_integer_s

    function get_integer_v(sz) result (rr)
        integer, intent(in) :: sz
        integer, dimension(sz) :: rr

        rr = i_current%data(1:sz,1)

    end function get_integer_v

    !
    ! Returns current string
    !
    function nis_get_string(sz) result (rr)
        integer, intent(in), dimension(:)    :: sz
        character, dimension(sz(1)) :: rr

        if (sz(2) > 1) then
            stop
        else
            rr = c_current%data
        end if

    end function nis_get_string


    !
    ! Send the NIS on screen or to other processes
    !
    subroutine nis_send
        integer :: err
        integer, dimension(2) :: sz
        integer, dimension(1,1) :: int_buff
        character(len=64)  :: form, slen
        integer :: len

        call nis_rewind()
        call nis_send_char("START_NIS_VER_1.0 ")

        do

            if (.not.nis_has_next()) then
                call nis_send_char("e")
                exit
            end if

            call nis_next(err)
            if (err < 0) then
                print *, "NIS Error: ", err
                stop
            end if

            if (nis_get_type() == "i") then
                call nis_send_char("m")
                sz = nis_get_size()
                call nis_send_int(sz(1))
                call nis_send_int(sz(2))
                call nis_send_char("i")
                call nis_send_int(nis_get_integer(sz))
            else if (nis_get_type() == "r") then
                call nis_send_char("m")
                sz = nis_get_size()
                call nis_send_int(sz(1))
                call nis_send_int(sz(2))
                call nis_send_char("r")
                call nis_send_real(nis_get_real(sz))
            else if (nis_get_type() == "c") then
                call nis_send_char("v")
                sz = nis_get_size()
                call nis_send_int(sz(1))
                call nis_send_char("c")
                write(slen,*) adjustl(nis_get_string(sz))
                call nis_send_char(trim(slen))
            end if

        end do

    end subroutine nis_send

!
! DISTRIBUTE
!

    !
    ! Sends the data from nis to all processes
    !
    subroutine nis_distribute
        integer :: err
        integer, dimension(2) :: sz
        integer, dimension(1,1) :: int_buff
        character(len=24)  :: form, slen
        integer :: len

        if (PARALLEL_USING_MPI) then

            call nis_set_mode(NIS_MODE_MPI)

            if (parallel_get_my_id()==0) then

                ! send
                call nis_send()

            else

                ! receive
                call nis_read()

            end if

        end if

    end subroutine nis_distribute


    !
    ! Distributes the seeds among all processes
    ! depending on their number, number of runs each of them
    ! is supposed to pass, and the estimated number of
    ! seeds needed per run
    !
    subroutine nis_distribute_seeds()


        integer(i4b) :: step, runs_per_proc, ssize
        integer :: i

        integer(i4b), dimension(:,:), allocatable :: seeds

        character(len=256) :: cbuff
        character(len=16) :: abuff

		real(dp), dimension(10) :: aaa

        call init_random()


        runs_per_proc = parallel_nr_runs_total/parallel_nr
        step = int(real(parallel_nr_uspr*runs_per_proc)*1.1)  ! step will be by 10% larger than the estimate

        ! size of the seed
        ssize = random_get_seedsize()
        allocate(rseed(ssize))

        call nis_set_mode(NIS_MODE_MPI)

        if (parallel_id == 0) then

            allocate(seeds(parallel_nr,ssize))

            ! prepare seeds
            do i = 1, parallel_nr
                call random_get_seed(seeds(i,:))
                call random_make_steps(step)
            end do

            ! distribute them
            do i = 1, ssize
                !print *, "spreading ", i
                call parallel_spread_int(seeds(:,i),rseed(i))
            end do

            deallocate(seeds)

        else

            do i = 1, ssize
                !print *, "catching ", i
                call parallel_catch_int(rseed(i))
            end do

        end if

        ! seeds in this process

        !write(abuff,'(i2)') ssize
        !write(cbuff,"(a,i3,a,"//trim(abuff)//"i12)") "seed of ",parallel_id," is ",rseed
        !call print_log_message(trim(cbuff),5)


    end subroutine nis_distribute_seeds


!
! READING
!

    !
    ! Reading input stream from the user layer or from the root
    !
    subroutine nis_read()

        character(len=36) :: char_buffer
        character(len=1)  :: char
        integer :: size1, size2

        ! local
        real, dimension(:,:), allocatable :: rmatrix
        integer, dimension(:,:), allocatable :: imatrix
        character(len=256) :: charb2
        character, dimension(:), allocatable :: cmatrix
        complex(dpc), dimension(:,:), allocatable :: zmatrix
        integer :: i

        integer :: isize1,isize2,rsize1,rsize2,csize1,csize2,zsize1,zsize2


		if (parallel_id == 0) then

			call nis_open_file(69,"nis.in") ! this does not really open it
			nis_mode = NIS_MODE_FILE
            open(unit=NIS_IO_UNIT,file=trim(output_file))

		end if

        if (((nis_mode==NIS_MODE_FILE).and.(parallel_id == 0)) &
           .or. (nis_mode==NIS_MODE_MPI)) then

            call nis_read_char(char_buffer)

			! First check if there is no switch to reading from file
			if (char_buffer == "SWITCH_MODE_TO_FILE") then

				call nis_read_char(char_buffer)
				write(*,*) char_buffer

				call nis_open_file(69,char_buffer) ! this does not really open it
				nis_mode = NIS_MODE_FILE
                open(unit=NIS_IO_UNIT,file=trim(output_file))


	            call nis_read_char(char_buffer)

			end if

            if (char_buffer == "START_NIS_VER_1.0") then

                ! default sizes
                isize1 = 10
                isize2 = 10
                rsize1 = 100
                rsize2 = 100
                zsize1 = 100
                zsize2 = 100
                csize1 = 64
                csize2 = 1 !64

                allocate(rmatrix(rsize1,rsize2),imatrix(isize1,isize2),cmatrix(csize1)) !,csize2))
                allocate(zmatrix(zsize1,zsize2))
                do

                    ! shape
                    call nis_read_char(char)

                    if (char == "e") then

                    	call nis_close_file()
                    	exit

					end if

                    select case (char)

                        case ("m")   ! matrix

                            call nis_read_int(size1)
                            call nis_read_int(size2)

                        case ("v")   ! vector

                            call nis_read_int(size1)
                            size2 = 1

                        case ("s")   ! scalar

                            size1 = 1
                            size2 = 1

                        case default

                            print '(a,a,a)', "NIS Input Error: unknown data shape -->", char, "<--"
                            stop

                    end select

                    ! type

                    call nis_read_char(char)

                    select case (char)

                        case ("r")   ! real

                            if ((rsize1 < size1).or.(rsize2 < size2)) then
                                deallocate(rmatrix)
                                allocate(rmatrix(size1,size2))
                                rsize1 = size1
                                rsize2 = size2
                            end if

                            call nis_read_real(rmatrix(1:size1,1:size2))
                            call nis_add(rmatrix(1:size1,1:size2))

                        case ("i")   ! integer

                            if ((isize1 < size1).or.(isize2 < size2)) then
                                deallocate(imatrix)
                                allocate(imatrix(size1,size2))
                                isize1 = size1
                                isize2 = size2
                            end if

                            call nis_read_int(imatrix(1:size1,1:size2))
                            call nis_add(imatrix(1:size1,1:size2))

                        case ("c")   ! character

                            if (size2 > 1) then
                                stop
                            end if

                            if ((csize1 < size1)) then !.or.(csize2 < size2)) then
                                deallocate(cmatrix)
                                allocate(cmatrix(size1)) !,size2))
                                csize1 = size1
                            end if

                            call nis_read_char(charb2)
                            charb2 = adjustl(charb2)
                            do i = 1, size1
                                cmatrix(i) = charb2(i:i)
                            end do
                            call nis_add(cmatrix(1:size1)) !,1:size2))

                        case ("z")   ! complex


                        case default

                            print *, "NIS Input Error: unknown data type"
                            stop

                    end select

                end do


                deallocate(imatrix,cmatrix,rmatrix)

            else

                print *, "NIS input error: wrong initial string!"
                stop

            end if

        end if


    end subroutine nis_read


!
!  READING
!

    !
    ! Reading character type from stdio or MPI message
    ! depending on the setting
    !
    subroutine nis_read_char(char)
        character(len=*), intent(out) :: char

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_recieve_char(char)

        else if (nis_mode == NIS_MODE_STD) then

            read(*,'(a)') char

        else if (nis_mode == NIS_MODE_FILE) then

        	read(NIS_IO_UNIT,'(a)') char

        end if

    end subroutine nis_read_char

    !
    ! Reads an integer matrix from NIS
    !
    subroutine nis_read_int_nm(int)
        integer, intent(out), dimension(:,:) :: int

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_recieve_integer(int)

        else if (nis_mode == NIS_MODE_STD) then

            read(*,*) int

        else if (nis_mode == NIS_MODE_FILE) then

        	read(NIS_IO_UNIT,*) int

        end if

    end subroutine nis_read_int_nm

   !
   ! Reads an integer scalar from NIS
   !
   subroutine nis_read_int_11(int)
        integer, intent(out) :: int

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_recieve_integer(int)

        else if (nis_mode == NIS_MODE_STD) then

            read(*,*) int

        else if (nis_mode == NIS_MODE_FILE) then

        	read(NIS_IO_UNIT,*) int

        end if

    end subroutine nis_read_int_11

    !
    ! Reads a real scalar from NIS
    !
    subroutine nis_read_real_11(rr)
        real, intent(out) :: rr

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_recieve_real(rr)

        else if (nis_mode == NIS_MODE_STD) then

            read(*,*) rr

        else if (nis_mode == NIS_MODE_FILE) then

        	read(NIS_IO_UNIT,*) rr

        end if

    end subroutine nis_read_real_11

    !
    ! Reads a real matrix from NIS
    !
    subroutine nis_read_real_nm(rr)
        real, intent(out), dimension(:,:) :: rr

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_recieve_real(rr)

        else if (nis_mode == NIS_MODE_STD) then

            read(*,*) rr

        else if (nis_mode == NIS_MODE_FILE) then

        	read(NIS_IO_UNIT,*) rr

        end if

    end subroutine nis_read_real_nm

!
!   WRITING/SENDING
!

    !
    ! Sends a character string
    !
    subroutine nis_send_char(char)
        character(len=*), intent(in) :: char

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_send_char(char)

        else if (nis_mode == NIS_MODE_STD) then

            write(*,*) char

        else if (nis_mode == NIS_MODE_FILE) then

            write(NIS_IO_UNIT,'(a)') char

        end if

    end subroutine nis_send_char

    !
    !
    !
    subroutine nis_send_int_nm(int)
        integer(i4b), intent(in), dimension(:,:) :: int

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_send_integer(int)

        else if (nis_mode == NIS_MODE_FILE) then

            write(NIS_IO_UNIT,*) int

        else if (nis_mode == NIS_MODE_STD) then

            write(*,*) int

        end if

    end subroutine nis_send_int_nm

   !
   !
   !
   subroutine nis_send_int_11(int)
        integer(i4b), intent(in) :: int

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_send_integer(int)

        else if (nis_mode == NIS_MODE_FILE) then

            write(NIS_IO_UNIT,*) int

        else if (nis_mode == NIS_MODE_STD) then

            write(*,*) int

        end if

    end subroutine nis_send_int_11

    !
    ! Sends real scalar
    !
    subroutine nis_send_real_11(rr)
        real, intent(in) :: rr

        if (nis_mode == NIS_MODE_MPI) then


        else if (nis_mode == NIS_MODE_FILE) then

            write(NIS_IO_UNIT,*) rr

        else if (nis_mode == NIS_MODE_STD) then


            write(*,*) rr

        end if

    end subroutine nis_send_real_11

    !
    ! Sends real matrix
    !
    subroutine nis_send_real_nm(rr)
        real, intent(in), dimension(:,:) :: rr

        if (nis_mode == NIS_MODE_MPI) then

            call parallel_send_real(rr)

        else if (nis_mode == NIS_MODE_FILE) then

            write(NIS_IO_UNIT,*) rr

        else if (nis_mode == NIS_MODE_STD) then


            write(*,*) rr

        end if

    end subroutine nis_send_real_nm


!
! SETTINGS
!

    !
    ! Sets the NIS mode
    !
    !   MPI ... reads and writes to and from other processes
    !   STD ... reads and writes to and from standard i/o
    !

    subroutine nis_set_mode(mode)
        character(len=*), intent(in) :: mode

        if (mode == NIS_MODE_MPI) then
            nis_mode = NIS_MODE_MPI
            call nis_close_file()
        else if (mode == NIS_MODE_STD) then
            nis_mode = NIS_MODE_STD
            call nis_close_file()
        else if (mode == NIS_MODE_FILE) then
            nis_mode = NIS_MODE_FILE
            if (output_file_set) then
                open(unit=NIS_IO_UNIT,file=trim(file_join(out_dir,output_file)))
            else
                write(cbuff,'(a,i1,a)') NIS_DEFAULT_FILENAME//"_",parallel_id,".log"
                write(cbuff,'(a,a)') trim(file_join(out_dir,trim(cbuff)))
                open(unit=NIS_IO_UNIT,file=trim(cbuff))
            end if
        else
            write(*,*) "Error: unknown NIS mode"
        end if

    end subroutine nis_set_mode


    subroutine clean_nis
		integer :: err
		type(integer_node), pointer :: ci
		type(real_node), pointer :: cr
		type(string_node), pointer :: cs
		call nis_rewind()

		err = 1

        if (1 == 2) then
		do

		if (p_current%data == "i") then

			ci => i_current
			if (nis_has_next()) then
			  	call nis_next(err)
			else
				err = -1
			end if
			deallocate(ci%data)
			deallocate(ci)
			print *, "deallocating i"

	    else if (p_current%data == "r") then

			cr => r_current
			if (nis_has_next()) then
			  	call nis_next(err)
			else
				err = -1
			end if
			deallocate(cr%data)
			deallocate(cr)
			print *, "deallocating r"

	    else if (p_current%data == "s") then

			cs => c_current
			if (nis_has_next()) then
			  	call nis_next(err)
			else
				err = -1
			end if
			deallocate(cs%data)
			deallocate(cs)
			print *, "deallocating s"

	    else

	    	call nis_next(err)

		end if

		if (err < 0) exit

		end do

        end if

    end subroutine clean_nis

end module nis
