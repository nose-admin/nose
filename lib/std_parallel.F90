!
! Support for a simple parallelism
!
module std_parallel

    use std_io

    implicit none

#ifdef HAVE_MPI

    include 'mpif.h'

    logical, public, parameter :: PARALLEL_USING_MPI = .true.

#else

    logical, public, parameter :: PARALLEL_USING_MPI = .false.

#endif

    integer :: parallel_id, parallel_nr, parallel_nodes

    interface parallel_recieve_real
        module procedure recieve_real_nm, recieve_real_11
    end interface

    interface parallel_send_real
        module procedure send_real_nm, send_real_11
    end interface

    interface parallel_send_integer
        module procedure send_integer_11 , send_integer_nm
    end interface

    interface parallel_recieve_integer
        module procedure recieve_integer_11 , recieve_integer_nm
    end interface

    private :: recieve_real_nm, recieve_real_11, send_real_11, send_real_nm, send_integer_11, recieve_integer_11
    private :: recieve_integer_nm, send_integer_nm
    character(len=256), private :: pcbuff

    integer :: parallel_nr_runs_total
    integer :: parallel_nr_uspr        ! Used Seeds Per Run


contains

    !
    ! Returns the id of current process
    !
    function parallel_get_my_id() result (id)
        integer :: id

        id = parallel_id

    end function parallel_get_my_id

    !
    ! Returns number of processors this program is run on
    !
    function parallel_get_procs() result (nr)
        integer :: nr
        nr = parallel_nr
    end function parallel_get_procs

    !
    ! Tests communication by sending simple message to
    ! the master process
    !
    subroutine parallel_communication_test
        character(len=32) :: message, str
        integer           :: dest, tag, ierr, source
#ifdef HAVE_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
#endif
        message = "Test message from "

        ierr = 0

        if (parallel_id > 0) then

            dest = 0
            tag  = 0

            write(str,'(i4)') parallel_id
            message = trim(message)//trim(str)
#ifdef HAVE_MPI
            ! sending message to 0
            call MPI_Send(message,len_trim(message),MPI_CHARACTER, &
                dest, tag, MPI_COMM_WORLD, ierr)
#endif

        else

            call print_log_message("                           Communication test ...",5)
            tag = 0

#ifdef HAVE_MPI
            ! receiving messages
            do source = 1, parallel_nr-1
                call MPI_Recv(message,32,MPI_CHARACTER,source,tag, &
                    MPI_COMM_WORLD,status,ierr)
            end do
#endif
            if (ierr == 0) then
                call print_log_message("                            ... succesful",5)
                write(pcbuff,'(a,i3,a,i3,a)') "              Working space: ", parallel_nr, " processor(s) on ", &
                                                                  parallel_nodes, " node(s)"
                call print_log_message(trim(pcbuff),5)
                call print_log_message(" ",5)
            else
                call print_error_message(ierr,"communication test return with error")
                stop
            end if

        end if

    end subroutine parallel_communication_test

!
! COLLECTING
!


    !
    ! Collect real scalar into a single array
    ! available to to the master process
    !
    subroutine parallel_collect_real(rr, ar)
        real, intent(in)   :: rr
        real, dimension(:) :: ar
        ! local
        integer           :: dest, tag, ierr, source
#ifdef HAVE_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
#endif
        real              :: r_in

        ar(1) = rr
        dest  = 0
        tag  = 0

#ifdef HAVE_MPI

        if (parallel_id == 0) then
            do source = 1, parallel_nr-1
                call MPI_Recv(r_in,1,MPI_REAL,source,tag, &
                    MPI_COMM_WORLD,status,ierr)
                ar(source+1) = r_in
            end do
        else
            call MPI_Send(rr,1,MPI_REAL, &
                dest, tag, MPI_COMM_WORLD, ierr)
        end if
#endif

    end subroutine parallel_collect_real


    !
    ! Performs parallel averaging of a given array. On each processor
    ! the array has alreadt beed summed N-times
    !
    subroutine parallel_average_real_array(rr,N)
        real, dimension(:), intent(inout) :: rr
        real, dimension(size(rr,1)) :: ra
        integer, intent(inout) :: N
        integer :: M

        integer :: count, ierr
        M  = 0
        ra = 0.0_dp

#ifdef HAVE_MPI

        count = size(rr,1)
        call MPI_Reduce(rr,ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(N,M,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (parallel_id == 0) then
             rr = ra/real(M)
             N = M
        end if

#else

        rr = rr/real(N)

#endif


    end subroutine parallel_average_real_array

    !
    ! Performs parallel averaging of a given array. On each processor
    ! the array has alreadt beed summed N-times
    !
    subroutine parallel_average_complex_array(rr,N)
        complex, dimension(:), intent(inout) :: rr
        complex, dimension(size(rr,1)) :: rc
        real, dimension(size(rr,1)) :: ra
        integer, intent(inout) :: N
        integer :: M

        integer :: count, ierr
        M  = 0
        ra = 0.0_dp
        rc = 0.0_dp

#ifdef HAVE_MPI

        count = size(rr,1)

        call MPI_Reduce(real(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = ra
        call MPI_Reduce(aimag(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = rc + (0.0_dp,1.0_dp)*ra

        call MPI_Reduce(N,M,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (parallel_id == 0) then
            rr = rc/real(M)
            N = M
        end if

#else

        rr = rr/real(N)

#endif


    end subroutine parallel_average_complex_array


    !
    ! Performs parallel averaging of a given array. On each processor
    ! the array has alreadt beed summed N-times
    !
    subroutine parallel_average_complex_array2(rr,N)
        complex, dimension(:,:), intent(inout) :: rr
        complex, dimension(size(rr,1),size(rr,2)) :: rc
        real, dimension(size(rr,1),size(rr,2)) :: ra
        integer, intent(inout) :: N

        integer :: M

        integer :: count, ierr
        M  = 0
        ra = 0.0_dp
        rc = 0.0_dp

#ifdef HAVE_MPI

        count = size(rr,1)*size(rr,2)

        call MPI_Reduce(real(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = ra
        call MPI_Reduce(aimag(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = rc + (0.0_dp,1.0_dp)*ra

        call MPI_Reduce(N,M,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (parallel_id == 0) then
            rr = rc/real(M)
            N = M
        end if

#else

        rr = rr/real(N)

#endif


    end subroutine parallel_average_complex_array2

    subroutine parallel_average_complex_array3(rr,N)
        complex, dimension(:,:,:), intent(inout) :: rr
        complex, dimension(size(rr,1),size(rr,2),size(rr,3)) :: rc
        real, dimension(size(rr,1),size(rr,2),size(rr,3)) :: ra
        integer, intent(inout) :: N

        integer :: M

        integer :: count, ierr
        M  = 0
        ra = 0.0_dp
        rc = 0.0_dp

#ifdef HAVE_MPI

        count = size(rr,1)*size(rr,2)*size(rr,3)

        call MPI_Reduce(real(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = ra
        call MPI_Reduce(aimag(rr),ra,count,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        rc = rc + (0.0_dp,1.0_dp)*ra

        call MPI_Reduce(N,M,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (parallel_id == 0) then
            rr = rc/real(M)
            N = M
        end if

#else

        rr = rr/real(N)

#endif


    end subroutine parallel_average_complex_array3


!
! SPREADING
!

    !
    ! Spreads an array from 0 to all other processes
    !
    subroutine parallel_spread_int(ai,int)
        integer, intent(in), dimension(:) :: ai
        integer, intent(out) :: int
        integer :: is, tag, ierr, dest, nr, md
#ifdef HAVE_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
#endif

#ifdef HAVE_MPI

        if (parallel_id == 0) then

            tag = 0
            nr = size(ai)/parallel_nr
            md = size(ai) - nr*parallel_nr

            do dest = 1, parallel_nr-1
                is = ai(dest+1)
                call MPI_Send(is,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr)
            end do

#endif

            int = ai(1)

#ifdef HAVE_MPI

        end if

#endif

    end subroutine parallel_spread_int

    !
    ! Catches what was spread
    !
    subroutine parallel_catch_int(int)
        integer, intent(out) :: int
        integer :: is, tag, ierr

#ifdef HAVE_MPI

        integer, dimension(MPI_STATUS_SIZE) :: status

        tag = 0

        if (parallel_id > 0) then

            call MPI_Recv(is,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
            int = is

        end if

#endif

    end subroutine parallel_catch_int


!
!  BROADCASTING
!

    !
    ! Send character to all processes
    !
    subroutine parallel_broadcast_char(char)
        character(len=*) char

#ifdef HAVE_MPI
        integer :: ierr

        call MPI_BCast(char,len(char),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

#endif

    end subroutine parallel_broadcast_char

!
!  SENDING
!

    !
    ! Sends characters from the root process
    !
    subroutine parallel_send_char(char)
        character(len=*), intent(in) :: char

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer :: dest
        integer :: length
        tag = 0

        length = len(char)
        do dest = 1, parallel_nr-1
            !print "(a,i3,a,i3,3a,i5)", "Sent from: ", parallel_id, " to ", &
            !    dest, ": ", char, ": length ",length
            call MPI_Send(length,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr)
            call MPI_Send(char,len(char),MPI_CHARACTER,dest,tag,MPI_COMM_WORLD,ierr)
        end do

#endif

    end subroutine parallel_send_char

    !
    ! Sends integers from the root process
    !
    subroutine send_integer_11(int)
        integer, intent(in) :: int

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer :: dest

        tag = 0

        do dest = 1, parallel_nr-1
            !print "(a,i3,a,i3,a,i3)", "Sent from: ", parallel_id, " to ", &
            !    dest, ": ", int
            call MPI_Send(int,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr)
        end do

#endif

    end subroutine send_integer_11

    !
    ! Sends reals from the root process
    !
    subroutine send_real_nm(rr)
        real, dimension(:,:), intent(in) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer :: dest

        tag = 0

        do dest = 1, parallel_nr-1
            !print *, "Sent from: ", parallel_id, " to ", &
            !    dest, ": ", rr
            call MPI_Send(rr,size(rr),MPI_REAL,dest,tag,MPI_COMM_WORLD,ierr)
        end do

#endif

    end subroutine send_real_nm

    !
    ! Sends reals from the root process
    !
    subroutine send_real_11(rr)
        real, intent(in) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer :: dest

        tag = 0

        do dest = 1, parallel_nr-1
            !print *, "Sent from: ", parallel_id, " to ", &
            !    dest, ": ", rr
            call MPI_Send(rr,1,MPI_REAL,dest,tag,MPI_COMM_WORLD,ierr)
        end do

#endif

    end subroutine send_real_11



    !
    ! Sends integers from the root process
    !
    subroutine send_integer_nm(rr)
        integer, dimension(:,:), intent(in) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer :: dest

        tag = 0

        do dest = 1, parallel_nr-1
            !print *, "Sent from: ", parallel_id, " to ", &
            !    dest, ": ", rr
            call MPI_Send(rr,size(rr),MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr)
        end do

#endif

    end subroutine send_integer_nm



!
!  RECEIVING
!

    !
    ! Receives character from the root process
    !
    subroutine parallel_recieve_char(char)
        character(len=*), intent(out) :: char

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer :: length
        character(len=256) :: cbuff

        tag = 0

        !print "(a,i3,a,i3)", "Receiving at ", parallel_id, " from ", 0
        call MPI_Recv(length,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Recv(cbuff(1:length),length,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,ierr)

        char = cbuff(1:length)
        !print *, "Received: ", trim(char), length

#endif

    end subroutine parallel_recieve_char


    !
    ! Receives character from the root process
    !
    subroutine recieve_integer_11(int)
        integer, intent(out) :: int

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        tag = 0

        !print "(a,i3,a,i3)", "Receiving at ", parallel_id, " from ", 0
        call MPI_Recv(int,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
        !print *, "Received: ", int

#endif

    end subroutine recieve_integer_11

    !
    ! Receives integers from the root process
    !
    subroutine recieve_integer_nm(rr)
        integer, dimension(:,:), intent(out) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        tag = 0

        !print "(a,i3,a,i3)", "Receiving at ", parallel_id, " from ", 0
        call MPI_Recv(rr,size(rr),MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
        !print *, "Received: ", rr

#endif

    end subroutine recieve_integer_nm


    !
    ! Receives character from the root process
    !
    subroutine recieve_real_11(rr)
        real, intent(out) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        tag = 0

        !print "(a,i3,a,i3)", "Receiving at ", parallel_id, " from ", 0
        call MPI_Recv(rr,1,MPI_REAL,0,tag,MPI_COMM_WORLD,status,ierr)
        !print *, "Received: ", int

#endif

    end subroutine recieve_real_11

    !
    ! Receives real from the root process
    !
    subroutine recieve_real_nm(rr)
        real, dimension(:,:), intent(out) :: rr

#ifdef HAVE_MPI

        integer :: ierr
        integer :: tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        tag = 0

        !print "(a,i3,a,i3)", "Receiving at ", parallel_id, " from ", 0
        call MPI_Recv(rr,size(rr),MPI_REAL,0,tag,MPI_COMM_WORLD,status,ierr)
        !print *, "Received: ", rr

#endif

    end subroutine recieve_real_nm






    !
    ! Initializes MPI and
    ! forks the process into multiple processes
    !
    subroutine init_parallel
        integer :: ierr, rank, size

#ifdef HAVE_MPI

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

#else

        rank = 0
        size = 1

#endif

        parallel_id = rank
        parallel_nr = size
        parallel_nodes = 1


    end subroutine init_parallel

    !
    ! Closes MPI session
    !
    subroutine finalize_parallel(err)
        integer, intent(out) :: err
        integer :: ierr

        err = 0

#ifdef HAVE_MPI
        call MPI_FINALIZE(ierr)
#endif

    end subroutine finalize_parallel


    subroutine parallel_synchronize()

#ifdef HAVE_MPI
        integer :: ierr
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

#endif

    end subroutine parallel_synchronize

end module std_parallel
