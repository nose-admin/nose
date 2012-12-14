#include "util_allocation.h"

!
! General Resources available to all entities in NOSE program
!
module resources

    use std_types
    use util_allocation

    use nis


    implicit none

    !
    ! Public resources
    !

	!
	! Constants
	!

	character(len=50), parameter :: NOSE_RDM_B01 = "dens_opt_coh" 	! "coherences"				!
	character(len=50), parameter :: NOSE_RDM_B11 = "dens_exc_block" 	! "populations_coherences"	!
	character(len=50), parameter :: NOSE_RDM_B12 = "dens_2_coh"		! "1to2_ex_coh"				!

	character(len=50), parameter :: NOSE_RDM_B01_ABS = "abs_dens_opt_coh" 	! ""
	character(len=50), parameter :: NOSE_RDM_B11_ABS = "abs_dens_exc_block" 	! ""
	character(len=50), parameter :: NOSE_RDM_B12_ABS = "abs_dens_2_coh" 		! ""

	character(len=50), parameter :: NOSE_RDM_B01_CONJG = "conjg_dens_opt_coh" 	! ""
	character(len=50), parameter :: NOSE_RDM_B11_CONJG = "conjg_dens_exc_block" 	! ""
	character(len=50), parameter :: NOSE_RDM_B12_CONJG = "conjg_dens_2_coh" 		! ""


	! regularization constant / Tmax
	double precision, parameter	:: NOSE_QME_NZ_REG_TIME_PER_TMAX = 0.250_dp


    !
    ! Types
    !


    !
    ! Block in the site basis (as described in the input file)
    !
    type site_block
        integer                                     :: id

        real(dp), dimension(:), pointer          :: en => NULL()			! site energies
        real(dp), dimension(:), pointer          :: en_orig => NULL()		! site energies without disorder
        integer , dimension(:), pointer          :: gindex => NULL()
        integer,  dimension(:), pointer          :: dindex => NULL()
        real(dp), dimension(:,:), pointer          :: dwidth => NULL()
        real(dp), dimension(:), pointer          :: ll => NULL()
        real(dp), dimension(:,:), pointer        :: dx => NULL()
        real(dp), dimension(:,:), pointer        :: dy => NULL()
        real(dp), dimension(:,:), pointer        :: dz => NULL()
        real(dp), dimension(:,:), pointer        :: dd => NULL()  	! transition dipole moment length
        real(dp), dimension(:), pointer          :: rs => NULL()  	! rotational strength
        real(dp), dimension(:,:), pointer        :: rr => NULL()  	! positions / or relaxation rates in case of inter block
        real(dp), dimension(:,:), pointer        :: J => NULL()   	! couplings between sites
        real(dp), dimension(:,:,:), pointer 		:: Jcell => NULL()	! coupling to neighbouring cells
        real(dp), dimension(:), pointer			:: cell_vec1 => NULL()
        real(dp), dimension(:), pointer			:: cell_vec2 => NULL()
        real(dp), dimension(:), pointer			:: cell_vec3 => NULL()
        integer, pointer                          :: N1 => NULL()
		integer										:: periodic
		integer										:: dimension
		integer										:: periodicity
		integer(i4b), dimension(:), pointer		:: QHO_lvls => NULL()	! number of levels
		real(dp), dimension(:), pointer			:: QHO_freq => NULL()	! fundamental transition energy
		real(dp), dimension(:), pointer			:: QHO_hrfact => NULL()

        integer                                     :: index
        type(site_block), pointer                  :: next => NULL()

    end type site_block

    !
    ! Block in the excition basis
    !
    type exciton_block
        integer                                  :: id

        ! single exicton part
        integer, pointer                          :: N1 => NULL()		! number of states
        real(dp), dimension(:), pointer          :: eng => NULL()		! ground state energies
        real(dp), dimension(:), pointer          :: en => NULL()		! excitonic energies
        real(dp), dimension(:,:), pointer        :: dx => NULL()
        real(dp), dimension(:,:), pointer        :: dy => NULL()
        real(dp), dimension(:,:), pointer        :: dz => NULL()
        real(dp), dimension(:,:), pointer        :: dd => NULL()
        real(dp), dimension(:,:), pointer        :: SS => NULL()
        real(dp), dimension(:,:), pointer        :: S1 => NULL()
        real(dp), dimension(:), pointer          :: rm => NULL()		! rotational strengths
        real(dp), dimension(:,:), pointer        :: rr => NULL()		! relaxation rates
        real(dp), dimension(:,:), pointer        :: dr => NULL()		! dephasing rates
        real(dp), dimension(:), pointer          :: dr0 => NULL()		! dephasing rates
        complex(dpc), dimension(:,:), pointer    :: gg => NULL()		! excitonic g(t)
        complex(dpc), dimension(:,:), pointer    :: cc => NULL()		! excitonic c(t)
        real(dp), dimension(:), pointer          :: ll => NULL()		! reorganization energy

        real(dp), dimension(:), pointer          :: fl_fac => NULL()	! steady state fluorescence factors

        ! related to two-excitons
        logical, pointer                           :: use_twoexcitons => NULL()  !
        integer, pointer                           :: N2 => NULL()        ! x
        real(dp), dimension(:), pointer          :: en_2 => NULL()      ! excitonic energies, x
        real(dp), dimension(:,:), pointer        :: dx_2 => NULL()      !
        real(dp), dimension(:,:), pointer        :: dy_2 => NULL()      !
        real(dp), dimension(:,:), pointer        :: dz_2 => NULL()      !
        real(dp), dimension(:,:), pointer        :: dd_2 => NULL()      !
        real(dp), dimension(:,:), pointer        :: SS_2 => NULL()      ! x
        real(dp), dimension(:,:), pointer        :: S1_2 => NULL()      !  x
        real(dp), dimension(:), pointer          :: drfe => NULL()      ! dephasing rates
        complex(dpc), dimension(:,:,:), pointer    :: gg_2 => NULL()      ! excitonic g(t)
        complex(dpc), dimension(:,:,:), pointer    :: cc_2 => NULL()      ! excitonic c(t)
        complex(dpc), dimension(:,:,:), pointer     :: gg_21 => NULL()
        real(dp), dimension(:), pointer          :: ll_2 => NULL()       ! reorganization energy
        integer, dimension(:,:), pointer         :: itwo => NULL()      ! x
        integer, dimension(:), pointer           :: ione1 => NULL()     ! x
        integer, dimension(:), pointer           :: ione2 => NULL()     ! x

        integer                                    :: index
        type(exciton_block), pointer              :: next => NULL()

    end type exciton_block




    !
    ! this is the new approach to storing blocks
    !

    ! matrix of blocks
    type site_block_matrix
        type(site_block)   , pointer :: sblock => NULL()
        type(exciton_block), pointer :: eblock => NULL()
    end type

    type(site_block_matrix), dimension(:,:), pointer :: iblocks => NULL()

    ! matrix of block evolution operators
    type block_ev_sop
        complex(dpc), dimension(:,:,:,:,:), pointer :: Ugg => NULL()
        complex(dpc), dimension(:,:,:,:,:), pointer :: Uee => NULL()
        complex(dpc), dimension(:,:,:,:,:), pointer :: Uff => NULL()
        complex(dpc), dimension(:,:,:,:,:), pointer :: Ueg => NULL()
        complex(dpc), dimension(:,:,:,:,:), pointer :: Ufe => NULL()
        complex(dpc), dimension(:,:,:,:,:), pointer :: Ufg => NULL()
        complex(dpc), dimension(:,:,:), pointer :: UfeS => NULL()	! secular version
    end type

    type(block_ev_sop), dimension(:,:), pointer :: evops => NULL()

    !
    ! Correlation function and the line shape function of a site
    !
    type site_goft
        integer                                      :: id
        integer                                      :: nr_modes
        character(len=32), dimension(:), pointer :: types => NULL()
        real(dp)                                     :: lambda  ! total reorganization energy
        real(dp), dimension(:,:), pointer         :: params => NULL()
        complex(dpc), dimension(:), pointer       :: gt => NULL()	! g(t)
        complex(dpc), dimension(:), pointer       :: ct => NULL()	! correlation function
        complex(dpc), dimension(:), pointer       :: ht => NULL()	! g dot of t
        integer                                      :: index
        type(site_goft), pointer                    :: next
    end type site_goft

    !
    ! Here we store gofts
    !
    type site_goft_vector
        type(site_goft), pointer :: goft
    end type
    type(site_goft_vector), dimension(:), pointer :: igofts


    !*****************************************************************
    ! type representing the whole resources
    !*****************************************************************
    type rsrcs
        type(site_block_matrix), pointer    :: blocks => NULL()	! blocks
        type(site_goft_vector), pointer     :: gofts => NULL()		! g(t)s
        type(block_ev_sop), pointer         :: evops => NULL()		! evolution operators

        character(len=64)                   :: modname
        real(dp)                             :: temp
        real(dp)                             :: dt
        integer(i4b), dimension(3)         :: gt
        integer(i4b), dimension(3)         :: Nt

        ! etc.

    end type rsrcs



    !*****************************************************************
    !  Configuration variables
    !*****************************************************************

    ! current module name
    character(len=64)          :: modname

    ! temperature
    real(dp)                   :: temp

    ! tau
    real(dp)                   :: tau_of_projector

    ! time step
    real(dp)                   :: dt

    ! time grid steps
    integer(i4b), dimension(3) :: gt

    ! time grid extent
    integer(i4b), dimension(3) :: Nt

    ! extended grid extent
    integer(i4b)                :: Nte

    ! polarization (output switch) and its allowed values
    character(len=64)          :: polarization
    character(len=64), dimension(3), parameter :: allowed_polarization &
        = (/"all-orders ","first-order","third-order"/)

    ! save NIS (output switch)
    character(len=12)          :: save_nis
    ! allowed values are yes_no

    ! save goft (output switch)
    character(len=12)          :: save_goft
    ! allowed values are yes_no

    ! parallel (execution switch)
    character(len=12)          :: parallel
    ! allowed values are yes_no

	character(len=12)			:: completeResultFile

    ! main input file name
    character(len=256)         :: inp_file

    ! module input file name
    character(len=256)         :: mod_inp_file

    integer                    :: restartFreq

    integer                    :: outnr, cur_outnr

    integer                    :: out_p

    character(len=32), dimension(:), allocatable :: outputs
    character(len=32), dimension(:,:), allocatable :: outppars

    ! yes and no values
    character(len=3), dimension(2), parameter :: RESOURCES_YES_NO      = (/"yes","no "/)
    character(len=*), parameter               :: RESOURCES_SITE_REP    = "SITE_REP"
    character(len=*), parameter               :: RESOURCES_EXCITON_REP = "EXCITON_REP"


    logical :: use_twoexcitons
    logical :: read_external_evops
    logical :: special_carotenoid_hack = .false.


    !****************************************************************
    ! Grid Parameters
    !****************************************************************

    ! number of basic steps
    integer                             :: grid_Nt
    real(dp), dimension(:), allocatable :: grid_t

    !****************************************************************
    ! Model Parameters
    !****************************************************************

    !
    ! Complete information about all blocks present in the input
    !
    type(site_block), pointer              :: current_s_block => NULL()
    type(site_goft), pointer               :: current_s_goft => NULL()

    ! interblock is a site block containing only coupling
    type(site_block), pointer              :: current_i_block => NULL()

    integer, pointer                        :: N1 => NULL()
    real(dp), dimension(:), pointer        :: en => NULL()
    real(dp), dimension(:,:), pointer      :: dx => NULL()
    real(dp), dimension(:,:), pointer      :: dy => NULL()
    real(dp), dimension(:,:), pointer      :: dz => NULL()
    real(dp), dimension(:,:), pointer      :: rr => NULL()
    real(dp), dimension(:,:), pointer      :: dr => NULL()
    real(dp), dimension(:), pointer        :: dr0 => NULL()
    real(dp), dimension(:,:), pointer      :: SS => NULL()
    real(dp), dimension(:,:), pointer      :: S1 => NULL()
    real(dp), dimension(:,:), pointer      :: jj => NULL()
    complex(dp), dimension(:,:), pointer  :: gg => NULL()
    complex(dp), dimension(:,:), pointer  :: cc => NULL()
    real(dp), dimension(:), pointer        :: ll => NULL()
    real(dp), dimension(:), pointer        :: rm => NULL()
    character(len=12)                       :: erepre


    real(dp), dimension(:), allocatable    :: PolA  ! laser field polarizations

    public :: Block_Sizes            ! Number of units in a block
    public :: Block_Sizes_Added      ! Sum of the block sizes up to a give index
    public :: Block_Sizes_Previous   ! Sum of the block sizes prior to a givem index
    public :: N_Excitons             ! Number of excitons in a block
    public :: N_Excitons_Added       ! Sum of the excitons up to a give index
    public :: N_Excitons_Previous    ! Sum of the excitons prior to a given index
    public :: Exciton_Order          ! Order of excitonic theory (=2)
    public :: N_ExOrder
    public :: N_ExOrder_Added
    public :: N_ExOrder_Previous

    !public :: bndx

    integer(i4b)                                   :: Exciton_Order
    integer(i4b), dimension(:), allocatable     :: Block_Sizes
    integer(i4b), dimension(:), allocatable     :: Block_Sizes_Added
    integer(i4b), dimension(:), allocatable     :: Block_Sizes_Previous
    integer(i4b), dimension(:), allocatable     :: N_Excitons
    integer(i4b), dimension(:), allocatable     :: N_Excitons_Added
    integer(i4b), dimension(:), allocatable     :: N_Excitons_Previous
    integer(i4b), dimension(:,:), allocatable   :: N_ExOrder
    integer(i4b), dimension(:,:), allocatable   :: N_ExOrder_Added
    integer(i4b), dimension(:,:), allocatable   :: N_ExOrder_Previous

    real(dp), dimension(:,:), allocatable       :: mun
    real(dp), dimension(:,:), allocatable       :: munS
    integer(i4b), dimension(:), allocatable     :: nmax

    integer(i4b), dimension(:), allocatable     :: NN


    !****************************************************************
    ! Derived Parameters
    !****************************************************************

    !
    ! Complete information about all blocks in eigen representation
    !
    type(exciton_block), pointer         :: current_e_block => NULL()



    integer                              :: nr_blocks  ! number of available blocks
    integer                              :: nr_gofts   ! number of available g(t)s
    integer                              :: nr_i_blocks ! number of interblocks

    real(dp)                             :: rwa        ! Rotating wave approximation frequency


    integer                              :: N_realizations  ! how many realizations will be used for averaging
    integer                              :: N_realizations_local
    integer, private                    :: N_rest
    integer                              :: i_realization
    integer								  :: N_since_saved

    !****************************************************************
    ! Management variables
    !****************************************************************
    type(site_block), pointer, private    :: root_s_block => NULL()
    type(site_block), pointer, private    :: root_i_block => NULL()
    type(exciton_block), pointer, private :: root_e_block => NULL()
    type(site_goft), pointer, private     :: root_s_goft => NULL()
    integer, parameter,private            :: MAX_NR_GOFT_MODES  = 20
    integer, parameter,private            :: MAX_NR_MODE_PARAMS = 4

	real(dp), dimension(:,:,:), allocatable :: Bands_1ex

	real(dp) :: ggC


	integer, parameter :: INI_ULTRA_FAST_PULSE = 1
	integer, parameter :: INI_POPULATION_EXCITON = 2
	integer, parameter :: INI_POPULATION_SITES = 3
	integer, parameter :: INI_DENSITY_MATRIX = 4

	integer :: ini_basis_type

	! initial density matrix for dynamics calculations
	complex(dpc), dimension(:,:), allocatable :: ini_dm

    !
    ! private
    !


contains

    !
    ! Initializes all resources
    !
    subroutine init_resources()

        ! initialize site blocks
        allocate(root_s_block)
        root_s_block%index = 0
        root_s_block%next  => null()

        allocate(root_i_block)
        root_i_block%index = 0
        root_i_block%next  => null()

        ! initialize excion blocks
        allocate(root_e_block)
        root_e_block%index = 0
        root_e_block%next  => null()

        allocate(root_s_goft)
        root_s_goft%index = 0
        root_s_goft%gt => null()
        root_s_goft%next => null()

        nr_blocks = 0
        nr_gofts  = 0
        nr_i_blocks = 0

        call resources_rewind_blocks()
        call resources_rewind_gofts()
        call resources_rewind_inter_blocks()

        ! set all pointers to null
        en => null()


        use_twoexcitons = .true. !.false.

        completeResultFile = 'no'

        ini_basis_type = INI_ULTRA_FAST_PULSE

    end subroutine init_resources

    !
    ! Sets resources that are known only after input has been read
    !
    subroutine resources_after_input()

        integer :: N_goft
        integer :: i, nma
        real(dp), dimension(:), allocatable    :: r
        character(len=64) :: cbuff


        !
        ! Handle realizations distribution
        !
        N_realizations_local = N_realizations/parallel_nr
        N_rest = N_realizations - parallel_nr*N_realizations_local

        !
        ! Message about the parallel task distribution
        !
        write(cbuff,'(a,i5,a)') "Total of ",N_realizations_local," on each processor"
        call print_log_message(trim(cbuff),5)
        write(cbuff,'(a,i5,a,i5)') "procesore inbalance: ", N_rest, " on ", &
            parallel_nr
        call print_log_message(trim(cbuff),5)
        if ((parallel_id > 0).and.(parallel_id <= N_rest)) then
            N_realizations_local = N_realizations_local + 1
        end if


        !
        ! Set the number of sites in each block
        !
        ALLOCATE(NN,(nr_blocks))

        call resources_rewind_blocks()
        do

            NN(current_s_block%id) = current_s_block%N1

            if (.not.resources_have_next_block()) exit
            call resources_next_block()

        end do

        !
        ! Number of different correlation function definitions
        !
        N_goft        = nr_gofts

        !
        ! We use at maximum the two-exciton states
        !
        Exciton_Order = 2


        !
        ! Initialization of usefull constants
        !
        ALLOCATE(Block_Sizes,(nr_blocks))
        ALLOCATE(Block_Sizes_Added,(nr_blocks))
        ALLOCATE(Block_Sizes_Previous,(nr_blocks))
        ALLOCATE(N_Excitons,(nr_blocks))
        ALLOCATE(N_Excitons_Added,(nr_blocks))
        ALLOCATE(N_Excitons_Previous,(nr_blocks))
        ALLOCATE(N_ExOrder,(nr_blocks,Exciton_Order))
        ALLOCATE(N_ExOrder_Added,(nr_blocks,Exciton_Order))
        ALLOCATE(N_ExOrder_Previous,(nr_blocks,Exciton_Order))
        do i = 1, nr_blocks
            Block_Sizes(i) = NN(i)
            N_ExOrder(i,1) = NN(i)
            N_ExOrder(i,2) = NN(i)*(NN(i)-1)/2
            N_Excitons(i)  = N_ExOrder(i,1) + N_ExOrder(i,2)
            if (i==1) then
                Block_Sizes_Added(i) = NN(i)
                Block_Sizes_Previous(i) = 0
                N_Excitons_Added(i) = NN(i) + NN(i)*(NN(i)-1)/2
                N_Excitons_Previous(i) = 0
                N_ExOrder_Added(i,1) = NN(i)
                N_ExOrder_Added(i,2) = NN(i)*(NN(i)-1)/2
                N_ExOrder_Previous(i,1) = 0
                N_ExOrder_Previous(i,2) = 0
            else
                Block_Sizes_Added(i) = Block_Sizes_Added(i-1) + NN(i)
                Block_Sizes_Previous(i) = Block_Sizes_Added(i-1)
                N_Excitons_Added(i) = N_Excitons_Added(i-1) + &
                    NN(i) + NN(i)*(NN(i)-1)/2
                N_Excitons_Previous(i) = N_Excitons_Added(i-1)
                N_ExOrder_Added(i,1) = N_ExOrder_Added(i-1,1) + NN(i)
                N_ExOrder_Added(i,2) = N_ExOrder_Added(i-1,2) + &
                    NN(i)*(NN(i)-1)/2
                N_ExOrder_Previous(i,1) = N_ExOrder_Added(i-1,1)
                N_ExOrder_Previous(i,2) = N_ExOrder_Added(i-1,2)
            end if
        end do

        ALLOCATE(nmax,(nr_blocks))
        nma = 0
        do i=1,nr_blocks
            nma=nma+N_ExOrder_Previous(i,1)*(N_ExOrder_Previous(i,2)+1)
            nmax(i)=nma
        end do

        if (.not.allocated(mun)) then
           ALLOCATE(mun,(nmax(nr_blocks)+N_ExOrder(nr_blocks,1)*(N_ExOrder(nr_blocks,2)+1),3))
        end if
        if (.not.allocated(munS)) then
           ALLOCATE(munS,(nmax(nr_blocks)+N_ExOrder(nr_blocks,1)*(N_ExOrder(nr_blocks,2)+1),3))
        end if

        !
        ! Set polarizations
        !

        allocate(PolA(4))
        PolA = 0.0_dp

        !
        ! Calculate and/or set RWA
        !


    end subroutine resources_after_input


!******************************************************************************
!    Input file blocks
!******************************************************************************

    !
    ! Creates a new empty s and e blocks and appends them to the chains
    ! Points current_e_block and current_s_block to them
    !
    subroutine resources_new_block()
        if ((current_s_block%index == 0)) then
            ! new block is the root block itself
            current_s_block%index = 1
            current_s_block%next => null()
        else
            ! completely new block
            allocate(current_s_block%next)
            current_s_block%next%index = current_s_block%index + 1
            current_s_block => current_s_block%next
            current_s_block%next => null()
        end if
        nr_blocks = nr_blocks + 1
    end subroutine resources_new_block

    !
    ! Create a new exciton block corresponding to the current site_block
    !
    subroutine resources_new_ex_block()

        if ((root_e_block%index == 0)) then
            current_e_block => root_e_block
            current_e_block%index = 1
            current_e_block%next => null()
        else
            allocate(current_e_block%next)
            current_e_block%next%index = current_e_block%index + 1
            current_e_block => current_e_block%next
            current_e_block%next => null()
        end if

        current_e_block%id = current_s_block%id
        allocate(current_e_block%N1)
        current_e_block%N1 = current_s_block%N1

    end subroutine resources_new_ex_block


    !
    ! Sets the pointer to the start of the block chain
    !
    subroutine resources_rewind_blocks()
        current_s_block => root_s_block
        current_e_block => root_e_block
        call resources_set_all_pointers(RESOURCES_SITE_REP)
    end subroutine resources_rewind_blocks

    !
    ! Returns true if the resources have next block
    !
    function resources_have_next_block() result(have)
        logical :: have
        have = .false.
        if (associated(current_s_block%next)) then
            have = .true.
        end if
    end function resources_have_next_block

    !
    ! Makes one step along to chain of gofts
    !
    subroutine resources_next_block()
        current_s_block => current_s_block%next
        if (associated(current_e_block%next)) then
            current_e_block => current_e_block%next
        end if
        call resources_set_all_pointers(RESOURCES_SITE_REP)
    end subroutine resources_next_block


!************************************************************************
!  Inter blocks
!************************************************************************

    !
    ! Creates a new empty s and e blocks and appends them to the chains
    ! Points current_e_block and current_s_block to them
    !
    subroutine resources_new_inter_block()
        if ((current_i_block%index == 0)) then
            ! new block is the root block itself
            current_i_block%index = 1
            current_i_block%next => null()
        else
            ! completely new block
            allocate(current_i_block%next)
            current_i_block%next%index = current_i_block%index + 1
            current_i_block => current_i_block%next
            current_i_block%next => null()
        end if
        nr_i_blocks = nr_i_blocks + 1
    end subroutine resources_new_inter_block

    !
    ! Sets the pointer to the start of the block chain
    !
    subroutine resources_rewind_inter_blocks()
        current_i_block => root_i_block
    end subroutine resources_rewind_inter_blocks


    !
    ! Returns true if the resources have next block
    !
    function resources_have_next_inter_block() result(have)
        logical :: have
        have = .false.
        if (associated(current_i_block%next)) then
            have = .true.
        end if
    end function resources_have_next_inter_block

    !
    ! Makes one step along to chain of gofts
    !
    subroutine resources_next_inter_block()
        current_i_block => current_i_block%next
    end subroutine resources_next_inter_block



!************************************************************************
!   Input file gofts
!************************************************************************

    !
    ! Sets the pointer to the start of the goft chain
    !
    subroutine resources_rewind_gofts()
        current_s_goft => root_s_goft
    end subroutine resources_rewind_gofts

    !
    ! Creates a new empty goft and appends it to the chain
    ! Points current_s_goft to it
    !
    subroutine resources_new_goft()
        if (current_s_goft%index == 0) then
            ! new goft is the root goft itself
            current_s_goft%index = 1

        else
            ! completely new block
            allocate(current_s_goft%next)
            current_s_goft%next%index = current_s_goft%index + 1
            current_s_goft => current_s_goft%next
            current_s_goft%next => null()
        end if
        current_s_goft%nr_modes = 0
!        call allocate_pointer(current_s_goft%params,MAX_NR_MODE_PARAMS,MAX_NR_GOFT_MODES)
        ALLOCATE(current_s_goft%params,(MAX_NR_MODE_PARAMS,MAX_NR_GOFT_MODES))
!        call allocate_pointer(current_s_goft%types,MAX_NR_GOFT_MODES)
        ALLOCATE(current_s_goft%types,(MAX_NR_GOFT_MODES))
        current_s_goft%params = 0.0_dp
        current_s_goft%types = "X"
        current_s_goft%gt => null()
        nr_gofts = nr_gofts + 1
    end subroutine resources_new_goft

    !
    ! Returns true if the resources have next goft
    !
    function resources_have_next_goft() result(have)
        logical :: have
        have = .false.
        if (associated(current_s_goft%next)) then
            have = .true.
        end if
    end function resources_have_next_goft

    !
    ! Makes one step along to chain of gofts
    !
    subroutine resources_next_goft()
        current_s_goft => current_s_goft%next
    end subroutine resources_next_goft


!*********************************************************************************
!
!*********************************************************************************

    !
    ! Sets all special pointers to point to the values in the current block
    !
    subroutine resources_set_all_pointers(repre)
        character(len=*), intent(in) :: repre
        if (trim(repre) == RESOURCES_SITE_REP) then

            en => current_s_block%en
            dx => current_s_block%dx
            dy => current_s_block%dy
            dz => current_s_block%dz
            !rr => current_s_block%rr  ! this here are the positions of the transitions
            jj => current_s_block%J
            N1 => current_s_block%N1
            ll => current_s_block%ll

            SS => null()
            S1 => null()
            gg => null()
            cc => null()
            rm => null()
            dr => null()
            dr0 => null()
            rr => null()

        else if (trim(repre) == RESOURCES_EXCITON_REP) then

            en => current_e_block%en
            dx => current_e_block%dx
            dy => current_e_block%dy
            dz => current_e_block%dz
            jj => current_s_block%J  ! still points to _s_, because there is no coupling in ex. rep
            N1 => current_e_block%N1
            gg => current_e_block%gg
            cc => current_e_block%cc
            ll => current_e_block%ll
            rm => current_e_block%rm
            SS => current_e_block%SS
            S1 => current_e_block%S1

            rr => current_e_block%rr
            dr => current_e_block%dr
            dr0 => current_e_block%dr0
            !jj => null()


        end if

        erepre = trim(repre)

    end subroutine resources_set_all_pointers

    !
    ! Checks the values agains the list of allowed values
    !
    subroutine resources_check_values(char, allowed_char,err)
        character(len=*), intent(in) :: char
        character(len=*), dimension(:), intent(in) :: allowed_char
        integer, intent(out) :: err
        ! local
        integer :: i

        err = 0

        do i = 1, size(allowed_char,1)
            if (trim(char)==trim(allowed_char(i))) return
        end do

        err = -1

    end subroutine resources_check_values

    subroutine resources_new_output(str)
        character(len=*), intent(in) :: str

        if (.not.allocated(outputs)) allocate(outputs(outnr))
        out_p = out_p + 1

        outputs(out_p) = str

    end subroutine resources_new_output

    subroutine resources_new_outppar(str)
        character(len=*), intent(in) :: str

        if (.not.allocated(outppars)) then
        	allocate(outppars(outnr,3))
        	outppars = ""
        end if
        !out_p = out_p + 1

        outppars(out_p,1) = str

    end subroutine resources_new_outppar


    !
    ! returns true if
    !
    function resources_output_contains(str) result(c)
        character(len=*), intent(in) :: str
        logical :: c
        integer :: i

        c = .false.


        i = 0
        do

            i = i + 1
            if (trim(outputs(i)) == trim(str)) then
                c = .true.
                cur_outnr = i
                exit
            end if

            if (i == outnr) then
            	cur_outnr = -1
            	exit
            end if

        end do

    end function resources_output_contains


    !
    ! Cleans all exciton blocks and resets root block
    !
    subroutine resources_clean_excitons

        type(exciton_block), pointer :: pp


        call resources_rewind_blocks()
        do

            if (associated(current_e_block%gg)) then


                DEALLOCATE(current_e_block%gg)
                DEALLOCATE(current_e_block%dx)
                DEALLOCATE(current_e_block%dd)
                DEALLOCATE(current_e_block%dy)
                DEALLOCATE(current_e_block%dz)
                !call deallocate_pointer(current_e_block%en)
                DEALLOCATE(current_e_block%en)
                !call deallocate_pointer(current_e_block%SS)
                DEALLOCATE(current_e_block%SS)
                !call deallocate_pointer(current_e_block%S1)
                DEALLOCATE(current_e_block%S1)
                deallocate(current_e_block%ll)
                deallocate(current_e_block%N1)
                !call deallocate_pointer(current_e_block%rm)
                DEALLOCATE(current_e_block%rm)
                !call deallocate_pointer(current_e_block%fl_fac)
                DEALLOCATE(current_e_block%fl_fac)
                !DEALLOCATE_S(current_e_block%N2)


            end if

            pp => current_e_block

           if (.not.resources_have_next_block()) exit
           call resources_next_block()

           deallocate(pp)

        end do

        deallocate(pp)

        ! initialize excion blocks
        allocate(root_e_block)
        root_e_block%index = 0
        root_e_block%id = 0
        root_e_block%next  => null()


    end subroutine resources_clean_excitons


    subroutine clean_resources()

		deallocate(root_s_block)

    end subroutine clean_resources

end module resources
