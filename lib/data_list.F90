module data_list

	use std_types

	implicit none
 
	private	

! public names
	public :: list_i4b, list_element_i4b
	public :: list_dp,  list_element_dp
	
	public :: list_init, list_ins_next
	public :: list_destroy
	public :: list_show
	
! types	
	type list_i4b
		integer(i4b)        :: lsize
		type(list_element_i4b), pointer :: head
		type(list_element_i4b), pointer :: tail
	end type list_i4b

	type list_element_i4b
		integer(i4b) :: value
		type(list_element_i4b), pointer :: next
	end type list_element_i4b

	type list_dp
		integer(i4b)        :: lsize
		type(list_element_dp), pointer :: head
		type(list_element_dp), pointer :: tail
	end type list_dp

	type list_element_dp
		real(dp) :: value
		type(list_element_dp), pointer :: next
	end type list_element_dp

! Constructors	
	interface list_init
		module procedure init_i4b, init_dp
	end interface
	
! Destructors
	interface list_destroy
		module procedure destroy_i4b, destroy_dp
	end interface list_destroy

! Methods
	interface list_ins_next
		module procedure ins_next_i4b, ins_next_dp
	end interface
	
	interface list_show
		module procedure show_i4b, show_dp
	end interface list_show
	
contains

!*******************************************************************************
!   Constructors
!*******************************************************************************
	subroutine init_i4b(ll)
		type(list_i4b), intent(out) :: ll
		ll%lsize = 0
		nullify(ll%head)
		nullify(ll%tail)
	end subroutine init_i4b

	subroutine init_dp(ll)
		type(list_dp), intent(out) :: ll
		ll%lsize = 0
		nullify(ll%head)
		nullify(ll%tail)
	end subroutine init_dp


!*******************************************************************************
!   Destructors
!*******************************************************************************
	subroutine destroy_i4b(ll)
		type(list_i4b), intent(inout) :: ll
		! go through all the list and delete everything
		type(list_element_i4b), pointer :: el, el_next
  		
  		if (associated(ll%head)) then
  		
  			el  => ll%head
 		
  			do 
  		
  				el_next => el%next
  				deallocate(el)		
  	
  				if (.not.associated(el_next)) exit
  		
  				el => el_next
  		
  			end do
  	
  		end if
  	
  		nullify(ll%head)
  		nullify(ll%tail)
  		ll%lsize = 0
		
	end subroutine destroy_i4b

	subroutine destroy_dp(ll)
		type(list_dp), intent(inout) :: ll
		! go through all the list and delete everything
		type(list_element_dp), pointer :: el, el_next
  		
  		if (associated(ll%head)) then
  		
  			el  => ll%head
 		
  			do 
  		
  				el_next => el%next
  				deallocate(el)		
  	
  				if (.not.associated(el_next)) exit
  		
  				el => el_next
  		
  			end do
  	
  		end if
  	
  		nullify(ll%head)
  		nullify(ll%tail)
  		ll%lsize = 0
		
	end subroutine destroy_dp

!*******************************************************************************
!   Methods
!*******************************************************************************

!******************
! insert next
!******************
! insert a new value behind a specified element. The pointer to the inserted
! element is returned.
!******************

	!
	! i4b version
	!
	subroutine ins_next_i4b(ll, lelem, value)
		type(list_i4b), intent(inout)   :: ll
		type(list_element_i4b), pointer :: lelem
		integer(i4b), intent(in)        :: value
		
		! new element
		type(list_element_i4b), pointer :: nelem
		
		! allocate new element
		allocate(nelem)
		
		! setting value of a new element 
		nelem%value = value
        
        if (ll%lsize == 0) lelem => null()        
		
		! in not associated element is sent we add to the head
		if (.not.associated(lelem)) then

			if (ll%lsize == 0) ll%tail => nelem
			nelem%next => ll%head
			ll%head    => nelem

		else
		
			! if the lelem%next is not assoc. -> it is the tail
			! we set a new tail
		 	if (.not.associated(lelem%next)) ll%tail => nelem
			
			nelem%next => lelem%next
			lelem%next => nelem
		
		end if
				
		lelem => nelem   ! we return a pointer to the last added element
		ll%lsize = ll%lsize + 1
	
	end subroutine ins_next_i4b


	!
	! dp version
	!
	subroutine ins_next_dp(ll, lelem, value)
		type(list_dp), intent(inout)   :: ll
		type(list_element_dp), pointer :: lelem
		real(dp), intent(in)           :: value
		
		! new element
		type(list_element_dp), pointer :: nelem
		
		! allocate new element
		allocate(nelem)
		
		! setting value of a new element 
		nelem%value = value
		        
        if (ll%lsize == 0) lelem => null()        
                
		! in not associated element is sent we add to the head
		if (.not.associated(lelem)) then

			if (ll%lsize == 0) ll%tail => nelem
			nelem%next => ll%head
			ll%head    => nelem

		else
		
			! if the lelem%next is not assoc. -> it is the tail
			! we set a new tail
		 	if (.not.associated(lelem%next)) ll%tail => nelem
			
			nelem%next => lelem%next
			lelem%next => nelem
		
		end if
				
		lelem => nelem   ! we return a pointer to the last added element
		ll%lsize = ll%lsize + 1
	   
        print *, ll%lsize
    
	end subroutine ins_next_dp



!******************
! show the elements of the list
!******************

	!
	! i4b version
	!
  	subroutine show_i4b(ll)
  		type(list_i4b) :: ll
  		type(list_element_i4b) :: el
  		
        if (.not.associated(ll%head)) return
            
  		print *, "list = ["
  		
  		el  = ll%head
 
  		do 
  		
  			print '(i10)', el%value
  		
  			if (.not.associated(el%next)) exit
  		
  			el = el%next
  		
  		end do
  		
  		print *, "]" 
  	
  	end subroutine show_i4b


	!
	! dp version
	!
  	subroutine show_dp(ll)
  		type(list_dp) :: ll
  		type(list_element_dp), pointer :: el
  		
        if (.not.associated(ll%head)) return
        
  		print *, "list = ["
  		
  		el  => ll%head
 
  		do 
  		
  			print '(a2,e18.12)', "  ", el%value
  		
  			if (.not.associated(el%next)) exit
  		
  			el => el%next
  		
  		end do
  		
  		print *, "]" 
  	
  	end subroutine show_dp


end module data_list
