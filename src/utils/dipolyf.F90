program dipolyf

implicit none 

   real, dimension(1:3) :: r1   !poloha 1. dipolu
   real, dimension(1:3) :: r2
   real, dimension(1:3) :: d1   !1. dipolovy moment 
   real, dimension(1:3) :: d2
   real, dimension(1:3) :: R    !relativna poloha 2. dipolu vzhladom k prvemu
   real, dimension(1:3) :: r_0
 
   integer :: errorcode,i
   real :: V

!konstanty
   real :: eps = 1e-20 !minimalna vzdialenost
   real :: debye = 3.33564e-30
   real :: angstrom = 1e-10
   real :: epsilon0 = 8.854e-12
   real :: K = 8.987551787e9
   real :: h = 6.626068e-34
   real :: c = 299792458


!nacitanie hodnot
  
  !nacitanie hodnot
! cteni r1
   do i = 1,3
      read(*,*) r1(i)   
      read(*,*) r2(i)
      read(*,*) d1(i)
      read(*,*) d2(i)
   end do
 
  
!vypocet

   call DipoleEnergy(r1,d1,r2,d2,errorcode,V);
   
!vypis    
   write(*,*) 'err', errorcode
   write(*,*) 'vvv', V
   


!funkcie
contains

   function V_Norm(a)
     real :: V_Norm
     real, dimension(1:3) :: a
     V_Norm = Sqrt(dot_product(a,a))
   end function V_Norm
   
   subroutine DipoleEnergy(r1,d1,r2,d2,errorcode,V)
     real, dimension(1:3) :: r1   
     real, dimension(1:3) :: d1  
     real, dimension(1:3) :: r2   
     real, dimension(1:3) :: d2       
     real, dimension(1:3) :: R 
     real, dimension(1:3) :: r0
     integer :: errorcode   
     real :: RR      
     real :: V
     
     
     !vypocet
     
     errorcode = -1
     R = r1 - r2  ! relativna poloha 2. vzhladom na prvy
     RR = V_Norm(R)
     
     if (RR > eps) then
       errorcode = 0        
       r0 = R/RR  !jednotkovy vektor
       V = 5034.1153*(dot_product(d1,d2)-3*dot_product(d1,r0)*dot_product(d2,r0))/(RR**3)
       !5034.1153 = c*1e-7*Debye**2/100/h/angstrom**3
     end if
     
  end subroutine DipoleEnergy    


end program dipolyf     
