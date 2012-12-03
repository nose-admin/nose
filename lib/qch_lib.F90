#include "util_allocation.h"

module qch_lib

	use std_types
	use std_io
	use numer_matrix
	use sci_misc

	use util_allocation

	implicit none

	integer(i4b), parameter, private				:: gaussians_per_slater = 3

	real(dp), private 								:: exponent_slaterP = 1.625_dp
	real(dp), private 								:: exponent_slaterS = 1.24_dp
	real(dp), private, parameter					:: shadow_factor   = 1.0_dp!3.350_dp
	real(dp), private, dimension(gaussians_per_slater) 	:: 	STOxiS, STOxiP, STOciP, STOciS

	integer(i4b), parameter, private				:: depth_of_recursion_of_4center = 25
	real(dp), parameter, private					:: hartree_in_eV = 27.21138386_dp
	real(dp), parameter, private					:: hartree_in_invCm = 219474.66_dp
	real(dp), parameter, private					:: bohr_radius_in_pm = 52.9177_dp
	real(dp), parameter, private					:: nuclei_shift_factor = 0.05_dp
	real(dp), parameter, private					:: Angstrom_in_Bohr_radii = 1.88972688_dp

	public::init
	private::overlap_SS
	private::overlap_PiS
	private::overlap_PiPj
	private::overlap
	private::density_S
	private::density_Pi
	private::density
	private::dipole_xi_SS
	private::dipole_xi_PiS
	private::dipole_xi_PiPj
	private::dipole_xi
	private::factorG
	private::dfactorG
	private::factorF
	private::dfactorF
	private::factorGcoul
	private::dfactorGcoul
	private::factorFcoul
	private::dfactorFcoul
	private::four_center_integral_inner
	private::four_center_integral
	private::kinetic_integral
	private::coulombic_integral
	public::four_center_integral_slater
	public::kinetic_integral_slater
	public::coulombic_integral_slater
	public::overlap_slater
	public::dipole_xi_slater
	public::density_slater
	public::four_center_integral_slater_directed
	public::kinetic_integral_slater_directed
	public::coulombic_integral_slater_directed
	public::overlap_slater_directed
	public::dipole_xi_slater_directed

	private::fill_four_center_tensor_and_H_S_slater_directed
	private::fill_four_center_tensor_and_H_S_slater
	private::fill_F_matrix
	private::fill_P_matrix
	private::fill_HOMO_state
	private::fill_LUMO_state
	public::solve_Roothan_equations
	public::full_energy
	public::full_energy_Huckel
	private::excitation_difference_code
	private::interaction_energy_inside_one_molecule
	public::interaction_energy
	public::calculate_Huckel_method
	public::fill_dipole_operator_elements_slater
	public::fill_dipole_operator_elements_slater_directed
	private::fill_intermolecular_four_center_tensor_and_H_slater_directed
	private::fill_intermolecular_four_center_tensor_and_H_slater
	public::test_hydrogen_molecule
	public::prepare_carotenoid
	public::prepare_alkene
	public::read_config_file
	private::fit_Huckel_parameters
	private::fitness_for_Huckel_parameters_function


 	contains

 	subroutine init(slater_expS, slater_expP)
 		real(dp), intent(in)		:: slater_expS, slater_expP
 		real(dp), dimension(3) 	:: vecA

		exponent_slaterS = slater_expS
		exponent_slaterP = slater_expP

 		STOxiS  = (/ 0.185742_dp, 1.12_dp,  13.5413_dp /) * exponent_slaterS**2
 		STOxiP  = (/ 0.185742_dp, 1.12_dp,  13.5413_dp /) * exponent_slaterP**2
 		STOciS = (/ 0.276651_dp,  0.433794_dp,  0.242884_dp /)

 		vecA = 0.0_dp

!		write(*,*) '  STOciS0', STOciS
!		write(*,*) '  overlap_slater', overlap_slater(vecA,vecA,0,0)

 		STOciS = STOciS / sqrt(overlap_slater(vecA,vecA,0,0))
 		STOciP = STOciS
 		STOciP = STOciP / sqrt(overlap_slater(vecA,vecA,1,1))

!		write(*,*) '  exp_slater', exponent_slater
!		write(*,*) '  STOciS', STOciS
!		write(*,*) '  STOciP', STOciP
 	end subroutine init

	real(dp) function overlap_SS(vecA, vecB, a, b) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b

		res = exp(-a*b/(a+b)*dot_product(vecA-vecB,vecA-vecB))
		res = res*((PI/(a+b))**(3.0/2.0))

 	end function overlap_SS

	real(dp) function overlap_PiS(vecA, vecB, a, b, i) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i

		if(i < 1 .or. i > 3) then
			call print_error_message(-1,"error in overlap_PiS()")
		end if

		res = overlap_SS(vecA, vecB, a, b)
		res = res*b/(a + b)*(vecB(i) - vecA(i))

 	end function overlap_PiS

	real(dp) function overlap_PiPj(vecA,vecB, a,b, i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j

		if(i < 1 .or. i > 3 .or. j < 1 .or. j > 3) then
			call print_error_message(-1,"error in overlap_PiS()")
		end if

		if(i == j) then
			res = overlap_SS(vecA, vecB, a, b)
			res = res*(a + b - 2*a*b*(vecB(i) - vecA(i))*(vecB(j) - vecA(j)))/(2*((a + b)**2))
		else
			res = overlap_SS(vecA, vecB, a, b)
			res = res*(-2*a*b*(vecB(i) - vecA(i))*(vecB(j) - vecA(j)))/(2*((a + b)**2))
		end if

 	end function overlap_PiPj

	real(dp) function overlap(vecA,vecB, a,b, i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j

		if(i < 0 .or. i > 3 .or. j < 0 .or. j > 3) then
			call print_error_message(-1,"error in overlap()")
		end if

		if(i < 1 .and. j < 1) then
			res = overlap_SS(vecA, vecB, a, b)
		elseif(i < 1) then
			res = overlap_PiS(vecB, vecA, b, a, j)
		elseif(j < 1) then
			res = overlap_PiS(vecA, vecB, a, b, i)
		else
			res = overlap_PiPj(vecA, vecB, a, b, i, j)
		end if

 	end function overlap

	real(dp) function density_S(vecA, vecB, a) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a

		res = exp(-a*dot_product(vecA-vecB,vecA-vecB))

 	end function density_S

	real(dp) function density_Pi(vecA,vecB, a, i) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a
		integer(i4b), intent(in)				:: i

		if(i < 1 .or. i > 3) then
			call print_error_message(-1,"error in density_PiS()")
		end if

		res = density_S(vecA, vecB, a)*(vecA(i)-vecB(i))


 	end function density_Pi

	real(dp) function density(vecA,vecB, a,i) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a
		integer(i4b), intent(in)				:: i

		if(i < 0 .or. i > 3) then
			call print_error_message(-1,"error in density()")
		end if

		if(i < 1) then
			res = density_S(vecA, vecB, a)
		else
			res = density_Pi(vecA, vecB, a, i)
		end if

 	end function density

	real(dp) function dipole_xi_SS(vecA, vecB, a, b, dipole_index) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b	), intent(in)				:: dipole_index

		if(dipole_index < 1 .or. dipole_index > 3) then
			call print_error_message(-1, "dimension error in dipole_xi_SS")
		end if

		res = overlap(vecA,vecB,a,b,dipole_index,0) + overlap(vecA,vecB,a,b,0,0)*vecA(dipole_index)

 	end function dipole_xi_SS

	real(dp) function dipole_xi_PiS(vecA, vecB, a, b, i, dipole_index) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i, dipole_index

		if(dipole_index < 1 .or. dipole_index > 3) then
			call print_error_message(-1, "dimension error in dipole_xi_PiS")
		end if

		res = overlap(vecA,vecB,a,b,i,dipole_index) + overlap(vecA,vecB,a,b,i,0)*vecB(dipole_index)

 	end function dipole_xi_PiS

	real(dp) function dipole_xi_PiPj(vecA,vecB, a,b, i,j, dipole_index) result(res)
		! not tested
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j, dipole_index

		integer(i4b) 	:: k

		if(i < 1 .or. i > 3 .or. j < 1 .or. j > 3) then
			call print_error_message(-1,"error in overlap_PiS()")
		end if

		res = overlap(vecA, vecB, a, b, i, j)
		if(i == j .and. i == dipole_index) then
			res = res*(-((b*(-3*b + a*(-3 + 2*b*(vecA(i) - vecB(i))**2))*(vecA(i) - vecB(i)))/&
					((a + b)*(-b + a*(-1 + 2*b*(vecA(i) - vecB(i))**2)))))
		elseif(i /= dipole_index .and. j /= dipole_index) then
			res = res*((b*(-vecA(dipole_index) + vecB(dipole_index)))/(a + b))
		elseif(i /= j .and. (i == dipole_index .or. j == dipole_index)) then
			if(i == dipole_index) then
				k = j
			else
				k = i
			end if

			res = ((b*(-b + a*(-1 + 2*b*(vecA(dipole_index) - vecB(dipole_index))**2))*(vecA(k) - vecB(k)))/(2.0_dp*(a + b)**3))
!			res = res*((a + b - 2*a*b*(vecA(dipole_index) - vecB(dipole_index))**2)/&
!					(2.0_dp*a*(a + b)*(vecA(dipole_index) - vecB(dipole_index))))
		end if

		res = res + vecA(dipole_index)*overlap(vecA, vecB, a, b, i, j)

 	end function dipole_xi_PiPj

	real(dp) function dipole_xi(vecA,vecB, a,b, i,j, dipole_index) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j,dipole_index

		if(i < 0 .or. i > 3 .or. j < 0 .or. j > 3) then
			call print_error_message(-1,"error in overlap()")
		end if

		if(i < 1 .and. j < 1) then
			res = dipole_xi_SS(vecA, vecB, a, b, dipole_index)
		elseif(i < 1) then
			res = dipole_xi_PiS(vecB, vecA, b, a, j, dipole_index)
		elseif(j < 1) then
			res = dipole_xi_PiS(vecA, vecB, a, b, i, dipole_index)
		else
			res = dipole_xi_PiPj(vecA, vecB, a, b, i, j, dipole_index)
		end if

 	end function dipole_xi

 	real(dp) function factorF(vecA,vecB,vecC,vecD, a,b,c,d) result(res)
 		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), intent(in)					:: a,b,c,d

		res = 	a*b/(a+b)*dot_product(vecA-vecB,vecA-vecB) + &
					c*d/(c+d)*dot_product(vecC-vecD,vecC-vecD)
 	end function factorF

 	real(dp) function dfactorF(vecA,vecB,vecC,vecD, a,b,c,d, dtf) result(res)
 		real(dp), dimension(3), intent(in) 		:: vecA, vecB, vecC, vecD
		real(dp), intent(in)						:: a,b,c,d
		integer(i4b), dimension(4,3), intent(in)	:: dtf

		integer(i4b)								:: i,j,k,l,  coord

		res = 0

		!if two of coordinate derivatives are nonzero, result is zero
		if(	sum(dtf(:,1)) > 0 .and. sum(dtf(:,2)) > 0 .or. &
			sum(dtf(:,1)) > 0 .and. sum(dtf(:,3)) > 0 .or. &
			sum(dtf(:,2)) > 0 .and. sum(dtf(:,3)) > 0 ) then

			res = 0.0_dp
			return
		end if

		if(sum(dtf) < 1) then
			res = 1.0_dp
			return
		end if

		do i=1,3
			if(sum(dtf(:,i)) > 0) then
				coord = i
			end if
		end do

		i = dtf(1,coord)
		j = dtf(2,coord)
		k = dtf(3,coord)
		l = dtf(4,coord)

		if(i+j+k+l > 2) then
			res = 0
			return
		else
			if(i == 2) then
				res = 2*a*b/(a+b)
				return
			else if(j == 2) then
				res = 2*a*b/(a+b)
				return
			else if(k == 2) then
				res = 2*c*d/(c+d)
				return
			else if(l == 2) then
				res = 2*c*d/(c+d)
				return
			else if(i == 1 .and. j == 1) then
				res = -2*a*b/(a+b)
				return
			else if(i == 1 .and. k == 1) then
				res = 0.0_dp
				return
			else if(i == 1 .and. l == 1) then
				res = 0.0_dp
				return
			else if(j == 1 .and. k == 1) then
				res = 0.0_dp
				return
			else if(j == 1 .and. l == 1) then
				res = 0.0_dp
				return
			else if(k == 1 .and. l == 1) then
				res = -2*c*d/(c+d)
				return

			else if(i == 1) then
				res = 2*a*b/(a+b)*(vecA(coord)-vecB(coord))
			else if(j == 1) then
				res = -2*a*b/(a+b)*(vecA(coord)-vecB(coord))
			else if(k == 1) then
				res = 2*c*d/(c+d)*(vecC(coord)-vecD(coord))
			else if(l == 1) then
				res = -2*c*d/(c+d)*(vecC(coord)-vecD(coord))
			end if
		end if
 	end function dfactorF

 	real(dp) function factorG(vecA,vecB,vecC,vecD, a,b,c,d) result(res)
 		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), intent(in)					:: a,b,c,d

		real(dp), dimension(3)					:: vecP, vecQ

		vecP = (a*vecA + b*vecB)/(a + b)
		vecQ = (c*vecC + d*vecD)/(c + d)
		res = (a+b)*(c+d)/(a+b+c+d)*dot_product(vecP-vecQ,vecP-vecQ)
 	end function factorG

 	real(dp) function dfactorG(vecA,vecB,vecC,vecD, a,b,c,d, dtg) result(res)
 		real(dp), dimension(3), intent(in) 		:: vecA, vecB, vecC, vecD
		real(dp), intent(in)						:: a,b,c,d
		integer(i4b), dimension(4,3), intent(in)	:: dtg

		integer(i4b)								:: i,j,k,l,  coord

		res = 0

		!if two of coordinate derivatives are nonzero, result is zero
		if(	sum(dtg(:,1)) > 0 .and. sum(dtg(:,2)) > 0 .or. &
			sum(dtg(:,1)) > 0 .and. sum(dtg(:,3)) > 0 .or. &
			sum(dtg(:,2)) > 0 .and. sum(dtg(:,3)) > 0) then

			res = 0
			return
		end if

		if(sum(dtg) < 1) then
			res = 1.0_dp
			return
		end if

		do i=1,3
			if(sum(dtg(:,i)) > 0) then
				coord = i
			end if
		end do

		i = dtg(1,coord)
		j = dtg(2,coord)
		k = dtg(3,coord)
		l = dtg(4,coord)

		if(i+j+k+l > 2) then
			res = 0
			return
		else
			if(i == 2) then
				res = (2*a**2*(c + d))/((a + b)*(a + b + c + d))
				return
			else if(j == 2) then
				res = (2*b**2*(c + d))/((a + b)*(a + b + c + d))
				return
			else if(k == 2) then
				res = (2*(a + b)*c**2)/((c + d)*(a + b + c + d))
				return
			else if(l == 2) then
				res = (2*(a + b)*d**2)/((c + d)*(a + b + c + d))
				return
			else if(i == 1 .and. j == 1) then
				res = (2*a*b*(c + d))/((a + b)*(a + b + c + d))
				return
			else if(i == 1 .and. k == 1) then
				res = (-2*a*c)/(a + b + c + d)
				return
			else if(i == 1 .and. l == 1) then
				res = (-2*a*d)/(a + b + c + d)
				return
			else if(j == 1 .and. k == 1) then
				res = (-2*b*c)/(a + b + c + d)
				return
			else if(j == 1 .and. l == 1) then
				res = (-2*b*d)/(a + b + c + d)
				return
			else if(k == 1 .and. l == 1) then
				res = (2*(a + b)*c*d)/((c + d)*(a + b + c + d))
				return

			else if(i == 1) then
				res = (2*a*(c + d)*((a*vecA(coord) + b*vecB(coord))/(a + b) - &
						(c*vecC(coord) + d*vecD(coord))/(c + d)))/(a + b + c + d)
			else if(j == 1) then
				res = (2*b*(c + d)*((a*vecA(coord) + b*vecB(coord))/(a + b) - &
						(c*vecC(coord) + d*vecD(coord))/(c + d)))/(a + b + c + d)
			else if(k == 1) then
				res = (-2*(a + b)*c*((a*vecA(coord) + b*vecB(coord))/(a + b) - &
						(c*vecC(coord) + d*vecD(coord))/(c + d)))/(a + b + c + d)
			else if(l == 1) then
				res = (-2*(a + b)*d*((a*vecA(coord) + b*vecB(coord))/(a + b) - &
						(c*vecC(coord) + d*vecD(coord))/(c + d)))/(a + b + c + d)
			end if
		end if
 	end function dfactorG

 	real(dp) function factorFcoul(vecA,vecB,vecC,vecD, a,b,c,d) result(res)
 		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), intent(in)					:: a,b,c,d

		res = 	a*b/(a+b)*dot_product(vecA-vecB,vecA-vecB)
 	end function factorFcoul

 	real(dp) function dfactorFcoul(vecA,vecB,vecC,vecD, a,b,c,d, dtf) result(res)
 		real(dp), dimension(3), intent(in) 		:: vecA, vecB, vecC, vecD
		real(dp), intent(in)						:: a,b,c,d
		integer(i4b), dimension(4,3), intent(in)	:: dtf

		integer(i4b)								:: i,j,k,l,  coord

		res = 0

		!this is the only difference between non-coulombic F
		if(sum(dtf(3,:)) + sum(dtf(4,:)) > 0) then
			res = 0.0_dp
			return
		end if

		!if two of coordinate derivatives are nonzero, result is zero
		if(	sum(dtf(:,1)) > 0 .and. sum(dtf(:,2)) > 0 .or. &
			sum(dtf(:,1)) > 0 .and. sum(dtf(:,3)) > 0 .or. &
			sum(dtf(:,2)) > 0 .and. sum(dtf(:,3)) > 0 ) then

			res = 0.0_dp
			return
		end if

		if(sum(dtf) < 1) then
			res = 1.0_dp
			return
		end if

		do i=1,3
			if(sum(dtf(:,i)) > 0) then
				coord = i
			end if
		end do

		i = dtf(1,coord)
		j = dtf(2,coord)
		k = dtf(3,coord)
		l = dtf(4,coord)

		if(i+j+k+l > 2) then
			res = 0
			return
		else
			if(i == 2) then
				res = 2*a*b/(a+b)
				return
			else if(j == 2) then
				res = 2*a*b/(a+b)
				return
			else if(k == 2) then
				res = 2*c*d/(c+d)
				return
			else if(l == 2) then
				res = 2*c*d/(c+d)
				return
			else if(i == 1 .and. j == 1) then
				res = -2*a*b/(a+b)
				return
			else if(i == 1 .and. k == 1) then
				res = 0.0_dp
				return
			else if(i == 1 .and. l == 1) then
				res = 0.0_dp
				return
			else if(j == 1 .and. k == 1) then
				res = 0.0_dp
				return
			else if(j == 1 .and. l == 1) then
				res = 0.0_dp
				return
			else if(k == 1 .and. l == 1) then
				res = -2*c*d/(c+d)
				return

			else if(i == 1) then
				res = 2*a*b/(a+b)*(vecA(coord)-vecB(coord))
			else if(j == 1) then
				res = -2*a*b/(a+b)*(vecA(coord)-vecB(coord))
			else if(k == 1) then
				res = 2*c*d/(c+d)*(vecC(coord)-vecD(coord))
			else if(l == 1) then
				res = -2*c*d/(c+d)*(vecC(coord)-vecD(coord))
			end if
		end if
 	end function dfactorFcoul

 	real(dp) function factorGcoul(vecA,vecB,vecC,vecD, a,b,c,d) result(res)
 		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), intent(in)					:: a,b,c,d

		real(dp), dimension(3)					:: vecP

		vecP = (a*vecA + b*vecB)/(a + b)
		res = (a+b)*dot_product(vecP-vecC,vecP-vecC)
 	end function factorGcoul

 	real(dp) function dfactorGcoul(vecA,vecB,vecC,vecD, a,b,c,d, dtg) result(res)
 		real(dp), dimension(3), intent(in) 		:: vecA, vecB, vecC, vecD
		real(dp), intent(in)						:: a,b,c,d
		integer(i4b), dimension(4,3), intent(in)	:: dtg

		integer(i4b)								:: i,j,k,l,  coord

		res = 0

		!! derivatives according to vecC (coulombic center) are not included !!

		!if two of coordinate derivatives are nonzero, result is zero
		if(	sum(dtg(:,1)) > 0 .and. sum(dtg(:,2)) > 0 .or. &
			sum(dtg(:,1)) > 0 .and. sum(dtg(:,3)) > 0 .or. &
			sum(dtg(:,2)) > 0 .and. sum(dtg(:,3)) > 0) then

			res = 0
			return
		end if

		if(sum(dtg) < 1) then
			res = 1.0_dp
			return
		end if

		do i=1,3
			if(sum(dtg(:,i)) > 0) then
				coord = i
			end if
		end do

		i = dtg(1,coord)
		j = dtg(2,coord)
		k = dtg(3,coord)
		l = dtg(4,coord)

		if(l > 0) then
			call print_error_message(-1,"dfactorGcoul - attempted to derive according to d")
		end if

		if(i+j+k+l > 2) then
			res = 0
			return
		else
			if(i == 2) then
				res = 2*a**2/(a+b)
				return
			else if(j == 2) then
				res = 2*b**2/(a+b)
				return
			else if(k == 2) then
				res = 2*(a+b)
				return
			else if(i == 1 .and. j == 1) then
				res = 2*a*b/(a+b)
				return
			else if(i == 1 .and. k == 1) then
				res = -2*a
				return
			else if(j == 1 .and. k == 1) then
				res = -2*b
				return

			else if(i == 1) then
				res = 2*a*((a*vecA(coord) + b*vecB(coord))/(a + b) - vecC(coord))
			else if(j == 1) then
				res = 2*b*((a*vecA(coord) + b*vecB(coord))/(a + b) - vecC(coord))
			else if(k == 1) then
				res = -2*(a*vecA(coord) + b*vecB(coord)) + 2*(a+b)*vecC(coord)
			end if
		end if
 	end function dfactorGcoul

	! represents element D^todo [(D^dtg g) (D^dtf f) (erff ? erf(sqrt(g)) : e^(-g) ) g^(gi/2)]
 	real(dp) recursive function four_center_integral_inner(vA,vB,vC,vD, a,b,c,d,&
 						gi,erff,dtg,dtf,todo,level,coulombic) result(res)
		real(dp), dimension(3), intent(in) 	:: vA, vB, vC, vD
		real(dp), intent(in)					:: a,b,c,d
		integer(i4b), intent(in)	:: gi, level
		logical, intent(in)		:: erff, coulombic
		integer(i4b), intent(inout), dimension(4,3)		:: todo
		integer(i4b), intent(inout), dimension(depth_of_recursion_of_4center,4,3)	&
														:: dtg, dtf ! power, center, {x,y,z}

		integer(i4b)	:: aa,bb,ll
		real(dp), dimension(3) 	:: vecA, vecB, vecC, vecD

		! we don't want to spoil them, but we need to move them in some occasions
		vecA = vA
		vecB = vB
		vecC = vC
		vecD = vD

		res = 0.0_dp

		! errors
		if(a <= 1e-6 .or. b <= 1e-6 .or. c <= 1e-6 .or. d <= 1e-6) then
			call print_error_message(-1,'gaussian xi is zero in four_center_integral_inner')
		endif

		!! special cases

		! factorF-dependence is exponential, we neglect the result  above some value
		if((factorF(vecA,vecB,vecC,vecD,a,b,c,d) > 15.0 .and. .not. coulombic) .or.&
			(factorFcoul(vecA,vecB,vecC,vecD,a,b,c,d) > 15.0 .and. coulombic)) then
			res = 0.0_dp
			return
		end if

		! We do not know how to handle cases when g ~ 0 and we need to mage derivatives
		! according to it. Therefore, if this happens, we simply move problematic nuclei
		! a bit further.
		do while(((.not. coulombic) .and. factorG(vecA,vecB,vecC,vecD,a,b,c,d) .le. 1.0e-5_dp)&
					.or. (coulombic .and. factorGcoul(vecA,vecB,vecC,vecD,a,b,c,d) .le. 1.0e-5_dp))

			vecC(3) = vecC(3) + nuclei_shift_factor/sqrt((a+b)*(c+d)/(a+b+c+d))
			vecD(3) = vecD(3) + nuclei_shift_factor/sqrt((a+b)*(c+d)/(a+b+c+d))

			if(coulombic .and. maxval(abs(vecA-vecB)) < 1e-5 .and. maxval(abs(vecA-vecC)) < 1e-5 .and. maxval(abs(vecA-vecD)) < 1e-5 .and. maxval(abs(vecA-vecB)) < 1e-5) then
				write(*,*) 'C'
				write(*,*) vecA
				write(*,*) vecB
				write(*,*) vecC
				write(*,*) vecD
				write(*,*) a,b,c,d
			end if

		end do


		! no further derivative to be done
		if(sum(todo) == 0) then
			if(coulombic) then
				res = 1.0_dp
				do ll=1,level-1
					res = res * dfactorGcoul(vecA,vecB,vecC,vecD,a,b,c,d,dtg(ll,:,:))
					res = res * dfactorFcoul(vecA,vecB,vecC,vecD,a,b,c,d,dtf(ll,:,:))
				end do

				if(factorGcoul(vecA,vecB,vecC,vecD,a,b,c,d) .le. 1.0e-5_dp) then
					! anomalous case
					if(.not. erff) then
						! These terms are actually very important for the result,
						! but all of them contribute to the final result. The result
						! is calculated in term that includes erf(), because there is
						! only one such term, while there is more other terms.
						res = 0.0_dp
						return
					else
						write(*,*) 'C'
						write(*,*) vecA
						write(*,*) vecB
						write(*,*) vecC
						write(*,*) vecD
						write(*,*) a,b,c,d
						write(*,*) factorG(vecA,vecB,vecC,vecD,a,b,c,d)
						call print_error_message(-1,"Anomalous case of iterative procedure not implemented - this should't happen")
						! we replace erf(sqrt(g))/sqrt(g) by its limit value for g->0
 						res = res*2.0_dp/sqrt(PI) ! works only for S-orbitals
					end if
				else
					! normal case
					res = res * factorGcoul(vecA,vecB,vecC,vecD,a,b,c,d)**(gi/2.0_dp)

					if(erff) then
						res = res * erf(sqrt(factorGcoul(vecA,vecB,vecC,vecD,a,b,c,d)))
					else
						res = res * exp(-factorGcoul(vecA,vecB,vecC,vecD,a,b,c,d))
					end if
				end if

				res = res * exp(-factorFcoul(vecA,vecB,vecC,vecD,a,b,c,d))
			else
				res = 1.0_dp
				do ll=1,level-1
					res = res * dfactorG(vecA,vecB,vecC,vecD,a,b,c,d,dtg(ll,:,:))
					res = res * dfactorF(vecA,vecB,vecC,vecD,a,b,c,d,dtf(ll,:,:))
				end do

				if(factorG(vecA,vecB,vecC,vecD,a,b,c,d) .le. 1.0e-5_dp) then
					! anomalous case
					if(.not. erff) then
						! These terms are actually very important for the result,
						! but all of them contribute to the final result. The result
						! is calculated in term that includes erf(), because there is
						! only one such term, while there is more other terms.
						res = 0.0_dp
						return
					else
						write(*,*) '4'
						write(*,*) vecA
						write(*,*) vecB
						write(*,*) vecC
						write(*,*) vecD
						write(*,*) a,b,c,d
						write(*,*) factorG(vecA,vecB,vecC,vecD,a,b,c,d)
						call print_error_message(-1,"Anomalous case of iterative procedure not implemented - this should't happen")
						! we replace erf(sqrt(g))/sqrt(g) by its limit value for g->0
 						res = res*2.0_dp/sqrt(PI) ! works only for S-orbitals
					end if
				else
					! normal case
					res = res * factorG(vecA,vecB,vecC,vecD,a,b,c,d)**(gi/2.0_dp)

					if(erff) then
						res = res * erf(sqrt(factorG(vecA,vecB,vecC,vecD,a,b,c,d)))
					else
						res = res * exp(-factorG(vecA,vecB,vecC,vecD,a,b,c,d))
					end if
				end if

				res = res * exp(-factorF(vecA,vecB,vecC,vecD,a,b,c,d))
			end if

			return

		else
			do aa=1,4
			do bb=1,3

				if(todo(aa,bb) >= 1) then
					todo(aa,bb) = todo(aa,bb) - 1
					if(erff) then
						!derivative of erf
						dtg(level,aa,bb) = dtg(level,aa,bb) + 1
						res = res + four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
									gi-1,.false.,dtg,dtf,todo,level+1,coulombic)/sqrt(PI)
						dtg(level,aa,bb) = dtg(level,aa,bb) - 1

					else
						!derivative of e^-g
						dtg(level,aa,bb) = dtg(level,aa,bb) + 1
						res = res - four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
									gi,.false.,dtg,dtf,todo,level+1,coulombic)
						dtg(level,aa,bb) = dtg(level,aa,bb) - 1

					end if

					!derivative of e^-f
					dtf(level,aa,bb) = dtf(level,aa,bb) + 1
					res = res - four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
								gi,erff,dtg,dtf,todo,level+1,coulombic)
					dtf(level,aa,bb) = dtf(level,aa,bb) - 1

					!derivative of g^i/2
					dtg(level,aa,bb) = dtg(level,aa,bb) + 1
					res = res + four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
								gi-2,erff,dtg,dtf,todo,level+1,coulombic)*gi/2.0_dp
					dtg(level,aa,bb) = dtg(level,aa,bb) - 1

					do ll=1,level-1
						!derivative of dtf
						if(sum(dtf(ll,:,:)) > 0) then
						dtf(ll,aa,bb) = dtf(ll,aa,bb) + 1
						res = res + four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
									gi,erff,dtg,dtf,todo,level+1,coulombic)
						dtf(ll,aa,bb) = dtf(ll,aa,bb) - 1
						end if

						!derivative of dtg
						if(sum(dtg(ll,:,:)) > 0) then
						dtg(ll,aa,bb) = dtg(ll,aa,bb) + 1
						res = res + four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
									gi,erff,dtg,dtf,todo,level+1,coulombic)
						dtg(ll,aa,bb) = dtg(ll,aa,bb) - 1
						end if
					end do

					todo(aa,bb) = todo(aa,bb) + 1
					return
				end if
			end do
			end do
		end if

 	end function four_center_integral_inner

 	real(dp) function four_center_integral(vecA,vecB,vecC,vecD, a,b,c,d, &
 									i,j,k,l) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), intent(in)					:: a,b,c,d
		integer(i4b), intent(in)				:: i,j,k,l

		integer(i4b), dimension(4,3)									:: todo
		integer(i4b), dimension(depth_of_recursion_of_4center,4,3)	:: dtg, dtf

		dtg = 0
		dtf = 0
		todo = 0
		res = 1.0_dp

		! 0 denotes s-orbital, 1,2,3 denote axis
		if(i > 0) then
			todo(1,i) = 1
			res = res / 2.0_dp / a
		end if
		if(j > 0) then
			todo(2,j) = 1
			res = res / 2.0_dp / b
		end if
		if(k > 0) then
			todo(3,k) = 1
			res = res / 2.0_dp / c
		end if
		if(l > 0) then
			todo(4,l) = 1
			res = res / 2.0_dp / d
		end if

		res = res * four_center_integral_inner(vecA,vecB,vecC,vecD,a,b,c,d,&
						-1,.true.,dtg,dtf,todo,1,.false.)

		res = res * PI**3/(a+b)/(c+d)/sqrt(a+b+c+d)

 	end function four_center_integral

 	real(dp) function coulombic_integral(vecA,vecB,vecCoul, a,b, i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecCoul
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j

		integer(i4b), dimension(4,3)									:: todo
		integer(i4b), dimension(depth_of_recursion_of_4center,4,3)	:: dtg, dtf

		real(dp), dimension(3) 	:: vecD

		dtg = 0
		dtf = 0
		todo = 0
		res = 1.0_dp

		vecD = 0.0_dp

		! 0 denotes s-orbital, 1,2,3 denote axis
		if(i > 0) then
			todo(1,i) = 1
			res = res / 2.0_dp / a
		end if
		if(j > 0) then
			todo(2,j) = 1
			res = res / 2.0_dp / b
		end if

		res = res * four_center_integral_inner(vecA,vecB,vecCoul,vecD,a,b,1.0_dp,1.0_dp,&
						-1,.true.,dtg,dtf,todo,1,.true.)
		res = res * 2*PI/(a+b) * PI**0.5_dp / 2.0_dp

 	end function coulombic_integral

 	real(dp) function kinetic_integral(vecA,vecB, a,b, i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), intent(in)					:: a,b
		integer(i4b), intent(in)				:: i,j

		real(dp)								:: FFF

		res = 1.0_dp

		if(i > 0) then
			res = res / 2.0_dp / a
		end if
		if(j > 0) then
			res = res / 2.0_dp / b
		end if

		! This function was "written" by copy-and-pasting the code from
		! Mathematica. It works well, but could be written in more human-readable
		! manner.

		if(i == j .and. i == 0) then
			res = (a * b * (1/(a + b))**2.5_dp * PI**1.5_dp *	&
				(3 - (2*a*b*((vecA(1) - vecB(1))**2 + 			&
				(vecA(2) - vecB(2))**2 + 						&
				(vecA(3) - vecB(3))**2))/(a + b)))/				&
				exp((a*b*((vecA(1) - vecB(1))**2 + 				&
				(vecA(2) - vecB(2))**2 + (vecA(3) - vecB(3))**2 &
				))/(a + b))  * res
		else if(i == 0 .or. j == 0) then
			! result is calculated for derivative according to x1

			res = (2*a**2*b**2*(1/(a + b))**4.5_dp*PI**1.5_dp *	&
				(vecA(1) - vecB(1))*							&
				(-5*b + a*(-5 + 								&
				2*b*((vecA(1) - vecB(1))**2 + 					&
				(vecA(2) - vecB(2))**2 + 						&
				(vecA(3) - vecB(3))**2))))/						&
				exp((a*b*((vecA(1) - vecB(1))**2 + 				&
				(vecA(2) - vecB(2))**2 + (vecA(3) - vecB(3))**2 &
				))/(a + b))  * res

			res = res * (vecA(max(i,j)) - vecB(max(i,j)))/(vecA(1) - vecB(1))

			if(j > i) then
				res = res*(-1)
			end if
		else
			FFF = dot_product(vecA-vecB,vecA-vecB)

			if(i == j) then
				res = (2*a**2*b**2*(1/(a + b))**5.5_dp*PI**1.5_dp*		&
					(5*b**2 - 2*a*b*									&
					(-5 + b*(FFF + 7*(vecA(i) - vecB(j))**2)) + 		&
					a**2*(5 - 2*b*(FFF + 7*(vecA(i) - vecB(j))**2) + 	&
					4*b**2*FFF*(vecA(i) - vecB(j))**2)))/				&
					exp((a*b*FFF)/(a + b))  * res
			else
				res = (4*a**3*b**3*(1/(a + b))**5.5_dp*			&
					(-7*(a + b) + 2*a*b*FFF)*PI**1.5_dp*		&
					(vecA(i) - vecB(i))*(vecA(j) - vecB(j)))/	&
					exp((a*b*FFF)/(a + b))  * res
			end if

			! conversion to all other results follow

		end if

 	end function kinetic_integral

 	real(dp) function four_center_integral_slater(vecA,vecB,vecC,vecD,i,j,k,l) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		integer(i4b), intent(in)				:: i,j,k,l

		integer(i4b)							:: r,s,t,u
		real(dp),dimension(gaussians_per_slater) :: STOci1, STOci2, STOci3, STOci4
		real(dp),dimension(gaussians_per_slater) :: STOxi1, STOxi2, STOxi3, STOxi4

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiP
		else
			call print_error_message(-1, "wrong i in four_center_integral_slater")
		end if
		if(j == 0) then
			STOci2 = STOciS
			STOxi2 = STOxiS
		elseif(j >= 1 .and. j <= 3) then
			STOci2 = STOciP
			STOxi2 = STOxiP
		else
			call print_error_message(-1, "wrong j in four_center_integral_slater")
		end if
		if(k == 0) then
			STOci3 = STOciS
			STOxi3 = STOxiS
		elseif(k >= 1 .and. k <= 3) then
			STOci3 = STOciP
			STOxi3 = STOxiP
		else
			call print_error_message(-1, "wrong k in four_center_integral_slater")
		end if
		if(l == 0) then
			STOci4 = STOciS
			STOxi4 = STOxiS
		elseif(l >= 1 .and. l <= 3) then
			STOci4 = STOciP
			STOxi4 = STOxiP
		else
			call print_error_message(-1, "wrong l in four_center_integral_slater")
		end if

		res = 0.0_dp

		do r=1,gaussians_per_slater
		do s=1,gaussians_per_slater
		do t=1,gaussians_per_slater
		do u=1,gaussians_per_slater

			res = res + four_center_integral(vecA,vecB,vecC,vecD, STOxi1(r),STOxi2(s),STOxi3(t),STOxi4(u), &
							i,j,k,l)  *  STOci1(r)*STOci2(s)*STOci3(t)*STOci4(u)

		end do
		end do
		end do
		end do

 	end function four_center_integral_slater

 	real(dp) function kinetic_integral_slater(vecA,vecB,i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		integer(i4b), intent(in)				:: i,j

		integer(i4b)							:: r,s
		real(dp),dimension(gaussians_per_slater) :: STOci1, STOci2, STOxi1, STOxi2

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiP
		else
			call print_error_message(-1, "wrong i in kinetic_integral_slater")
		end if
		if(j == 0) then
			STOci2 = STOciS
			STOxi2 = STOxiS
		elseif(j >= 1 .and. j <= 3) then
			STOci2 = STOciP
			STOxi2 = STOxiP
		else
			call print_error_message(-1, "wrong j in kinetic_integral_slater")
		end if
		res = 0.0_dp

		do r=1,gaussians_per_slater
		do s=1,gaussians_per_slater

			res = res + kinetic_integral(vecA,vecB, STOxi1(r),STOxi2(s), i,j) * STOci1(r)*STOci2(s)

		end do
		end do

 	end function kinetic_integral_slater

 	real(dp) function coulombic_integral_slater(vecA,vecB,vecCoul,i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecCoul
		integer(i4b), intent(in)				:: i,j

		integer(i4b)							:: r,s
		real(dp),dimension(gaussians_per_slater) :: STOci1, STOci2, STOxi1, STOxi2

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiS
		else
			call print_error_message(-1, "wrong i in coulombic_integral_slater")
		end if
		if(j == 0) then
			STOci2 = STOciS
			STOxi2 = STOxiS
		elseif(j >= 1 .and. j <= 3) then
			STOci2 = STOciP
			STOxi2 = STOxiP
		else
			call print_error_message(-1, "wrong j in coulombic_integral_slater")
		end if

		res = 0.0_dp

		do r=1,gaussians_per_slater
		do s=1,gaussians_per_slater

			res = res + coulombic_integral(vecA,vecB,vecCoul, STOxi1(r),STOxi2(s), i,j) * STOci1(r)*STOci2(s)

		end do
		end do

 	end function coulombic_integral_slater

 	real(dp) function overlap_slater(vecA,vecB,i,j) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		integer(i4b), intent(in)				:: i,j

		integer(i4b)							:: r,s
		real(dp), dimension(gaussians_per_slater) :: STOci1, STOci2, STOxi1, STOxi2

		res = 0.0_dp

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiP
		else
			call print_error_message(-1, "wrong i in overlap_slater")
		end if
		if(j == 0) then
			STOci2 = STOciS
			STOxi2 = STOxiS
		elseif(j >= 1 .and. j <= 3) then
			STOci2 = STOciP
			STOxi2 = STOxiP
		else
			call print_error_message(-1, "wrong i in overlap_slater")
		end if


		do r=1,gaussians_per_slater
		do s=1,gaussians_per_slater

			res = res + overlap(vecA,vecB, STOxi1(r),STOxi2(s),i,j) * STOci1(r)*STOci2(s)

		end do
		end do

 	end function overlap_slater

 	real(dp) function density_slater(vecA,vecB,i) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		integer(i4b), intent(in)				:: i

		integer(i4b)							:: r
		real(dp), dimension(gaussians_per_slater) :: STOci1,STOxi1

		res = 0.0_dp

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiP
		else
			call print_error_message(-1, "wrong i in overlap_slater")
		end if

		do r=1,gaussians_per_slater

			res = res + density(vecA,vecB, STOxi1(r),i) * STOci1(r)

		end do

 	end function density_slater

 	real(dp) function dipole_xi_slater(vecA,vecB,i,j,dipole_index) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		integer(i4b), intent(in)				:: i,j,dipole_index

		integer(i4b)							:: r,s
		real(dp), dimension(gaussians_per_slater) :: STOci1, STOci2,STOxi1, STOxi2

		res = 0.0_dp

		if(i == 0) then
			STOci1 = STOciS
			STOxi1 = STOxiS
		elseif(i >= 1 .and. i <= 3) then
			STOci1 = STOciP
			STOxi1 = STOxiP
		else
			call print_error_message(-1, "wrong i in overlap_slater")
		end if
		if(j == 0) then
			STOci2 = STOciS
			STOxi2 = STOxiS
		elseif(j >= 1 .and. j <= 3) then
			STOci2 = STOciP
			STOxi2 = STOxiP
		else
			call print_error_message(-1, "wrong i in overlap_slater")
		end if


		do r=1,gaussians_per_slater
		do s=1,gaussians_per_slater

			res = res + dipole_xi(vecA,vecB, STOxi1(r),STOxi2(s),i,j,dipole_index) * STOci1(r)*STOci2(s)

		end do
		end do

 	end function dipole_xi_slater

 	real(dp) function four_center_integral_slater_directed(vecA,vecB,vecC,vecD,&
				i,j,k,l, nA,nB,nC,nD) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecC, vecD
		real(dp), dimension(3), intent(in) 	:: nA, nB, nC, nD
		logical, intent(in)					:: i,j,k,l

		integer(i4b)							:: r,s,t,u

		res = 0.0_dp

		if((.not. i) .or. (.not. j)) then
			call print_error_message(-1,"not implemented yet! (four_center_integral_slater_directed)")
		end if

		do r=1,3
		do s=1,3
		do t=1,3
		do u=1,3
			! we don't calculate too small contributions
			if(nA(r)*nB(s)*nC(t)*nD(u) <= 1e-6) then
				cycle
			end if

			res = res + four_center_integral_slater(vecA,vecB,vecC,vecD,r,s,t,u)*nA(r)&
					*nB(s)*nC(t)*nD(u)
		end do
		end do
		end do
		end do

	end function four_center_integral_slater_directed

 	real(dp) function kinetic_integral_slater_directed(vecA,vecB,i,j,nA,nB) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), dimension(3), intent(in) 	:: nA, nB
		logical, intent(in)					:: i,j

		integer(i4b)							:: r,s

		res = 0.0_dp

		if(size(nA) /= 3 .or. size(nB) /= 3 .or. size(vecA) /= 3 .or. size(vecB) /= 3) then
			call print_error_message(-1, "dimension_error in kinetic_integral_slater_directed")
		end if

		if((.not. i) .or. (.not. j)) then
			call print_error_message(-1,"not implemented yet! (kinetic_integral_slater_directed)")
		end if

		do r=1,3
		do s=1,3
			res = res + kinetic_integral_slater(vecA,vecB,r,s)*nA(r)*nB(s)
		end do
		end do

	end function kinetic_integral_slater_directed

 	real(dp) function coulombic_integral_slater_directed(vecA,vecB,vecCoul,&
				i,j, nA,nB) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB, vecCoul
		real(dp), dimension(3), intent(in) 	:: nA, nB
		logical, intent(in)					:: i,j

		integer(i4b)							:: r,s

		res = 0.0_dp

		if(size(nA) /= 3 .or. size(nB) /= 3 .or. size(vecA) /= 3 .or. size(vecB) /= 3 .or. size(vecCoul) /= 3) then
			call print_error_message(-1, "dimension_error in coulombic_integral_slater_directed")
		end if

		if((.not. i) .or. (.not. j)) then
			call print_error_message(-1,"not implemented yet! (coulombic_integral_slater_directed)")
		end if

		do r=1,3
		do s=1,3
			res = res + coulombic_integral_slater(vecA,vecB,vecCoul,r,s)*nA(r)*nB(s)
		end do
		end do

	end function coulombic_integral_slater_directed

 	real(dp) function overlap_slater_directed(vecA,vecB,i,j,nA,nB) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), dimension(3), intent(in) 	:: nA, nB
		logical, intent(in)					:: i,j

		integer(i4b)							:: r,s

		res = 0.0_dp

		if(size(nA) /= 3 .or. size(nB) /= 3 .or. size(vecA) /= 3 .or. size(vecB) /= 3) then
			call print_error_message(-1, "dimension_error in overlap_slater_directed")
		end if

		if((.not. i) .or. (.not. j)) then
			call print_error_message(-1,"not implemented yet! (overlap_slater_directed)")
		end if

		do r=1,3
		do s=1,3
			res = res + overlap_slater(vecA,vecB,r,s)*nA(r)*nB(s)
		end do
		end do

	end function overlap_slater_directed

 	real(dp) function dipole_xi_slater_directed(vecA,vecB,i,j,nA,nB,dipole_index) result(res)
		real(dp), dimension(3), intent(in) 	:: vecA, vecB
		real(dp), dimension(3), intent(in) 	:: nA, nB
		logical, intent(in)					:: i,j

		integer(i4b)							:: r,s,dipole_index

		if((.not. i) .or. (.not. j)) then
			call print_error_message(-1,"not implemented yet! (dipole_xi_slater_directed)")
		end if

		res = 0.0_dp

		do r=1,3
		do s=1,3
			res = res + dipole_xi_slater(vecA,vecB,r,s,dipole_index)*nA(r)*nB(s)
		end do
		end do

	end function dipole_xi_slater_directed

	subroutine fill_four_center_tensor_and_H_S_slater_directed(four_center_tensor,H,S,vec,n,Z,include_four_center_integrals)
		real(dp), dimension(:,:,:,:), intent(out) :: four_center_tensor
		real(dp), dimension(:,:), intent(out)	:: H,S
		real(dp), dimension(:,:), intent(in)	:: vec,n
		real(dp), dimension(:), intent(in)	:: Z
		logical, intent(in)					:: include_four_center_integrals

		real(dp)	  :: shadowing
		integer(i4b) :: a,b,c,d, mu,nu,la,si

		if(.not.(size(four_center_tensor,1) == size(four_center_tensor,2) .and. &
			size(four_center_tensor,2) == size(four_center_tensor,3) .and. &
			size(four_center_tensor,3) == size(four_center_tensor,4) .and. &
			size(four_center_tensor,1) == size(H,1) .and. size(H,1) == size(H,2) .and. &
			size(S,1) == size(S,2) .and. size(H,1) == size(S,2) .and. &
			size(vec,1) == 3 .and. size(vec,2) == size(Z) )) then

			call print_error_message(-1, "wrong dimensions in fill_four_center_tensor_and_H()")
		end if

		write(*,*) 'fill_intermolecular_four_center_tensor_and_H_slater_directed'

		four_center_tensor = 0.0_dp

		H = 0.0_dp
		S = 0.0_dp

		do mu=1,size(S,1)
		do nu=1,size(S,2)

			S(mu,nu) = overlap_slater_directed(vec(:,mu),vec(:,nu),.true.,.true.,n(:,mu),n(:,nu))

		end do
		end do

		do mu=1,size(H,1)
		do nu=1,size(H,2)

			H(mu,nu) = H(mu,nu) + kinetic_integral_slater_directed(vec(:,mu),vec(:,nu),.true.,.true.,n(:,mu),n(:,nu))

			do a=1,size(H,1)
				! both on the same nucleus => nucleus not perfectly shadowed
!				if(	dot_product(vec(:,mu)-vec(:,a),vec(:,mu)-vec(:,a)) <= 1e-6 .and. &
!				 	dot_product(vec(:,nu)-vec(:,a),vec(:,nu)-vec(:,a)) <= 1e-6) then
!					shadowing = shadow_factor
!				elseif(	dot_product(vec(:,mu)-vec(:,a),vec(:,mu)-vec(:,a)) > 1e-6 .and. &
!				 	dot_product(vec(:,nu)-vec(:,a),vec(:,nu)-vec(:,a)) > 1e-6) then
!
!				 	shadowing = 1.0_dp
!				else
!					shadowing = sqrt(shadow_factor)
!				end if
				shadowing = abs(S(mu,nu))*shadow_factor + (1.0_dp - abs(S(mu,nu)))*1

				H(mu,nu) = H(mu,nu) - Z(a)*shadowing*coulombic_integral_slater_directed&
											(vec(:,mu),vec(:,nu),vec(:,a),.true.,.true.,n(:,mu),n(:,nu))
			end do

		if(include_four_center_integrals) then

		do la=1,size(four_center_tensor,3)
		do si=1,size(four_center_tensor,4)

			!write(*,*) 'NEZOBRAZOVAT',mu,nu,la,si

			four_center_tensor(mu,nu,la,si) = four_center_integral_slater_directed(&
				vec(:,mu),vec(:,nu),vec(:,la),vec(:,si),.true.,.true.,.true.,.true., &
				n(:,mu),n(:,nu),n(:,la),n(:,si) )

		end do
		end do

		end if

		end do
		end do

	end subroutine fill_four_center_tensor_and_H_S_slater_directed

	subroutine fill_four_center_tensor_and_H_S_slater(four_center_tensor,H,S,vec,n,Z,T,V1,V2,slater_index)
		real(dp), dimension(:,:,:,:), intent(out) :: four_center_tensor
		real(dp), dimension(:,:), intent(out)	:: H,S,  T,V1,V2
		real(dp), dimension(:,:), intent(in)	:: vec,n
		real(dp), dimension(:), intent(in)	:: Z
		integer(i4b), intent(in)	:: slater_index

		real(dp)	  :: shadowing
		integer(i4b) :: a,b,c,d, mu,nu,la,si

		if(.not.(size(four_center_tensor,1) == size(four_center_tensor,2) .and. &
			size(four_center_tensor,2) == size(four_center_tensor,3) .and. &
			size(four_center_tensor,3) == size(four_center_tensor,4) .and. &
			size(four_center_tensor,1) == size(H,1).and. size(H,1) == size(H,2)) ) then

			call print_error_message(-1, "wrong dimensions in fill_four_center_tensor_and_H()")
		end if

		four_center_tensor = 0.0_dp

		H = 0.0_dp

		T = 0.0_dp
		V1 = 0.0_dp
		V2 = 0.0_dp

		S = 0.0_dp

		do mu=1,size(S,1)
		do nu=1,size(S,2)

			S(mu,nu) = overlap_slater(vec(:,mu),vec(:,nu),slater_index,slater_index)

		end do
		end do

		do mu=1,size(four_center_tensor,1)
		do nu=1,size(four_center_tensor,2)

			! kinetic_integral_slater gives directly the kinetic energy of particular orbital
			H(mu,nu) = H(mu,nu) + kinetic_integral_slater(vec(:,mu),vec(:,nu),slater_index,slater_index)
			T(mu,nu) = T(mu,nu) + kinetic_integral_slater(vec(:,mu),vec(:,nu),slater_index,slater_index)

			do a=1,size(H,1)
				! both on the same nucleus => nucleus not perfectly shadowed
				if(	dot_product(vec(:,mu)-vec(:,a),vec(:,mu)-vec(:,a)) <= 1e-6 .and. &
				 	dot_product(vec(:,nu)-vec(:,a),vec(:,nu)-vec(:,a)) <= 1e-6) then
					shadowing = shadow_factor
				elseif(	dot_product(vec(:,mu)-vec(:,a),vec(:,mu)-vec(:,a)) > 1e-6 .and. &
				 	dot_product(vec(:,nu)-vec(:,a),vec(:,nu)-vec(:,a)) > 1e-6) then

				 	shadowing = 1.0_dp
				else
					shadowing = sqrt(shadow_factor)
				end if

				H(mu,nu) = H(mu,nu) - Z(a)*shadowing*coulombic_integral_slater(vec(:,mu),vec(:,nu),vec(:,a),slater_index,slater_index)

				if(a==1)then
					V1(mu,nu) = V1(mu,nu) - Z(a)*shadowing*coulombic_integral_slater(vec(:,mu),vec(:,nu),vec(:,a),slater_index,slater_index)
				elseif(a==2)then
					V2(mu,nu) = V2(mu,nu) - Z(a)*shadowing*coulombic_integral_slater(vec(:,mu),vec(:,nu),vec(:,a),slater_index,slater_index)
				endif
			end do

		do la=1,size(four_center_tensor,3)
		do si=1,size(four_center_tensor,4)

			four_center_tensor(mu,nu,la,si) = four_center_integral_slater(&
				vec(:,mu),vec(:,nu),vec(:,la),vec(:,si),slater_index,slater_index,slater_index,slater_index)

		end do
		end do
		end do
		end do

	end subroutine fill_four_center_tensor_and_H_S_slater

	subroutine fill_F_matrix(Fa,Fb,Pa,Pb,four_center_tensor,H)
		real(dp), dimension(:,:), intent(out)	:: Fa,Fb
		real(dp), dimension(:,:), intent(in)	:: H,Pa,Pb
		real(dp), dimension(:,:,:,:), intent(in) :: four_center_tensor


		integer(i4b) :: a, mu,nu,la,si

		Fa = 0.0_dp
		Fb = 0.0_dp

		do mu=1,size(Fa,1)
		do nu=1,size(Fa,2)

			Fa(mu,nu) = H(mu,nu)
			Fb(mu,nu) = H(mu,nu)

		do la=1,size(Fa,1)
		do si=1,size(Fa,2)

			Fa(mu,nu) = Fa(mu,nu) + (Pa(la,si)+Pb(la,si))*four_center_tensor(mu,nu,la,si)
			Fb(mu,nu) = Fb(mu,nu) + (Pa(la,si)+Pb(la,si))*four_center_tensor(mu,nu,la,si)

			Fa(mu,nu) = Fa(mu,nu) - Pa(la,si)*four_center_tensor(mu,si,la,nu)
			Fb(mu,nu) = Fb(mu,nu) - Pb(la,si)*four_center_tensor(mu,si,la,nu)

		end do
		end do
		end do
		end do

!		write(*,*) '  fill_F H',H
!		write(*,*) '  fill_F Pa',Pa
!		write(*,*) '  fill_F Pb',Pb
!		write(*,*) '  fill_F Fa',Fa
!		write(*,*) '  fill_F Fb',Fb

	end subroutine fill_F_matrix

	subroutine fill_P_matrix(Pa,Pb,ca,cb,excMol)
		real(dp), dimension(:,:), intent(in)	:: ca,cb
		logical, dimension(:,:), intent(in)	:: excMol
		real(dp), dimension(:,:), intent(out)	:: Pa,Pb

		integer(i4b) :: i,mu,nu

		Pa = 0.0_dp
		Pb = 0.0_dp

		do mu=1,size(ca,1)
		do nu=1,size(ca,1)

		do i=1,size(ca,1)
			! conjg should be here if c are complex
			if(excMol(1,i)) then
				Pa(mu,nu) = Pa(mu,nu) + ca(mu,i)*ca(nu,i)
			end if

			if(excMol(2,i)) then
				Pb(mu,nu) = Pb(mu,nu) + cb(mu,i)*cb(nu,i)
			end if

		end do
		end do
		end do

!		write(*,*) '  excMol', excMol
!		write(*,*) '  ca', ca
!		write(*,*) '  cb', cb
!		write(*,*) '  Pa', Pa
!		write(*,*) '  Pb', Pb

	end subroutine fill_P_matrix

	subroutine fill_HOMO_state(excMol, electrons)
		logical, dimension(:,:), intent(out)	:: excMol
		integer(i4b), intent(in) :: electrons

		integer(i4b)	::	i, e, j

		if(size(excMol,1) /= 2 .or. electrons > 2*size(excMol,2)) then
			call print_error_message(-1, "error in fill_HOMO_state")
		end if

		excMol = .false.
		e = 0

		do i = 1, size(excMol,2)
		do j = 1, 2
			if(e >= electrons) then
				cycle
			end if

			excMol(j,i) = .true.
			e = e + 1
		end do
		end do
	end subroutine fill_HOMO_state

	subroutine fill_LUMO_state(excMol, electrons)
		logical, dimension(:,:), intent(out)	:: excMol
		integer(i4b), intent(in) :: electrons

		integer(i4b)	::	i, j

		if(size(excMol,1) /= 2 .or. electrons > 2*size(excMol,2)) then
			call print_error_message(-1, "error in fill_LUMO_state")
		end if

		call fill_HOMO_state(excMol, electrons)

		do i = size(excMol,2) - 1, 1, -1
		do j = 1, 2 ! HACK
			if(excMol(j,i) .and. .not. excMol(j,i+1)) then
				excMol(j,i+1) = .true.
				excMol(j,i) = .false.
				return
			end if
		end do
		end do
	end subroutine fill_LUMO_state

	subroutine solve_Roothan_equations(vec,n,Z,S,H,four_center_tensor, ca,cb,Fa,Fb,Ea,Eb)
		real(dp), dimension(:,:), intent(in)			:: vec, n, S, H
		real(dp), dimension(:,:,:,:), intent(in)		:: four_center_tensor
		real(dp), dimension(:), intent(in)			:: Z
		real(dp), dimension(:,:), intent(out)			:: ca, cb, Fa, Fb
		real(dp), dimension(:), intent(out)			:: Ea, Eb


		real(dp), parameter							:: lambda = 1.0_dp
		real(dp), dimension(size(vec,2),size(vec,2))	:: Wa, Wb, Pa, Pb, Fa_drive, Fb_drive
		logical, dimension(2,size(Z))					:: excMolGround, excMolExc

		real(dp), dimension(size(vec,2),size(vec,2))	:: ca_prev, cb_prev, ca_pprev, cb_pprev
		real(dp)										:: E_full, E_full_prev, E_full_pprev

		integer(i4b)	:: i,j,iteration,interation_with_the_same_energy

		if(.not.(size(Z) == size(vec,2).and. size(Z) == size(n,2).and. &
				size(vec,1) == size(n,1) .and. size(vec,1) == 3)) then

				call print_error_message(-1, 'dimension error in solve_Rootham_equations')
		end if

!		write(*,*) 'tens:', four_center_tensor
!		write(*,*) 'H:', H
!		write(*,*) 'S:', S

		call fill_HOMO_state(excMolGround, size(Z))

		excMolExc = .false.
		excMolExc(1,2) = .true.
		excMolExc(2,2) = .true.

		Fa = 0.0_dp
		Fb = 0.0_dp
		Fa_drive = 0.0_dp
		Fb_drive = 0.0_dp
		ca = 0.0_dp
		cb = 0.0_dp
		Ea = 0.0_dp
		Eb = 0.0_dp

		! orthogonal initial condition favours unrestricted solution
		ca(:,1) = (/ 1.0, 0.0 /)
		cb(:,1) = (/ 0.0, 1.0 /)
		ca_prev = 1e+6_dp
		cb_prev = 1e+6_dp
		ca_pprev = 1e+6_dp
		cb_pprev = 1e+6_dp
		E_full_prev = 1e9_dp
		E_full_pprev = 1e9_dp

		interation_with_the_same_energy = 0

		do i=1, 500
			call fill_F_matrix(Fa_drive,Fb_drive,Pa,Pb,four_center_tensor,H)
			Fa = Fa*(1-lambda) + Fa_drive*lambda
			Fb = Fb*(1-lambda) + Fb_drive*lambda

			call spec_generalized(Fa,S,ca,Wa)
			call spec_generalized(Fb,S,cb,Wb)

			E_full = full_energy(vec,n,Z,four_center_tensor,H,Pa,Pb)

			if(.false.) then
				write(*,*) '*************>  ca      ',ca
				write(*,*) '*************>  cb      ',cb
				write(*,*) '+++++++++++++>  Pa      ',Pa
				write(*,*) '+++++++++++++>  Pb      ',Pb
				write(*,*) '------------->  Fa      ',Fa
				write(*,*) '------------->  Fb      ',Fb
				write(*,*) '#############>  Fa_drive',Fa_drive
				write(*,*) '/////////////>  Wa',Wa
				write(*,*) '\\\\\\\\\\\\\>  E', E_full, E_full_prev
			end if

			if(abs(E_full-E_full_prev) <= 1e-4_dp .and. &
				maxval(abs(Fa-Fa_drive)) <= 1e-2_dp) then
				interation_with_the_same_energy = interation_with_the_same_energy + 1
			else
				interation_with_the_same_energy = 0
			end if

			if(interation_with_the_same_energy >= 3) then
				exit
			end if

			! detection of oscillations between 2 points
			if(maxval(abs(ca-ca_prev))/maxval(abs(ca)+1e-21) >= 1e-2 .and. &
				maxval(abs(ca_pprev-ca))/maxval(abs(ca)+1e-21) < 1e-2) then
				write(*,*) 'OSCILLATION DETECTED!'
				call random_number(ca)
				call random_number(cb)

				write(*,*) '>>>>>>>>>>>>>>  ca      ',ca
				write(*,*) '>>>>>>>>>>>>>>  cb      ',cb
				ca_prev = 1e6_dp
				cb_prev = 1e6_dp
				ca_pprev = 1e6_dp
				cb_pprev = 1e6_dp
				E_full_prev = 1e9_dp
				E_full_pprev = 1e9_dp

				call fill_P_matrix(Pa,Pb,ca,cb,excMolGround)
				cycle
			end if

			call fill_P_matrix(Pa,Pb,ca,cb,excMolGround)

			ca_pprev = ca_prev
			cb_pprev = cb_prev
			E_full_pprev = E_full_prev
			ca_prev = ca
			cb_prev = cb
			E_full_prev = E_full

		end do

!		call fill_P_matrix(Pa,Pb,ca,cb,excMolGround)
		call fill_F_matrix(Fa,Fb,Pa,Pb,four_center_tensor,H)

!			write(*,*) '  excMol',excMolGround
!			write(*,*)
!			write(*,*) '  S',S
!			write(*,*)
!			write(*,*) '  Pa',Pa
!			write(*,*)
!			write(*,*) '  Pb',Pb
!			write(*,*)
!			write(*,*) '  Fa',Fa
!			write(*,*)
!			write(*,*) '  Fb',Fb
!			write(*,*)
!			write(*,*) '  ca',ca
!			write(*,*)
!			write(*,*) '  cb',cb
!			write(*,*)
!			write(*,*) '  Wa',Wa
!			write(*,*)
!			write(*,*) '  Wb',Wb
!			write(*,*)

!			write(*,*) '  Hm',H
!			write(*,*)
!			write(*,*) '  Tens',four_center_tensor
!			write(*,*)

!			write(*,*) ' S_munu P_munu = ', trace(matmul(Pa+Pb,S))
!			write(*,*)

!			write(*,*) ' Eg2',full_energy(vec,n,Z,four_center_tensor,H,Pa,Pb)
			call fill_P_matrix(Pa,Pb,ca,cb,excMolExc)
!			write(*,*) ' Ee2',full_energy(vec,n,Z,four_center_tensor,H,Pa,Pb)
!			write(*,*) '          Pa',Pa

!			write(*,*) ' test_vlastnich_vektoru',matmul(Fa,ca(:,1)),Wa(1,1)*matmul(S,ca(:,1))
!			write(*,*) ' test_vlastnich_vektoru',matmul(Fa,ca(:,2)),Wa(2,2)*matmul(S,ca(:,2))
!			write(*,*) ' test_normalizace',dot_product(ca(:,1),matmul(S,ca(:,1)))
!			write(*,*) ' test_normalizace',dot_product(ca(:,2),matmul(S,ca(:,2)))

			! If convergence not reached, set values to some obvious nonsense
			if(.not. maxval(abs(matmul(Fa,ca(:,1))-Wa(1,1)*matmul(S,ca(:,1)))) < 1e-1) then
				write(*,*) 'CONVERGENCE NOT REACHED'
				ca = 1e6_dp
				cb = 1e6_dp
				Wa = 1e6_dp
				Wb = 1e6_dp
				Fa = 1e6_dp
				Fb = 1e6_dp
			end if

		do i=1,size(Wa,1)
			Ea(i) = Wa(i,i)
			Eb(i) = Wb(i,i)
		end do

	end subroutine solve_Roothan_equations

	real(dp) function full_energy(vec,n,Z,four_center_tensor,H,Pa,Pb) result(E)
		real(dp), dimension(:,:,:,:), intent(in)		:: four_center_tensor
		real(dp), dimension(:,:), intent(in)			:: H,Pa,Pb
		real(dp), dimension(:,:), intent(in)			:: vec, n
		real(dp), dimension(:), intent(in)			:: Z


		real(dp), dimension(size(H,1),size(H,2))	:: Fa,Fb
		integer(i4b)	:: a,b,mu

		! Calculates only energy of introduced electrons, contributions electron-electron, electron-nucleus.
		! Energy of type nucleus-nucleus is not accounted.

		call fill_F_matrix(Fa,Fb,Pa,Pb,four_center_tensor,H)

!		write(*,*)
!		write(*,*) 'FaR',Fa
!		write(*,*) '(FaR+H)/2',(H+Fa)/2

		E = 0.0_dp
		E = E + trace(matmul(Pa,H+Fa))/2
		E = E + trace(matmul(Pb,H+Fb))/2

		do a=1,size(Z)
		do b=a+1,size(Z)

			E = E + Z(a)*Z(b)/sqrt(dot_product(vec(:,a)-vec(:,b),vec(:,a)-vec(:,b)))
		end do
		end do

	end function full_energy

	real(dp) function full_energy_Huckel(vec,n,Z,H,Pa,Pb) result(E)
		real(dp), dimension(:,:), intent(in)			:: H,Pa,Pb
		real(dp), dimension(:,:), intent(in)			:: vec, n
		real(dp), dimension(:), intent(in)			:: Z


		real(dp), dimension(size(H,1),size(H,2))	:: Fa,Fb
		integer(i4b)	:: a,b,mu

		! Calculates only energy of introduced electrons, contributions electron-electron, electron-nucleus.
		! Energy of type nucleus-nucleus is not accounted.

		E = 0.0_dp
		E = E + trace(matmul(Pa+Pb,H))

	end function full_energy_Huckel

	subroutine calculate_Huckel_method(vec,n,Z,S, c,F,E, parameter_A,parameter_B)
		real(dp), dimension(:,:), intent(in)			:: vec, n, S
		real(dp), dimension(:), intent(in)			:: Z
		real(dp)										:: parameter_A, parameter_B
		real(dp), dimension(:,:), intent(out)			:: c, F
		real(dp), dimension(:), intent(out)			:: E

		real(dp), parameter							:: WH_constant = 1.75_dp
		real(dp), parameter							:: C_core = -11.4_dp/hartree_in_eV
		real(dp), dimension(size(S,2),size(S,2))		:: W

		integer(i4b)	:: i,j, mu,nu

		if(.not.(size(Z) == size(vec,2).and. size(Z) == size(n,2).and. &
				size(vec,1) == size(n,1) .and. size(vec,1) == 3)) then

				call print_error_message(-1, 'dimension error in calculate_Huckel_method')
		end if

		c = 0.0_dp
		F = 0.0_dp
		E = 0.0_dp

		do mu=1,size(S,1)
		do nu=1,size(S,2)
			if(mu == nu) then
				F(mu,nu) = parameter_A
			else
				F(mu,nu) = S(mu,nu)*parameter_B
			endif
		end do
		end do

		call spec_generalized(F,S,c,W)

		do mu=1,size(E)
			E(mu) = W(mu,mu)
		end do

	end subroutine calculate_Huckel_method


	subroutine excitation_difference_code(excState1, excState2, first1, first_spin1, first2, first_spin2, second1, second_spin1, second2, second_spin2, n)
		logical, dimension(:,:), intent(in)   :: excState1,excState2

		integer(i4b), intent(out)	:: first1, first_spin1, first2, first_spin2, second1, second_spin1, second2, second_spin2, n
		integer(i4b)				:: i,j,k,l

		if(.not.(size(excState1,1) == 2 .and. size(excState1,2) == size(excState2,2) &
			.and. size(excState1,1) == 2)) then

			call print_error_message(-1,'Dimension error in excitation_difference_code')
		end if

		first1 = 0
		first_spin1 = 0
		first2 = 0
		first_spin2 = 0
		second1 = 0
		second_spin1 = 0
		second2 = 0
		second_spin2 = 0

		n = 0
		do i=1,2
		do j=1,size(excState1,2)
			if(.not. (excState1(i,j) .eqv. excState2(i,j))) then

				if(n < 2) then
					if(excState1(i,j)) then
						first1 = j
						first_spin1 = i
					else
						first2 = j
						first_spin2 = i
					end if
				else
					if(excState1(i,j)) then
						second1 = j
						second_spin1 = i
					else
						second2 = j
						second_spin2 = i
					end if
				end if

				n = n + 1
			end if
		end do
		end do

		n = n / 2

	end subroutine excitation_difference_code

	real(dp) function dipole_moment_between_states(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,&
							ca11, cb11, ca12, cb12, DD) result(Dip)
		real(dp), dimension(:,:), intent(in)	:: vec1,n1,ca11,cb11,ca12,cb12,DD
		real(dp), dimension(:), intent(in)	:: Z1,Z2
		real(dp), dimension(:,:), intent(in)	:: vec2,n2

		logical, dimension(:,:), intent(in)   :: excMol1_1,excMol1_2

		! the same states (0), differing in one(i), differ in two(10^6i + j), differ in more(2*10^6)
		integer(i4b)				:: firstDifIndMol1_1,firstDifIndMol1_2
		integer(i4b)				:: secondDifIndMol1_1,secondDifIndMol1_2
		integer(i4b)				:: firstDifSpinMol1_1,firstDifSpinMol1_2
		integer(i4b)				:: secondDifSpinMol1_1,secondDifSpinMol1_2

		integer(i4b)				:: i,j,k,l,n,mu,nu,a,b,la,si,c,d,s,differing_excitationsMol1


		call excitation_difference_code(excMol1_1,excMol1_2,&
			firstDifIndMol1_1,firstDifSpinMol1_1,firstDifIndMol1_2,firstDifSpinMol1_2,&
			secondDifIndMol1_1,secondDifSpinMol1_1,secondDifIndMol1_1,secondDifSpinMol1_1,&
			differing_excitationsMol1)

		Dip = 0.0_dp

		! Simple application of Slater-Condon rule

		if(differing_excitationsMol1 < 2) then
			do mu=1,size(Z1)
			do nu=1,size(Z1)

			if(differing_excitationsMol1 == 0) then
			! no differing spinorbitals
			do i=1,size(Z1)
				! we may not want H1 there, because it is not interaction energy
				if(excMol1_1(1,i)) then
					Dip = Dip + ca11(mu,i)*ca11(nu,i)*DD(mu,nu)
				end if
				if(excMol1_1(2,i)) then
					Dip = Dip + cb11(mu,i)*cb11(nu,i)*DD(mu,nu)
				end if
			end do
			elseif(differing_excitationsMol1 == 1) then
				! one differing spinorbital
				if(firstDifSpinMol1_1 == 1 .and. firstDifSpinMol1_2 == 1) then
					Dip = Dip + ca11(mu,firstDifIndMol1_1)*ca12(nu,firstDifIndMol1_2)*DD(mu,nu)
				elseif(firstDifSpinMol1_1 == 1 .and. firstDifSpinMol1_2 == 2) then
					Dip = Dip + ca11(mu,firstDifIndMol1_1)*cb12(nu,firstDifIndMol1_2)*DD(mu,nu)
				elseif(firstDifSpinMol1_1 == 2 .and. firstDifSpinMol1_2 == 1) then
					Dip = Dip + cb11(mu,firstDifIndMol1_1)*ca12(nu,firstDifIndMol1_2)*DD(mu,nu)
				else
					Dip = Dip + cb11(mu,firstDifIndMol1_1)*cb12(nu,firstDifIndMol1_2)*DD(mu,nu)
				end if
			end if

			end do
			end do
		end if

	end function dipole_moment_between_states

	real(dp) function interaction_energy_inside_one_molecule(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,&
							four_center_tensor1, H1, ca11, cb11, ca12, cb12, HMol1_2) result(Eint)
		real(dp), dimension(:,:), intent(in)	:: vec1,n1,ca11,cb11,ca12,cb12,H1
		real(dp), dimension(:), intent(in)	:: Z1,Z2
		real(dp), dimension(:,:), intent(in)	:: vec2,n2

		real(dp), dimension(:,:,:,:), intent(in)	:: four_center_tensor1

		real(dp), dimension(:,:), intent(in) 	:: HMol1_2

		logical, dimension(:,:), intent(in)   :: excMol1_1,excMol1_2

		! the same states (0), differing in one(i), differ in two(10^6i + j), differ in more(2*10^6)
		integer(i4b)				:: firstDifIndMol1_1,firstDifIndMol1_2
		integer(i4b)				:: secondDifIndMol1_1,secondDifIndMol1_2
		integer(i4b)				:: firstDifSpinMol1_1,firstDifSpinMol1_2
		integer(i4b)				:: secondDifSpinMol1_1,secondDifSpinMol1_2

		integer(i4b)				:: i,j,k,l,n,mu,nu,a,b,la,si,c,d,s,differing_excitationsMol1

		real(dp) :: tmp, tmp2


		call excitation_difference_code(excMol1_1,excMol1_2,&
			firstDifIndMol1_1,firstDifSpinMol1_1,firstDifIndMol1_2,firstDifSpinMol1_2,&
			secondDifIndMol1_1,secondDifSpinMol1_1,secondDifIndMol1_1,secondDifSpinMol1_1,&
			differing_excitationsMol1)

		Eint = 0.0_dp

!		HMol1_2 = 0.0_dp

		! Terms of H_core consist of u_i(1)|f|u_i(1) terms. Since there are no
		! exchange effects between molecules and Slater orbitals have completely
		! different coordinates, only off-diagonal elements of H

		if(differing_excitationsMol1 < 2) then
			do mu=1,size(Z1)
			do nu=1,size(Z1)

!			do a=1,size(Z2)
!				HMol1_2(mu,nu) = HMol1_2(mu,nu) - Z2(a)*coulombic_integral_slater_directed&
!							(vec1(:,mu),vec1(:,nu),vec2(:,a),.true.,.true., n1(:,mu),n1(:,nu))
!
!			end do

			if(differing_excitationsMol1 == 0) then
			! no differing spinorbitals
			do i=1,size(Z1)
				! we may not want H1 there, because it is not interaction energy
				if(excMol1_1(1,i)) then
					Eint = Eint + ca11(mu,i)*ca11(nu,i)*(H1(mu,nu) + HMol1_2(mu,nu))
				end if
				if(excMol1_1(2,i)) then
					Eint = Eint + cb11(mu,i)*cb11(nu,i)*(H1(mu,nu) + HMol1_2(mu,nu))
				end if
			end do
			elseif(differing_excitationsMol1 == 1) then
				! one differing spinorbital
				if(firstDifSpinMol1_1 == 1 .and. firstDifSpinMol1_2 == 1) then
					Eint = Eint + ca11(mu,firstDifIndMol1_1)*ca12(nu,firstDifIndMol1_2)*(H1(mu,nu) + HMol1_2(mu,nu))
				elseif(firstDifSpinMol1_1 == 1 .and. firstDifSpinMol1_2 == 2) then
					Eint = Eint + ca11(mu,firstDifIndMol1_1)*cb12(nu,firstDifIndMol1_2)*(H1(mu,nu) + HMol1_2(mu,nu))
				elseif(firstDifSpinMol1_1 == 2 .and. firstDifSpinMol1_2 == 1) then
					Eint = Eint + cb11(mu,firstDifIndMol1_1)*ca12(nu,firstDifIndMol1_2)*(H1(mu,nu) + HMol1_2(mu,nu))
				else
					Eint = Eint + cb11(mu,firstDifIndMol1_1)*cb12(nu,firstDifIndMol1_2)*(H1(mu,nu) + HMol1_2(mu,nu))
				end if
			end if

			end do
			end do
		end if

		! 2-electron contributions inside one molecule
		if(differing_excitationsMol1 < 3) then
			do mu=1,size(Z1)
			do nu=1,size(Z1)
			do la=1,size(Z1)
			do si=1,size(Z1)

			if(differing_excitationsMol1 == 0) then
			do i=1,size(Z1)-1
			do j=i+1,size(Z1)
				!no difference between c11 and c12 anyway
				if(excMol1_1(1,i) .and. excMol1_1(1,j)) then
					Eint = Eint + (	ca11(mu,i)*ca11(nu,j)*ca11(la,i)*ca11(si,j) - &
									ca11(mu,i)*ca11(nu,j)*ca11(si,i)*ca11(la,j))*&
						four_center_tensor1(mu,nu,la,si)
				end if
				if(excMol1_1(2,i) .and. excMol1_1(2,j)) then
					Eint = Eint + (	ca11(mu,i)*ca11(nu,j)*ca11(la,i)*ca11(si,j) - &
									ca11(mu,i)*ca11(nu,j)*ca11(si,i)*ca11(la,j))*&
						four_center_tensor1(mu,nu,la,si)
				end if
			end do
			end do
			elseif(differing_excitationsMol1 == 1) then
				! one differing spinorbital
				do j=1,size(Z1)
				do s=1,2
					tmp  =  four_center_tensor1(mu,nu,la,si)
					tmp2 = -four_center_tensor1(mu,nu,la,si)

					if(firstDifSpinMol1_1 == 1) then
						tmp  = tmp*ca11(mu,firstDifIndMol1_1)
						tmp2 = tmp*ca11(mu,firstDifIndMol1_1)
					else
						tmp  = tmp*cb12(mu,firstDifIndMol1_1)
						tmp2 = tmp*cb12(mu,firstDifIndMol1_1)
					end if
					if(firstDifSpinMol1_2 == 1) then
						tmp  = tmp*ca11(la,firstDifIndMol1_2)
						tmp2 = tmp*ca11(si,firstDifIndMol1_2)
					else
						tmp  = tmp*cb12(la,firstDifIndMol1_2)
						tmp2 = tmp*cb12(si,firstDifIndMol1_2)
					end if
					if(excMol1_2(s,j)) then
						if(s == 1) then
							tmp  = tmp*ca11(nu,j)*ca11(si,j)
							tmp2 = tmp*ca11(nu,j)*ca11(la,j)
						else
							tmp  = tmp*cb11(nu,j)*cb11(si,j)
							tmp2 = tmp*cb11(nu,j)*cb11(la,j)
						end if
					else
						cycle
					end if

					Eint = Eint + tmp + tmp2

				end do
				end do
			else
					tmp  =  four_center_tensor1(mu,nu,la,si)
					tmp2 = -four_center_tensor1(mu,nu,la,si)

					if(firstDifSpinMol1_1 == 1) then
						tmp  = tmp*ca11(mu,firstDifIndMol1_1)
						tmp2 = tmp*ca11(mu,firstDifIndMol1_1)
					else
						tmp  = tmp*cb11(mu,firstDifIndMol1_1)
						tmp2 = tmp*cb11(mu,firstDifIndMol1_1)
					end if
					if(firstDifSpinMol1_2 == 1) then
						tmp  = tmp*ca12(la,firstDifIndMol1_2)
						tmp2 = tmp*ca12(si,firstDifIndMol1_2)
					else
						tmp  = tmp*cb12(la,firstDifIndMol1_2)
						tmp2 = tmp*cb12(si,firstDifIndMol1_2)
					end if
					if(secondDifSpinMol1_1 == 1) then
						tmp  = tmp*ca11(nu,secondDifIndMol1_1)
						tmp2 = tmp*ca11(nu,secondDifIndMol1_1)
					else
						tmp  = tmp*cb11(nu,secondDifIndMol1_1)
						tmp2 = tmp*cb11(nu,secondDifIndMol1_1)
					end if
					if(secondDifSpinMol1_2 == 1) then
						tmp  = tmp*ca12(si,secondDifIndMol1_2)
						tmp2 = tmp*ca12(la,secondDifIndMol1_2)
					else
						tmp  = tmp*cb12(si,secondDifIndMol1_2)
						tmp2 = tmp*cb12(la,secondDifIndMol1_2)
					end if
					Eint = Eint + tmp + tmp2
			end if

			end do
			end do
			end do
			end do
		end if


	end function interaction_energy_inside_one_molecule

	real(dp) function interaction_energy(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,excMol2_1,excMol2_2, &
							four_center_tensor1, H1, four_center_tensor2, H2, &
							ca11, cb11, ca12, cb12, ca21, cb21, ca22, cb22, &
							four_center_tensor_int_12, Hint1_2, Hint2_1) result(Eint)
!							ca1, cb1, ca2, cb2) result(Eint)
		real(dp), dimension(:,:), intent(in)	:: vec1,n1,H1, Hint1_2!, ca1, cb1, ca2, cb2
		real(dp), dimension(:), intent(in)	:: Z1
		real(dp), dimension(:,:), intent(in)	:: vec2,n2,H2, Hint2_1
		real(dp), dimension(:,:), intent(in)	:: ca21,cb21,ca22,cb22,ca11,cb11,ca12,cb12
		real(dp), dimension(:), intent(in)	:: Z2

		real(dp), dimension(:,:,:,:), intent(in)	:: four_center_tensor1, four_center_tensor2,&
													   four_center_tensor_int_12

		! excMol1_1(spin, level)
		logical, dimension(:,:), intent(in)   :: excMol1_1,excMol1_2,excMol2_1,excMol2_2

		! the same states (0), differing in one(i), differ in two(10^6i + j), differ in more(2*10^6)
		integer(i4b)				:: firstDifIndMol1_1,firstDifIndMol1_2
		integer(i4b)				:: firstDifIndMol2_1,firstDifIndMol2_2
		integer(i4b)				:: secondDifIndMol1_1,secondDifIndMol1_2
		integer(i4b)				:: secondDifIndMol2_1,secondDifIndMol2_2
		integer(i4b)				:: firstDifSpinMol1_1,firstDifSpinMol1_2
		integer(i4b)				:: firstDifSpinMol2_1,firstDifSpinMol2_2
		integer(i4b)				:: secondDifSpinMol1_1,secondDifSpinMol1_2
		integer(i4b)				:: secondDifSpinMol2_1,secondDifSpinMol2_2
		integer(i4b)				:: i,j,k,l,n,mu,nu,a,b,la,si,c,d,s,s1,s2,differing_excitationsMol1,differing_excitationsMol2

!		real(dp), dimension(size(H1,1),size(H1,2)) 	:: HMol1_2
!		real(dp), dimension(size(H2,1),size(H2,2)) 	:: HMol2_1

		real(dp)	:: tmp

		if(	size(four_center_tensor_int_12,1) /= size(H1,1) .or. &
			size(four_center_tensor_int_12,2) /= size(H1,2) .or. &
			size(four_center_tensor_int_12,3) /= size(H2,1) .or. &
			size(four_center_tensor_int_12,4) /= size(H2,2) .or. &
			size(four_center_tensor1,1) /= size(H1,1) .or. &
			size(four_center_tensor1,2) /= size(H1,1) .or. &
			size(four_center_tensor1,3) /= size(H1,1) .or. &
			size(four_center_tensor1,4) /= size(H1,1) .or. &
			size(four_center_tensor2,1) /= size(H2,1) .or. &
			size(four_center_tensor2,2) /= size(H2,1) .or. &
			size(four_center_tensor2,3) /= size(H2,1) .or. &
			size(four_center_tensor2,4) /= size(H2,1) ) then

			call print_error_message(-1, 'dimension error in interaction_energy()')
		endif


		call excitation_difference_code(excMol1_1,excMol1_2,&
			firstDifIndMol1_1,firstDifSpinMol1_1,firstDifIndMol1_2,firstDifSpinMol1_2,&
			secondDifIndMol1_1,secondDifSpinMol1_1,secondDifIndMol1_1,secondDifSpinMol1_1,&
			differing_excitationsMol1)
		call excitation_difference_code(excMol2_1,excMol2_2,&
			firstDifIndMol2_1,firstDifSpinMol2_1,firstDifIndMol2_2,firstDifSpinMol2_2,&
			secondDifIndMol2_1,secondDifSpinMol2_1,secondDifIndMol2_1,secondDifSpinMol2_1,&
			differing_excitationsMol2)

		Eint = 0.0_dp

!		HMol1_2 = 0.0_dp
!		HMol2_1 = 0.0_dp

		! nuclei - nuclei - comes into game only if both molecules stay unchanged
		if(differing_excitationsMol1 == 0 .and. differing_excitationsMol2 == 0) then
			do i=1,size(Z1)
			do j=1,size(Z2)
				Eint = Eint + Z1(i)*Z2(j)/sqrt(dot_product(vec1(:,i),vec2(:,j)))
			end do
			end do
		end if


		! comes into game only if the other molecule stays unchanged
		if(differing_excitationsMol2 == 0) then
		Eint = Eint + interaction_energy_inside_one_molecule(&
				vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,four_center_tensor1, &
				H1, ca11, cb11, ca12, cb12, Hint1_2)
		end if
		if(differing_excitationsMol1 == 0) then
		Eint = Eint + interaction_energy_inside_one_molecule(&
				vec2,n2,Z2,vec1,n1,Z1, excMol2_1,excMol2_2,four_center_tensor2, &
				H2, ca21, cb21, ca22, cb22, Hint2_1)
		end if

!		write(*,*)
!		write(*,*) 'internalo energy consistency check:',Eint
		do mu=1,size(H1,1)
		do nu=1,size(H1,1)
		do c=1,size(H2,1)
		do d=1,size(H2,1)
		if(differing_excitationsMol1 == 0 .and. differing_excitationsMol2 == 0) then
			call print_warning_message("Interaction energy when one molecule's state doesn't change not tested!", -1)
			do i=1,size(H1,1)
			do s1=1,2
			do j=1,size(H2,1)
			do s2=1,2
				if((.not.(excMol1_1(s1,i))).or.(.not.(excMol2_1(s2,j))) ) then
					cycle
				end if

				tmp = four_center_tensor_int_12(mu,nu,c,d)
				if(s1 == 1) then
					tmp = tmp*ca11(mu,i)*ca11(nu,i)
				else
					tmp = tmp*cb11(mu,i)*cb11(nu,i)
				end if

				if(s2 == 1) then
					tmp = tmp*ca21(c,j)*ca21(d,j)
				else
					tmp = tmp*cb21(c,j)*cb21(d,j)
				end if

				Eint = Eint + tmp
			end do
			end do
			end do
			end do
		elseif(differing_excitationsMol1 == 1 .and. differing_excitationsMol2 == 0) then
			call print_warning_message("Interaction energy when one molecule's state doesn't change not tested!", -1)
			do i=1,size(H2,1)
			do s=1,2
				if(.not.(excMol2_1(s,i))) then
					cycle
				end if

				tmp = four_center_tensor_int_12(mu,nu,c,d)
				if(firstDifSpinMol1_1 == 1) then
					tmp = tmp*ca11(mu,firstDifIndMol1_1)
				else
					tmp = tmp*cb11(mu,firstDifIndMol1_1)
				end if
				if(firstDifSpinMol1_2 == 1) then
					tmp = tmp*ca12(nu,firstDifIndMol1_2)
				else
					tmp = tmp*cb12(nu,firstDifIndMol1_2)
				end if
				if(s == 1) then
					tmp = tmp*ca21(c,i)*ca21(d,i)
				else
					tmp = tmp*cb21(c,i)*cb21(d,i)
				end if

				Eint = Eint + tmp
			end do
			end do
		elseif(differing_excitationsMol1 == 0 .and. differing_excitationsMol2 == 1) then
			call print_warning_message("Interaction energy when one molecule's state doesn't change not tested!", -1)
			do i=1,size(H1,1)
			do s=1,2
				if(.not.(excMol1_1(s,i))) then
					cycle
				end if

				tmp = four_center_tensor_int_12(mu,nu,c,d)
				if(s == 1) then
					tmp = tmp*ca11(mu,i)*ca11(nu,i)
				else
					tmp = tmp*cb11(mu,i)*cb11(nu,i)
				end if
				if(firstDifSpinMol2_1 == 1) then
					tmp = tmp*ca21(c,firstDifIndMol2_1)
				else
					tmp = tmp*cb21(c,firstDifIndMol2_1)
				end if
				if(firstDifSpinMol2_2 == 1) then
					tmp = tmp*ca22(d,firstDifIndMol2_2)
				else
					tmp = tmp*cb22(d,firstDifIndMol2_2)
				end if

				Eint = Eint + tmp
			end do
			end do
		elseif(differing_excitationsMol1 == 1 .and. differing_excitationsMol2 == 1) then
			tmp = four_center_tensor_int_12(mu,nu,c,d)

			if(firstDifSpinMol1_1 == 1) then
				tmp = tmp*ca11(mu,firstDifIndMol1_1)
!				write(*,'(f12.6) ',advance='no') ca11(mu,firstDifIndMol1_1)
			else
				tmp = tmp*cb11(mu,firstDifIndMol1_1)
!				write(*,'(f12.6) ',advance='no') cb11(mu,firstDifIndMol1_1)
			end if
			if(firstDifSpinMol1_2 == 1) then
				tmp = tmp*ca12(nu,firstDifIndMol1_2)
!				write(*,'(f12.6) ',advance='no') ca12(nu,firstDifIndMol1_2)
			else
				tmp = tmp*cb12(nu,firstDifIndMol1_2)
!				write(*,'(f12.6) ',advance='no') cb12(nu,firstDifIndMol1_2)
			end if
			if(firstDifSpinMol2_1 == 1) then
				tmp = tmp*ca21(c,firstDifIndMol2_1)
!				write(*,'(f12.6) ',advance='no') ca21(c,firstDifIndMol2_1)
			else
				tmp = tmp*cb21(c,firstDifIndMol2_1)
!				write(*,'(f12.6) ',advance='no') cb21(c,firstDifIndMol2_1)
			end if
			if(firstDifSpinMol2_2 == 1) then
				tmp = tmp*ca22(d,firstDifIndMol2_2)
!				write(*,'(f12.6) ',advance='no') ca22(d,firstDifIndMol2_2)
			else
				tmp = tmp*cb22(d,firstDifIndMol2_2)
!				write(*,'(f12.6)') cb22(d,firstDifIndMol2_2)
			end if
			Eint = Eint + tmp

!			write(*,*)
!			write(*,*) 'testing int: :',tmp
		end if
		end do
		end do
		end do
		end do

	end function interaction_energy

	subroutine fill_dipole_operator_elements_slater(vec,n,Z,S,DD,MolSlaterType)
		real(dp), dimension(:,:), intent(in)	:: vec,n,S
		real(dp), dimension(:), intent(in)	:: Z
		real(dp), dimension(:,:,:)	, intent(out)	:: DD
		integer(i4b), intent(in)					:: MolSlaterType

		real(dp), dimension(3)			:: DDnuc
		integer(i4b) 					:: dip_ind,a,mu,nu,c,d

		DD = 0.0_dp
		DDnuc = 0.0_dp

		! dipole of the nuclei
		do dip_ind = 1,3
		do a=1, size(Z)
			DDnuc(dip_ind) = DDnuc(dip_ind) + vec(dip_ind,a)/size(Z)
		end do
		end do

		do mu=1,size(DD,2)
		do nu=1,size(DD,3)
		do dip_ind=1,3
			DD(dip_ind, mu,nu) = dipole_xi_slater(vec(:,mu),vec(:,nu),MolSlaterType,MolSlaterType,dip_ind)
		end do
		end do
		end do

		do dip_ind=1,3
			DD(dip_ind,:,:) = -DD(dip_ind,:,:) + S*DDnuc(dip_ind)
		end do

	end  subroutine fill_dipole_operator_elements_slater

	subroutine fill_dipole_operator_elements_slater_directed(vec,n,Z,S,DD)
		real(dp), dimension(:,:), intent(in)	:: vec,n,S
		real(dp), dimension(:), intent(in)	:: Z
		real(dp), dimension(:,:,:)	, intent(out)	:: DD

		real(dp), dimension(3)			:: DDnuc
		integer(i4b) 					:: dip_ind,a,mu,nu,c,d

		if(size(vec,1) /= 3 .or. size(DD,1) /= 3 .or. size(DD,3) /= size(DD,2) .or.&
			size(dd,2) /= size(S,1) .or. size(S,1) /= size(S,2) .or. size(n,1) /= 3 .or.&
			size(n,2) /= size(Z)) then
			call print_error_message(-1, "dimension error in fill_dipole_operator_elements_slater_directed")
		end if

		DD = 0.0_dp
		DDnuc = 0.0_dp

		! dipole of the nuclei
		do dip_ind = 1,3
		do a=1, size(Z)
			DDnuc(dip_ind) = DDnuc(dip_ind) + vec(dip_ind,a)/size(Z)
		end do
		end do

		do mu=1,size(DD,2)
		do nu=1,size(DD,3)
		do dip_ind=1,3
			DD(dip_ind, mu,nu) = dipole_xi_slater_directed(vec(:,mu),vec(:,nu),.true.,.true.,n(:,mu),n(:,nu),dip_ind)
		end do
		end do
		end do

		do dip_ind=1,3
			DD(dip_ind,:,:) = -DD(dip_ind,:,:) + S*DDnuc(dip_ind)
		end do

	end  subroutine fill_dipole_operator_elements_slater_directed

	subroutine fill_intermolecular_four_center_tensor_and_H_slater_directed(vec1,n1,Z1,vec2,n2,Z2,HMol1_2,HMol2_1,four_center_tensor_int_12)
		real(dp), dimension(:,:), intent(in)	:: vec1,n1!,H1
		real(dp), dimension(:), intent(in)	:: Z1
		real(dp), dimension(:,:), intent(in)	:: vec2,n2!,H2
		real(dp), dimension(:), intent(in)	:: Z2

		real(dp), dimension(:,:), intent(out)	:: HMol1_2, HMol2_1
		real(dp), dimension(:,:,:,:), intent(out)	:: four_center_tensor_int_12

		integer(i4b)	:: a,b,c,d,mu,nu,la,si,i,j,k,l

		if(size(vec1,1) /= 3  .or. size(n1,1) /= 3 .or. size(n1,2) /= size(Z1) .or.&
			size(vec2,1) /= 3  .or. size(n2,1) /= 3 .or. size(n2,2) /= size(Z2) ) then
			call print_error_message(-1, "dimension error in fill_intermolecular_four_center_tensor_and_H_slater_directed")
		end if


		four_center_tensor_int_12 = 0.0_dp
		HMol1_2 = 0.0_dp
		HMol2_1 = 0.0_dp

		write(*,*) 'fill_intermolecular_four_center_tensor_and_H_slater_directed'
		do mu=1,size(vec1,2)
		do nu=1,size(vec1,2)
		do c=1,size(vec2,2)
		do d=1,size(vec2,2)
			!write(*,*) 'NEZOBRAZOVAT', mu,nu,c,d

			four_center_tensor_int_12(mu,nu,c,d) = four_center_integral_slater_directed(&
				vec1(:,mu),vec1(:,nu),vec2(:,c),vec2(:,d),.true.,.true.,.true.,.true., &
				n1(:,mu),n1(:,nu),n2(:,c),n2(:,d) )
		end do
		end do
		end do
		end do

		do mu=1,size(vec1,2)
		do nu=1,size(vec1,2)
		do a=1,size(Z2)
			HMol1_2(mu,nu) = HMol1_2(mu,nu) - Z2(a)*coulombic_integral_slater_directed&
						(vec1(:,mu),vec1(:,nu),vec2(:,a),.true.,.true., n1(:,mu),n1(:,nu))
		end do
		end do
		end do

		do c=1,size(vec2,2)
		do d=1,size(vec2,2)
		do a=1,size(Z1)
			HMol2_1(c,d) = HMol2_1(c,d) - Z1(a)*coulombic_integral_slater_directed&
						(vec2(:,c),vec2(:,d),vec1(:,a),.true.,.true., n2(:,c),n2(:,d))
		end do
		end do
		end do


	end  subroutine fill_intermolecular_four_center_tensor_and_H_slater_directed

	subroutine fill_intermolecular_four_center_tensor_and_H_slater(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12,&
			Mol1SlaterType,Mol2SlaterType)
		real(dp), dimension(:,:), intent(in)	:: vec1,n1!,H1
		real(dp), dimension(:), intent(in)	:: Z1
		real(dp), dimension(:,:), intent(in)	:: vec2,n2!,H2
		real(dp), dimension(:), intent(in)	:: Z2
		integer(i4b), intent(in)				:: Mol1SlaterType,Mol2SlaterType

		real(dp), dimension(:,:), intent(out)	:: Hint1_2, Hint2_1
		real(dp), dimension(:,:,:,:), intent(out)	:: four_center_tensor_int_12


		integer(i4b)	:: a,b,c,d,mu,nu,la,si,i,j,k,l

		four_center_tensor_int_12 = 0.0_dp
		Hint1_2 = 0.0_dp
		Hint2_1 = 0.0_dp

		do mu=1,size(vec1,2)
		do nu=1,size(vec1,2)
		do c=1,size(vec2,2)
		do d=1,size(vec2,2)
			four_center_tensor_int_12(mu,nu,c,d) = four_center_integral_slater(&
				vec1(:,mu),vec1(:,nu),vec2(:,c),vec2(:,d),Mol1SlaterType,Mol1SlaterType,Mol2SlaterType,Mol2SlaterType)
		end do
		end do
		end do
		end do

		do mu=1,size(vec1,2)
		do nu=1,size(vec1,2)
		do a=1,size(Z2)
			Hint1_2(mu,nu) = Hint1_2(mu,nu) - Z2(a)*coulombic_integral_slater&
						(vec1(:,mu),vec1(:,nu),vec2(:,a),Mol1SlaterType,Mol1SlaterType)
		end do
		end do
		end do

		do c=1,size(vec2,2)
		do d=1,size(vec2,2)
		do a=1,size(Z1)
			Hint2_1(c,d) = Hint2_1(c,d) - Z1(a)*coulombic_integral_slater&
						(vec2(:,c),vec2(:,d),vec1(:,a),Mol2SlaterType,Mol2SlaterType)
		end do
		end do
		end do

	end  subroutine fill_intermolecular_four_center_tensor_and_H_slater

	subroutine test_hydrogen_molecule(four_center_tensor1, S, H, T, V, E)
		integer(i4b), parameter	:: Dim = 2

		real(dp), dimension(:,:,:,:), intent(out)	:: four_center_tensor1
		real(dp), dimension(:,:), intent(out) 	:: S, H, T, V
		real(dp), intent(out) 						:: E

		real(dp), dimension(3,Dim)		:: vec1, n1
		real(dp), dimension(Dim)		:: Z1, Ea1, Eb1
		real(dp), dimension(Dim,Dim) 	:: ca1, cb1, Fa1, Fb1, Pa11, Pb11, Pa12, Pb12, V2

		logical, dimension(2,size(Z1)) :: excMol1_1,excMol1_2

		integer(i4b) :: i,j,k,l

		if(.not. (size(four_center_tensor1,1) == Dim .and. size(four_center_tensor1,2) == Dim .and. size(four_center_tensor1,3) == Dim .and. size(four_center_tensor1,4) == Dim .and. &
		size(S,1) == Dim .and. size(S,2) == Dim .and. &
		size(H,1) == Dim .and. size(H,2) == Dim .and. &
		size(T,1) == Dim .and. size(T,2) == Dim .and. &
		size(V,1) == Dim .and. size(V,2) == Dim )) then

			call print_error_message(-1, "dimension error in test_hydrogen_molecule")

		end if

		n1(:,1) = (/ 1.0, 0.0, 0.0 /)
		n1(:,2) = (/ 1.0, 0.0, 0.0 /)

		vec1(:,1) = (/ 0.0,  0.0, 0.0 /)
		vec1(:,2) = (/ 0.0,  0.0, 0.0 /)

		Z1 = (/ 1.0, 1.0 /)

		excMol1_1 = .false.
		excMol1_2 = .false.

		excMol1_1(1,1) = .true.
		excMol1_1(2,1) = .true.
		excMol1_2(1,1) = .true.
		excMol1_2(1,2) = .true.

		call init(1.24_dp,1.24_dp)
		vec1(2,1) = 1.40_dp

		call fill_four_center_tensor_and_H_S_slater(four_center_tensor1,H,S, vec1,n1,Z1,T,V,V2,0)
		call solve_Roothan_equations(vec1,n1,Z1,S,H,four_center_tensor1, ca1,cb1,Fa1,Fb1,Ea1,Eb1)
		call fill_P_matrix(Pa11,Pb11,ca1,cb1,excMol1_1)
		call fill_P_matrix(Pa12,Pb12,ca1,cb1,excMol1_2)

		E = full_energy(vec1,n1,Z1,four_center_tensor1,H,Pa11,Pb11)*hartree_in_eV
	end subroutine test_hydrogen_molecule

	subroutine prepare_alkene(vec,n,Z,Dim)
		! prepares linear carotenoid with given dimension along x-axis with
		! normal vectors along z-axis. Carotenoid is centered at the center
		! of coordinate system
		real(dp), dimension(:,:), intent(out)		:: vec, n
		real(dp), dimension(:)	, intent(out)		:: Z
		integer(i4b), intent(in)					:: Dim

		real(dp), parameter				:: alpha = 120.0_dp/180*PI, cc_bond = 1.4*Angstrom_in_Bohr_radii
		real(dp), dimension(3)				:: center

		integer(i4b)	:: i,j,k

		if(size(vec,2) /= size(n,2) .or. size(vec,1) /= 3 .or. size(n,1) /= 3 .or.&
			size(vec,2) /= Dim) then

			call print_error_message(-1, "dimension error in prepare_carotenoid")
		end if

		do i=1, Dim
			Z(i) = 1.0_dp
			vec(3,i) = 0.0_dp

			if(mod(i,2) == 1) then
				vec(2,i) = 0.0_dp
			else
				vec(2,i) = cc_bond*cos(alpha/2)
			end if

			vec(1,i) = cc_bond*sin(alpha/2)*i

			n(1,i) = 0.0_dp
			n(2,i) = 0.0_dp
			n(3,i) = 1.0_dp
		end do

		center = 0.0_dp
		do i=1, Dim
			center = center + vec(:,i)/Dim
		end do

		call rotate_and_move_molecule(vec,n,Z, 0.0_dp,0.0_dp,0.0_dp,-center(1),-center(2),-center(3))
	end subroutine prepare_alkene

	subroutine prepare_carotenoid(vec,n,Z,Dim)
		! prepares linear carotenoid with given dimension along x-axis with
		! normal vectors along z-axis. Carotenoid is centered at the center
		! of coordinate system
		real(dp), dimension(:,:), intent(out)		:: vec, n
		real(dp), dimension(:)	, intent(out)		:: Z
		integer(i4b), intent(in)					:: Dim

!		real(dp), dimension(Dim-2,Dim-2)		:: vvec, nn
!		real(dp), dimension(Dim-2)				:: ZZ
!		integer(i4b)							:: DDim

		real(dp), parameter				:: alpha = 120.0_dp/180*PI, cc_bond = 1.4*Angstrom_in_Bohr_radii
		real(dp), dimension(3)				:: center
		real(dp), dimension(3,3)			:: rot
		integer(i4b)	:: i,j,k

		rot = 0.0
		rot(1,1) = cos(alpha)
		rot(1,2) = sin(alpha)
		rot(2,2) = cos(alpha)
		rot(2,1) = -sin(alpha)
		rot(3,3) = 1.0

		if(size(vec,2) /= size(n,2) .or. size(vec,1) /= 3 .or. size(n,1) /= 3 .or.&
			size(vec,2) /= Dim .or. Dim < 3) then

			call print_error_message(-1, "dimension error in prepare_carotenoid")
		end if

		do i=2, Dim-1
			Z(i) = 1.0_dp
			vec(3,i) = 0.0_dp

			if(mod(i,2) == 0) then
				vec(2,i) = 0.0_dp
			else
				vec(2,i) = cc_bond*cos(alpha/2)
			end if

			vec(1,i) = cc_bond*sin(alpha/2)*i

			n(1,i) = 0.0_dp
			n(2,i) = 0.0_dp
			n(3,i) = 1.0_dp
		end do

		Z(1) = 1.0_dp

		vec(:,1) = vec(:,2) + matmul(rot,vec(:,3)-vec(:,2))

		n(1,1) = 0.0_dp
		n(2,1) = 0.0_dp
		n(3,1) = 1.0_dp

		Z(Dim) = 1.0_dp

		vec(:,Dim) = vec(:,Dim-1) + matmul(rot,vec(:,Dim-2)-vec(:,Dim-1))

		n(1,Dim) = 0.0_dp
		n(2,Dim) = 0.0_dp
		n(3,Dim) = 1.0_dp

		center = 0.0_dp
		do i=1, Dim
			center = center + vec(:,i)/Dim
		end do

		call rotate_and_move_molecule(vec,n,Z, 0.0_dp,0.0_dp,0.0_dp,-center(1),-center(2),-center(3))

!		do i=1,Dim
!			write(*,*) vec(:,i)
!		end do
!		write(*,*)

	end subroutine prepare_carotenoid


	subroutine rotate_and_move_molecule(vec,n,Z, rot1xy, rot2xz, rot3xy, movx, movy, movz)
		real(dp), dimension(:,:), intent(inout)	:: vec, n
		real(dp), dimension(:)	, intent(inout)		:: Z
		real(dp), intent(in)						:: rot1xy, rot2xz, rot3xy, movx, movy, movz

		real(dp), dimension(3,3)					:: A, B, C
		integer(i4b)								:: i,j,k

		if(size(vec,1) /= size(n,1) .or. size(vec,2) /= size(n,2) .or.&
			size(vec,1) /= 3) then

			call print_error_message(-1, "dimension error in rotate_and_move_molecule")
		end if

		A = 0.0_dp
		B = 0.0_dp
		C = 0.0_dp

		A(1,1) = cos(rot1xy)
		A(2,2) = A(1,1)
		A(3,3) = 1.0_dp
		A(1,2) = sin(rot1xy)
		A(2,1) = -A(1,2)

		B(1,1) = cos(rot2xz)
		B(3,3) = B(1,1)
		B(2,2) = 1.0_dp
		B(1,3) = sin(rot2xz)
		B(3,1) = -B(1,3)

		C(1,1) = cos(rot3xy)
		C(2,2) = C(1,1)
		C(3,3) = 1.0_dp
		C(1,2) = sin(rot3xy)
		C(2,1) = -C(1,2)

		do i=1, size(vec,2)
			vec(:,i) 	= matmul(C,matmul(B,matmul(A,vec(:,i))))
			n(:,i) 		= matmul(C,matmul(B,matmul(A,n(:,i))))

			vec(1,i) 	= vec(1,i) + movx
			vec(2,i) 	= vec(2,i) + movy
			vec(3,i) 	= vec(3,i) + movz
		end do

	end subroutine rotate_and_move_molecule

	subroutine read_config_file(vec,n,moleculeCode,status)
		real(dp), dimension(:,:), intent(out) :: vec, n
		integer(i4b), intent(in)	:: moleculeCode
		integer(i4b), intent(out)	:: status

		integer(i4b)	:: i,j,k, Code
		character(len=256)	:: name

		real(dp), dimension(3,0:size(vec,2)+1+1)		:: ReadArrayVec

		open(unit=42,file=trim(file_join(out_dir,"input_carotenoid.dat")))

		status = 0
		vec = 0.0_dp
		n = 0.0_dp
		Code = -42

		i = 0
		do
			if(code == moleculeCode .and. name(1:1) /= '#') then
				i = i + 1
			end if

			! -1 because it starts from 0
			if(i > size(ReadArrayVec,2) - 1) then
				status = -1
				close(42)
				goto 120
			end if

			if(status /= 0) then
				call print_error_message(-1,"should not happen")
			end if

!			write(*,*) i, size(ReadArrayVec,2)

			read(42,*, err=110, end=120) name, Code, ReadArrayVec(1,i), ReadArrayVec(2,i), ReadArrayVec(3,i)

!			write(*,*) 'po cteni'

		end do

	110		call print_error_message(-1, "error reading config file in qch_lib")
	120		close(42)

		do j=1,i-2
			vec(:,j) = ReadArrayVec(:,j)*Angstrom_in_Bohr_radii !HACK
			n(:,j) = cross_product(ReadArrayVec(:,j-1)-ReadArrayVec(:,j),ReadArrayVec(:,j+1)-ReadArrayVec(:,j))

			if(maxval(abs(n(:,j))) > 1e-5) then
				n(:,j) = n(:,j)/sqrt(dot_product(n(:,j),n(:,j)))
			end if
		end do

	end subroutine read_config_file

	subroutine read_dim_from_config_file(Dim,moleculeCode)
		integer(i4b), intent(out)	:: Dim
		integer(i4b), intent(in)	:: moleculeCode

		integer(i4b) :: i, status

		real(dp), dimension(3,100) :: vec, n

		vec = 0.0_dp
		n = 0.0_dp
		Dim = 0

		call read_config_file(vec,n,moleculeCode,status)

		if(status /= 0) then
			call print_error_message(-1, "error in read_dim_from_config_file")
		end if

		do i=1, size(n,2)
			if(dot_product(n(:,i),n(:,i)) < 1e-5) then
				cycle
			end if

			Dim = Dim + 1
		end do

	end subroutine read_dim_from_config_file

	real(dp) function fitness_for_Huckel_parameters_function(vec,n,Z,S, targetHOMO_LUMO,targetD, parA,parB, DD, excMolHOMO,excMolLUMO, includeD) result(fitness)
		real(dp), intent(in)						:: targetHOMO_LUMO
		real(dp), dimension(:), intent(in)		:: targetD
		real(dp), dimension(:,:), intent(in)		:: vec, n
		real(dp), dimension(:), intent(in)		:: Z
		real(dp), dimension(:,:), intent(in)	 	:: S
		real(dp), intent(in)						:: parA, parB
		real(dp), dimension(:,:,:), intent(in)	:: DD
		logical, dimension(:,:), intent(in)	 	:: excMolHOMO,excMolLUMO
		logical, intent(in)						:: includeD

		real(dp), dimension(3)		:: fitD
		real(dp), dimension(size(S,1),size(S,1))		:: H, ca, cb, PaHOMO, PbHOMO, PaLUMO, PbLUMO
		real(dp), dimension(size(S,1))					:: Ea, Eb
		integer(i4b) :: dip_ind

		! we calculate the fitness
		call calculate_Huckel_method(vec,n,Z,S, ca,H,Ea, parA,parB)
		call calculate_Huckel_method(vec,n,Z,S, cb,H,Eb, parA,parB)
		call fill_P_matrix(PaHOMO,PbHOMO,ca,cb,excMolHOMO)
		call fill_P_matrix(PaLUMO,PbLUMO,ca,cb,excMolLUMO)

		do dip_ind=1,3
		fitD(dip_ind) = dipole_moment_between_states(vec,n,Z,vec,n,Z, excMolHOMO,excMolLUMO,&
						ca, cb, ca, cb, DD(dip_ind,:,:))
		end do

		if(fitD(1) < 0.0_dp .or. (fitD(1) == 0.0_dp .and. fitD(2) < 0.0_dp) .or. &
			(fitD(1) == 0.0_dp .and. fitD(2) == 0.0_dp .and. fitD(3) < 0.0_dp)) then
			fitD = -fitD
		end if

		write(*,*) '  actual H-L',full_energy_Huckel(vec,n,Z,H,PaLUMO,PbLUMO) - full_energy_Huckel(vec,n,Z,H,PaHOMO,PbHOMO)
		write(*,*) '  actual D',fitD
		write(*,*) '  target H-L',targetHOMO_LUMO
		write(*,*) '  target D',targetD

		fitness = 	((full_energy_Huckel(vec,n,Z,H,PaLUMO,PbLUMO) - full_energy_Huckel(vec,n,Z,H,PaHOMO,PbHOMO) - targetHOMO_LUMO)/(abs(targetHOMO_LUMO)+1e-5))**2

		if(includeD) then
			fitness = fitness + abs((dot_product(targetD-fitD,targetD-fitD)/(dot_product(targetD,targetD)+1e-5)))
		end if

	end function fitness_for_Huckel_parameters_function

	subroutine fit_Huckel_parameters(vec,n,Z,S, targetHOMO_LUMO,targetDin,excMolHOMO,excMolLUMO, parA,parB)
		real(dp), intent(in)						:: targetHOMO_LUMO
		real(dp), dimension(:), intent(in)		:: targetDin
		real(dp), dimension(:,:), intent(in)		:: vec, n
		real(dp), dimension(:), intent(in)		:: Z
		real(dp), dimension(:,:), intent(in)	 	:: S
		real(dp), intent(out)						:: parA, parB
		logical, dimension(:,:), intent(in) :: excMolHOMO,excMolLUMO

		real(dp), parameter							:: delta = 0.01_dp, gamma = 0.1_dp

		real(dp), dimension(size(S,1),size(S,1))		:: H, ca, cb, Pa1, Pb1, Pa2, Pb2
		real(dp), dimension(size(S,1))					:: Ea, Eb
		real(dp), dimension(3,size(S,1),size(S,1))		:: DD
		real(dp), dimension(3)							:: targetD

		integer(i4b) :: i,j,k,l, c,d,mu,nu, a, dip_ind, Dim

		logical, parameter :: includeD = .false.

		real(dp)				:: fitness, grad_fitness_A, grad_fitness_B
		real(dp)				:: best_fitness, best_parB

		if(.not.(size(vec,1) == 3 .and. size(n,1) == 3) .and. size(targetD) == 3) then
			call print_error_message(-1, "error in fit_Huckel_parameters")
		end if

		if(targetDin(1) < 0.0_dp .or. (targetDin(1) == 0.0_dp .and. targetDin(2) < 0.0_dp) .or. &
			(targetDin(1) == 0.0_dp .and. targetDin(2) == 0.0_dp .and. targetDin(3) < 0.0_dp)) then
			targetD = -targetDin
		else
			targetD = targetDin
		end if

		write(*,*)
		write(*,*) 'target H-L', targetHOMO_LUMO
		write(*,*) 'target D', targetD
		write(*,*), 'HOMO',excMolHOMO
		write(*,*), 'LUMO',excMolLUMO

		Dim = size(S,1)
		parA = -11.4/Hartree_in_eV
		parB = -11.4/Hartree_in_eV*1.75

		call fill_dipole_operator_elements_slater_directed(vec,n,Z,S,DD)


		do i=1, 1000
			fitness = fitness_for_Huckel_parameters_function(vec,n,Z,S, targetHOMO_LUMO,targetD, parA,parB, DD, excMolHOMO,excMolLUMO,includeD)
			grad_fitness_A = (fitness_for_Huckel_parameters_function(vec,n,Z,S, targetHOMO_LUMO,targetD, parA+delta,parB, DD, excMolHOMO,excMolLUMO,includeD)&
								- fitness)/delta
			grad_fitness_B = (fitness_for_Huckel_parameters_function(vec,n,Z,S, targetHOMO_LUMO,targetD, parA,parB+delta, DD, excMolHOMO,excMolLUMO,includeD)&
								- fitness)/delta

			!parA = parA - gamma*grad_fitness_A
			parA = 0.0_dp
			parB = parB - gamma*grad_fitness_B

			write(*,*) 'iteration',i,':', fitness,grad_fitness_A,grad_fitness_B, parA, parB
		end do

!		best_fitness = 1e10_dp
!		best_parB = 1e10_dp
!		do i=1, 1000
!			parA = 0.0_dp
!			parB = (i-500)*0.002_dp
!
!			fitness = fitness_for_Huckel_parameters_function(vec,n,Z,S, targetHOMO_LUMO,targetD, parA,parB, DD, excMolHOMO,excMolLUMO, &
!				includeD)
!
!			if(fitness < best_fitness) then
!				best_fitness = fitness
!				best_parB = parB
!			end if
!
!			write(*,*) 'iteration',i,':', best_fitness, best_parB
!		end do
!
!		parA = 0.0_dp
!		parB = best_parB

	end subroutine fit_Huckel_parameters

	real(dp) function model_surface_of_two_blocks(phi,x,y,Xlength,Ylength,Zlength) result(S)
		real(dp), intent(in)		:: phi,x,y
		real(dp), intent(in)		:: Xlength,Ylength,Zlength

		real(dp), dimension(3,3)	:: RM ! rot move matrix
		real(dp), dimension(3)		:: vec ! rot move matrix
		real(dp)					:: commonXY

		integer(i4b)	:: i,j

		S = 0.0_dp

		! the rest is obtaining areal common to both blocks
		RM = 0.0_dp
		RM(1,1) = cos(phi)
		RM(2,2) = cos(phi)
		RM(1,2) = sin(phi)
		RM(2,1) = -sin(phi)
		RM(1,3) = x
		RM(2,3) = y
		RM(3,3) = 1.0_dp

		vec = (/ 1, 0, 1 /)

		vec = matmul(RM,vec)

!		write(*,*) "TESTING model_surface_of_two_blocks", vec

		commonXY = 0.0_dp
		do i=-100, 100
		do j=-100, 100
			vec = 0.0_dp
			vec(1) = i/201.0_dp*Xlength
			vec(2) = j/201.0_dp*Ylength
			vec(3) = 1
			vec = matmul(RM, vec)

			if(vec(1) > -Xlength/2 .and. vec(1) <= Xlength/2 .and. &
			   vec(2) > -Ylength/2 .and. vec(2) <= Ylength/2) then

				commonXY = commonXY + 1.0_dp/201/201*Xlength*Ylength
			end if

		end do
		end do

!		write(*,*) "TESTING model_surface_of_two_blocks, commonS",  commonXY
!		stop

		! surface
		S = 2*(Xlength*Ylength+2*Xlength*Zlength+2*Ylength*Zlength) &
			+ 2*Xlength*Ylength - commonXY

	end function model_surface_of_two_blocks

	subroutine pokus()
		integer(i4b), parameter	:: Dim = 2

		real(dp), dimension(3,Dim)		:: vec1, n1
		real(dp), dimension(Dim)		:: Z1, Ea1, Eb1
		real(dp), dimension(Dim,Dim) 	:: S1, H1, ca1, cb1, Fa1, Fb1, Pa11, Pb11, Pa12, Pb12, V1, V2, T

		real(dp), dimension(3,Dim)		:: vec2, n2
		real(dp), dimension(Dim)		:: Z2, Ea2, Eb2
		real(dp), dimension(Dim,Dim) 	:: S2, H2, ca2, cb2, Fa2, Fb2, Pa21, Pb21, Pa22, Pb22

		real(dp), dimension(Dim,Dim,Dim,Dim)	:: four_center_tensor1, four_center_tensor2

		logical, dimension(2,size(Z1)) :: excMol1_1,excMol1_2
		logical, dimension(2,size(Z2)) :: excMol2_1,excMol2_2

		integer(i4b) :: i,j,k,l
		logical :: tmp_bool
		real(dp) :: tmp_dp

		n1(:,1) = (/ 1.0, 0.0, 0.0 /)
		n1(:,2) = (/ 1.0, 0.0, 0.0 /)
		n2(:,1) = (/ 1.0, 0.0, 0.0 /)
		n2(:,2) = (/ 1.0, 0.0, 0.0 /)

		vec1(:,1) = (/ 0.0,  0.0, 0.0 /)
		vec1(:,2) = (/ 0.0,  0.0, 0.0 /)
		vec2(:,1) = (/ 1000.0, -1.0, 0.0 /)
		vec2(:,2) = (/ 1000.0,  1.0, 0.0 /)

		Z1 = (/ 1.0, 1.0 /)
		Z2 = (/ 1.0, 1.0 /)

!		Z1 = (/ 1.0 /)

		excMol1_1 = .false.
		excMol1_2 = .false.
		excMol2_1 = .false.
		excMol2_2 = .false.

		excMol1_1(1,1) = .true.
		excMol1_1(2,1) = .true.
		excMol1_2(1,1) = .true.
		excMol1_2(1,2) = .true.
		excMol2_1(1,1) = .true.
		excMol2_1(2,1) = .true.
		excMol2_2(1,2) = .true.
		excMol2_2(2,1) = .true.

		write(*,*) 'POKUS called'

	call init_random_seed()
	call random_number(tmp_dp)

!!   plotting subunit
!	open(unit=112,file='/home/olsij4am/prace/nose-2Dplot-overlap.dat')
!
!		do i=0,99
!		do j=1,100
!			write(112,'(f12.6 )', advance='no') model_surface_of_two_blocks((i-50)*PI/100.0_dp,(j-1)*10/100.0_dp,0.0_dp,10.0_dp,1.0_dp,1.0_dp)
!		end do
!			write(112,*)
!		end do
!	close(112)
!	stop


!	if(tmp_dp > 0.5) then
!		call pokus4(24)
!	else
		call pokus3()
!	end if

!	call pokus2()
	stop

	IF(.TRUE.) THEN
		call pokus2()
		call pokus3()
	ELSE
		open(unit=42,file='/home/olsij4am/prace/nose-molecular-energies-C-C.dat')
		open(unit=43,file='/home/olsij4am/prace/nose-molecular-energies-C-C--H.dat')
		open(unit=44,file='/home/olsij4am/prace/nose-molecular-energies-C-C--c.dat')
		open(unit=45,file='/home/olsij4am/prace/nose-molecular-energies-C-C--v.dat')
		open(unit=52,file='/home/olsij4am/prace/nose-molecular-energies-C-C2.dat')
		open(unit=53,file='/home/olsij4am/prace/nose-molecular-energies-C-C--H2.dat')
		open(unit=54,file='/home/olsij4am/prace/nose-molecular-energies-C-C--c2.dat')
		open(unit=55,file='/home/olsij4am/prace/nose-molecular-energies-C-C--v2.dat')

		call init(1.24_dp,1.24_dp)

		do i=1,300
			vec1(2,1) = i/40.0_dp
!			vec1(2,1) = 1.40_dp
!			exponent_slater = 1.24 !1/(i/150.0_dp)

!			STOxi = (/ 0.185742_dp, 1.12_dp,  13.5413_dp /) * exponent_slater**2
!			STOciS = (/ 0.276651_dp,  0.433794_dp,  0.242884_dp /)
!			call init()

!			call pokus2()
!			stop
!			cycle

			call fill_four_center_tensor_and_H_S_slater(four_center_tensor1,H1,S1, vec1,n1,Z1,T,V1,V2,0)
			call solve_Roothan_equations(vec1,n1,Z1,S1,H1,four_center_tensor1, ca1,cb1,Fa1,Fb1,Ea1,Eb1)
!			call solve_Roothan_equations(vec1,n1,Z1,S1,H1,four_center_tensor1, ca11,cb11,Fa11,Fb11,Ea11,Eb11,excMol1_1)
!			call solve_Roothan_equations(vec1,n1,Z1,S1,H1,four_center_tensor1, ca12,cb12,Fa12,Fb12,Ea12,Eb12,excMol1_2)
			call fill_P_matrix(Pa11,Pb11,ca1,cb1,excMol1_1)
			call fill_P_matrix(Pa12,Pb12,ca1,cb1,excMol1_2)

			write(*,*)
			write(*,*) '   T > ', T
			write(*,*) '   V1> ', V1
			write(*,*) '   V2> ', V2
			write(*,*)

			write(*,*) 'molecule 1 solved:'
			write(*,*) '---------->',vec1(2,1),full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)

!			call calculate_molecule(vec2,n2,Z2, S2,H2,four_center_tensor2)
!			call solve_Roothan_equations(vec2,n2,Z2,S2,H2,four_center_tensor2, ca21,cb21,Fa21,Fb21,Ea21,Eb21,excMol2_1)
!			call solve_Roothan_equations(vec2,n2,Z2,S2,H2,four_center_tensor2, ca22,cb22,Fa22,Fb22,Ea22,Eb22,excMol2_2)
!			write(*,*) 'molecule 2 solved'

!			write(42,*), vec1(2,1),sqrt(dot_product(matmul(S1,ca1(:,1)),matmul(S1,ca1(:,1)))),sqrt(dot_product(matmul(S1,ca1(:,2)),matmul(S1,ca1(:,2))))!,full_energy(four_center_tensor1,H1,Pa11,Pb11), S1(1,2), Ea1
			write(42,*), vec1(2,1)*bohr_radius_in_pm,min(max(full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)*hartree_in_eV,-20000.0_dp),20000.0_dp)
			write(43,*), vec1(2,1)*bohr_radius_in_pm!,tmp_full_energy1(four_center_tensor1,H1,Pa11,Pb11)*hartree_in_eV
			write(44,*), vec1(2,1)*bohr_radius_in_pm,overlap_slater(vec1(:,1),vec1(:,1),0,0)
!			write(45,*), vec1(2,1)*bohr_radius_in_pm,density_slater(vec1(:,1),vec1(:,2),0),&
!						exp(-exponent_slater*sqrt(dot_product(vec1(:,1)-vec1(:,2),vec1(:,1)-vec1(:,2))))/sqrt(PI/exponent_slater**3)
			call flush(42)
			call flush(43)
			call flush(44)
			call flush(45)
			write(52,*), vec1(2,1)*bohr_radius_in_pm,min(max(full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa12,Pb12)*hartree_in_eV,-20000.0_dp),20000.0_dp)
			write(53,*), vec1(2,1)*bohr_radius_in_pm!,tmp_full_energy1(four_center_tensor1,H1,Pa12,Pb12)*hartree_in_eV
			write(54,*), vec1(2,1)*bohr_radius_in_pm!,(tmp_full_energy2(four_center_tensor1,H1,Pa12,Pb12) + tmp_full_energy3(four_center_tensor1,H1,Pa12,Pb12))*hartree_in_eV
			write(55,*), vec1(2,1)*bohr_radius_in_pm!,(tmp_full_energy3(four_center_tensor1,H1,Pa12,Pb12) - tmp_full_energy3(four_center_tensor1,H1,Pa12,Pb12))*hartree_in_eV
			call flush(52)
			call flush(53)
			call flush(54)
			call flush(55)

		end do
		close(42)
		close(43)
		close(44)
		close(45)
		close(52)
		close(53)
		close(54)
		close(55)
	END IF

	end subroutine pokus

	subroutine pokus2()

		real(dp), dimension(3) 	:: vecA,vecB,vecC,vecD, nA,nB,nC,nD, vecE, vecF, vecG
		real(dp), dimension(2,2)	:: tmp, S_, eigvec, values

		integer(i4b)				:: i,j,k

		call init(exponent_slaterS,exponent_slaterP)

		vecA = 0.0_dp
		vecB = (/ 0.0_dp, 0.0_dp, 4.1_dp*Angstrom_in_Bohr_radii /)
		vecC = (/ 1.4_dp*Angstrom_in_Bohr_radii, 0.0_dp, 0.0_dp /)
		vecD = (/ 0.0_dp, 0.0_dp,    0.0_dp /)
		vecE = (/ 0.0_dp, 1.0_dp,    0.0_dp /)
		vecF = (/ 0.0_dp, 0.0_dp, 1000.0_dp /)
		vecG = (/ 0.0_dp, 1.0_dp, 1000.0_dp /)

		nA = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
		nB = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
		nC = (/ 1.0_dp, 2.0_dp, 3.0_dp /)
		nD = (/ 1.0_dp, 1.0_dp, 1.0_dp /)

		nA = nA/sqrt(dot_product(nA,nA))
		nB = nB/sqrt(dot_product(nB,nB))
		nC = nC/sqrt(dot_product(nC,nC))
		nD = nD/sqrt(dot_product(nD,nD))

		write(*,*)
		write(*,*) 'overlap 00', overlap(vecB,vecB,2.6_dp,2.6_dp,0,0)
		write(*,*) 'overlap 11', overlap(vecB,vecB,2.6_dp,2.6_dp,1,1)

		write(*,*) 'four_center',four_center_integral(vecB,vecB,nA,nA, &
						2.6_dp,2.6_dp,2.6_dp,2.6_dp, 1,1,1,1)/&
						overlap(vecB,vecB, 2.6_dp,2.6_dp, 1,1)/&
						overlap(vecD,vecD, 2.6_dp,2.6_dp, 1,1)
		write(*,*) 'four_center0',four_center_integral(vecB,vecB,vecB,vecB, &
						2.6_dp,2.6_dp,2.6_dp,2.6_dp, 1,1,1,1)/&
						overlap(vecB,vecB, 2.6_dp,2.6_dp, 1,1)/&
						overlap(vecD,vecD, 2.6_dp,2.6_dp, 1,1)
		write(*,*) 'four_centerS',four_center_integral(vecB,vecB,nA,nA, &
						2.6_dp,2.6_dp,2.6_dp,2.6_dp, 0,0,0,0)/&
						overlap(vecB,vecB, 2.6_dp,2.6_dp, 0,0)/&
						overlap(vecD,vecD, 2.6_dp,2.6_dp, 0,0)
		write(*,*) 'four_centerS0',four_center_integral(vecB,vecB,vecB,vecB, &
						2.6_dp,2.6_dp,2.6_dp,2.6_dp, 0,0,0,0)/&
						overlap(vecB,vecB, 2.6_dp,2.6_dp, 0,0)/&
						overlap(vecD,vecD, 2.6_dp,2.6_dp, 0,0)

		write(*,*) 'coulombic 00',coulombic_integral(vecB,vecB,vecD, &
						1.0_dp,1.0_dp, 0,0)/&
						overlap(vecB,vecB, 1.0_dp,1.0_dp, 0,0)
		write(*,*) 'coulombic 10',coulombic_integral(vecB,vecB,vecD, &
						1.0_dp,1.0_dp, 1,0)/&
						overlap(vecB,vecB, 1.0_dp,1.0_dp, 1,0)
		write(*,*) 'coulombic 11',coulombic_integral(vecB,vecB,vecD, &
						2.0_dp,2.0_dp, 1,1)/&
						overlap(vecB,vecB, 2.0_dp,2.0_dp, 1,1)
		write(*,*) 'kinetic 00',kinetic_integral(vecB,vecB, &
						2.0_dp,2.0_dp, 0,0)/&
						overlap(vecB,vecB, 2.0_dp,2.0_dp, 0,0)
		write(*,*) 'kinetic 10',kinetic_integral(vecB,vecB, &
						2.0_dp,2.0_dp, 1,0)/&
						overlap(vecB,vecB, 1.0_dp,1.0_dp, 1,0)
		write(*,*) 'kinetic 11',kinetic_integral(vecB,vecB, &
						2.0_dp,2.0_dp, 1,1)/&
						overlap(vecB,vecB, 2.0_dp,2.0_dp, 1,1)

		write(*,*)
		write(*,*) 'overlap_slater 00', overlap_slater(vecB,vecB,0,0)
		write(*,*) 'overlap_slater 11', overlap_slater(vecB,vecB,1,1)

		write(*,*) 'dipole_slater 01-1', dipole_xi_slater(vecB,vecC,0,1,1)
		write(*,*) 'dipole_slater 01-1', dipole_xi_slater(vecB,vecC,0,1,1)
		write(*,*) 'dipole_slater 10-1', dipole_xi_slater(vecB,vecC,1,0,1)
		write(*,*) 'dipole_slater 10-1', dipole_xi_slater(vecB,vecC,1,0,1)
		write(*,*) 'dipole_slater 00-1', dipole_xi_slater(vecB,vecC,0,0,1)
		write(*,*) 'dipole_slater 01-2', dipole_xi_slater(vecB,vecC,0,1,2)
		write(*,*) 'dipole_slater 01-2', dipole_xi_slater(vecB,vecC,0,1,2)
		write(*,*) 'dipole_slater 02-1', dipole_xi_slater(vecB,vecC,0,2,1)
		write(*,*) 'dipole_slater 12-1', dipole_xi_slater(vecB,vecC,1,2,1)
		write(*,*) 'dipole_slater 21-1', dipole_xi_slater(vecB,vecC,2,1,1)
		write(*,*) 'dipole_slater 11-1', dipole_xi_slater(vecB,vecC,1,1,1)
		write(*,*) 'dipole_slater 11-2', dipole_xi_slater(vecB,vecC,1,1,2)
		write(*,*) 'dipole_slater 12-3', dipole_xi_slater(vecB,vecC,1,2,3)
		write(*,*) 'dipole_slater 33-3', dipole_xi_slater(vecD,vecE,3,3,3)
		write(*,*) 'dipole_slater 33-3', dipole_xi_slater(vecF,vecG,3,3,3)
		write(*,*) 'dipole_slater 33-3', dipole_xi_slater(vecG,vecG,3,3,3)


		write(*,*) 'overlap_slater', overlap_slater(vecD,vecD, 0,0)
		write(*,*) 'four_center_slater',four_center_integral_slater(vecB,vecB,vecD,vecD, &
						0,0,0,0)
		write(*,*) 'coulombic_slater',coulombic_integral_slater(vecB,vecB,vecD, &
						0,0)
		write(*,*) 'kinetic_slater', kinetic_integral_slater(vecD,vecD, 0,0)
		write(*,*)

		call init(1.65_dp,1.65_dp)
		write(*,*) 'overlapD CC-aggregate', overlap_slater_directed(vecA,vecB,.true.,.true., nB,nB)
		write(*,*) 'overlapD CC-bond', overlap_slater_directed(vecA,vecC,.true.,.true., nB,nB)
		write(*,*) 'overlapD CC-self', overlap_slater_directed(vecA,vecA,.true.,.true., nB,nB)
		write(*,*) 'overlap CC-aggregate', overlap_slater(vecA,vecB,3,3)
		write(*,*) 'overlap CC-bond', overlap_slater(vecA,vecC,3,3)
		write(*,*) 'overlap CC-self', overlap_slater(vecA,vecA,3,3)


		write(*,*) 'tests of Skala book'
		write(*,*) '----------------------------------------'

!		do i=1,10

!		nuclei_shift_factor = i/20.0_dp

!		write(*,*)
!		write(*,*) 'nuclei_shift_factor', nuclei_shift_factor

		call init(2.6_dp,2.6_dp)
		write(*,*) '(3|3)', overlap_slater(vecB,vecB,1,1)
		write(*,*) '(3+4+5|3+4+5)', overlap_slater_directed(vecB,vecB,.true.,.true., nD,nD)
		write(*,*) '(33|33)', four_center_integral_slater(vecB,vecB,vecB,vecB, &
						1,1,1,1)
		write(*,*) '(33+44+55|33+44+55)', four_center_integral_slater_directed(&
					vecB,vecB,vecB,vecB, .true.,.true.,.true.,.true., nD,nD,nD,nD)

		call init(8.7_dp,8.7_dp)
		write(*,*) '(11|11)', four_center_integral_slater(vecB,vecB,vecB,vecB, &
						0,0,0,0)

		call init(1.0_dp,1.0_dp)
		write(*,*) '(66|66)', four_center_integral_slater(vecB,vecB,vecB,vecB, &
						0,0,0,0)

!		enddo

!		stop


	end subroutine pokus2

	subroutine pokus3()
		integer(i4b)	:: Dim1, Dim2, Dim

		real(dp), dimension(:,:), allocatable		:: vec1, n1
		real(dp), dimension(:), allocatable		:: Z1, Ea1, Eb1
		real(dp), dimension(:,:), allocatable 	:: S1, H1, ca1, cb1, Fa1, Fb1, Pa11, Pb11, Pa12, Pb12, Hint1_2

		real(dp), dimension(:,:), allocatable		:: vec2, n2
		real(dp), dimension(:), allocatable		:: Z2, Ea2, Eb2
		real(dp), dimension(:,:), allocatable 	:: S2, H2, ca2, cb2, Fa2, Fb2, Pa21, Pb21, Pa22, Pb22, Hint2_1, HHuckel

		real(dp), dimension(:,:,:,:), allocatable	:: four_center_tensor1, four_center_tensor2, four_center_tensor_int_12

		real(dp), dimension(:,:,:), allocatable		:: DD1, DD2
		real(dp)						:: tmp1, tmp2, tmpRabs, tmpDipE, tmpNeDipE
		real(dp), dimension(3)			:: tmpRR, tmpDD1, tmpDD2, DDnuc1, DDnuc2

		logical, dimension(:,:), allocatable :: excMol1_1,excMol1_2, excMolHOMO, excMolLUMO
		logical, dimension(:,:), allocatable :: excMol2_1,excMol2_2

		real(dp)	:: parA, parB
		real(dp)	:: targetHOMO_LUMO
		real(dp), dimension(3)	:: targetD

		integer(i4b) 	:: i,j,k,l, c,d,mu,nu, a, dip_ind, phi, x
		real(dp)		:: phi_D, x_D

		!construction of linear carotenoid and calculation through extended Huckel

		write(*,*) 'POKUS3 called'

		call read_dim_from_config_file(Dim1,1)
		write(*,*) 'DIM=',Dim1
		call read_dim_from_config_file(Dim2,2)
		write(*,*) 'DIM=',Dim2

		ALLOCATE(vec1, (3,Dim1))
		ALLOCATE(n1, (3,Dim1))
		ALLOCATE(Z1, (Dim1))
		ALLOCATE(Ea1, (Dim1))
		ALLOCATE(Eb1, (Dim1))
		ALLOCATE(S1, (Dim1,Dim1))
		ALLOCATE(H1, (Dim1,Dim1))
		ALLOCATE(ca1, (Dim1,Dim1))
		ALLOCATE(cb1, (Dim1,Dim1))
		ALLOCATE(Fa1, (Dim1,Dim1))
		ALLOCATE(Fb1, (Dim1,Dim1))
		ALLOCATE(Pa11, (Dim1,Dim1))
		ALLOCATE(Pb11, (Dim1,Dim1))
		ALLOCATE(Pa12, (Dim1,Dim1))
		ALLOCATE(Pb12, (Dim1,Dim1))
		ALLOCATE(Hint1_2, (Dim1,Dim1))
		ALLOCATE(four_center_tensor1, (Dim1,Dim1,Dim1,Dim1))
		ALLOCATE(excMol1_1, (2,Dim1))
		ALLOCATE(excMol1_2, (2,Dim1))
		ALLOCATE(DD1, (3,Dim1,Dim1))

		ALLOCATE(four_center_tensor_int_12,  (Dim1,Dim1,Dim2,Dim2))

		ALLOCATE(vec2, (3,Dim2))
		ALLOCATE(n2, (3,Dim2))
		ALLOCATE(Z2, (Dim2))
		ALLOCATE(Ea2, (Dim2))
		ALLOCATE(Eb2, (Dim2))
		ALLOCATE(S2, (Dim2,Dim2))
		ALLOCATE(H2, (Dim2,Dim2))
		ALLOCATE(ca2, (Dim2,Dim2))
		ALLOCATE(cb2, (Dim2,Dim2))
		ALLOCATE(Fa2, (Dim2,Dim2))
		ALLOCATE(Fb2, (Dim2,Dim2))
		ALLOCATE(Pa21, (Dim2,Dim2))
		ALLOCATE(Pb21, (Dim2,Dim2))
		ALLOCATE(Pa22, (Dim2,Dim2))
		ALLOCATE(Pb22, (Dim2,Dim2))
		ALLOCATE(Hint2_1, (Dim2,Dim2))
		ALLOCATE(four_center_tensor2, (Dim2,Dim2,Dim2,Dim2))
		ALLOCATE(excMol2_1, (2,Dim2))
		ALLOCATE(excMol2_2, (2,Dim2))
		ALLOCATE(DD2, (3,Dim2,Dim2))

		ALLOCATE(excMolHOMO, (2,Dim2))
		ALLOCATE(excMolLUMO, (2,Dim2))

		ALLOCATE(HHuckel, (Dim2,Dim2))

		call read_config_file(vec1,n1,1,i)
		write(*,*) 'STATUS=',i
		write(*,*) 'VEC=',vec1
		write(*,*) 'N=',n1

		call read_config_file(vec2,n2,2,i)
		write(*,*) 'STATUS=',i
		write(*,*) 'VEC=',vec2
		write(*,*) 'N=',n2

		Z1 = 1.0_dp
		Z2 = 1.0_dp



		if(Dim1 == Dim2) then
			Dim = Dim1
		else
			call print_error_message(-1,"function not prepared for unequal dimensions")
		end if

		call init(1.65_dp,1.65_dp)

		call fill_four_center_tensor_and_H_S_slater_directed(four_center_tensor1,H1,S1, vec1,n1,Z1,.true.)
		call fill_four_center_tensor_and_H_S_slater_directed(four_center_tensor2,H2,S2, vec2,n2,Z2,.true.)


		write(*,*) 'VEC1=',vec1
		write(*,*) 'VEC2=',vec2
		write(*,*) 'N1=',n1
		write(*,*) 'N2=',n2
		write(*,*) 'S1=',S1
		write(*,*) 'S2=',S2
		write(*,*) 'Dim', Dim

		write(*,*) 'Ea1',Ea1*hartree_in_invCm
		stop

		call fill_HOMO_state(excMolHOMO, Dim)
		call fill_LUMO_state(excMolLUMO, Dim)

		call fill_intermolecular_four_center_tensor_and_H_slater_directed(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12)

		call fill_dipole_operator_elements_slater_directed(vec1,n1,Z1,S1,DD1)
		call fill_dipole_operator_elements_slater_directed(vec2,n2,Z2,S2,DD2)


		targetD = 0.0_dp !(/ -1.07, 0.0, 2.85 /)
		targetHOMO_LUMO = 4.7484_dp/Hartree_in_eV	! 6
		parA = -11.4_dp/Hartree_in_eV
		parB = parA*1.75_dp*3.5

		call calculate_Huckel_method(vec1,n1,Z1,S1, ca1,HHuckel,Ea1, parA,parB)
		call calculate_Huckel_method(vec1,n1,Z1,S1, cb1,HHuckel,Eb1, parA,parB)
		call calculate_Huckel_method(vec2,n2,Z2,S2, ca2,HHuckel,Ea2, parA,parB)
		call calculate_Huckel_method(vec2,n2,Z2,S2, cb2,HHuckel,Eb2, parA,parB)


		DO k=Dim/2-1,Dim/2
		DO a=1,Dim/2+1+1

		call fill_HOMO_state(excMol1_1, Dim)
		call fill_HOMO_state(excMol2_1, Dim)

		call fill_HOMO_state(excMol1_2, Dim)
		call fill_HOMO_state(excMol2_2, Dim)

		if(.not.((excMol1_2(1,k) .eqv. .true.) .and. (excMol1_2(1,a) .eqv. .false.))) then
			cycle
		end if

		excMol1_2(1,k) = .false.
		excMol1_2(1,a) = .true.
		excMol2_2(1,k) = .false.
		excMol2_2(1,a) = .true.

		write(*,*) excMol1_1,k,a
		write(*,*) excMol1_2,k,a
		write(*,*) excMol2_1,k,a
		write(*,*) excMol2_2,k,a

		call fill_P_matrix(Pa11,Pb11,ca1,cb1,excMol1_1)
		call fill_P_matrix(Pa12,Pb12,ca1,cb1,excMol1_2)
		call fill_P_matrix(Pa21,Pb21,ca2,cb2,excMol2_1)
		call fill_P_matrix(Pa22,Pb22,ca2,cb2,excMol2_2)

		write(*,*) 'trace P----', trace(matmul(Pa11+Pb11,S1)),trace(matmul(Pa12+Pb12,S1)),trace(matmul(Pa21+Pb21,S2)),trace(matmul(Pa22+Pb22,S2))
		write(*,*) 'HR-----------------', H1
		write(*,*) 'Hü-----------------', HHuckel
		write(*,*) '----'
		write(*,*) 'ca1----', ca1(:,k)
		write(*,*) 'ca2----', ca1(:,a)
		write(*,*)  '----'
		write(*,*) 'ca----', ca1(:,1)
		write(*,*) 'ca----', ca1(:,2)
		write(*,*) 'ca----', ca1(:,3)
		write(*,*) 'ca----', ca1(:,4)
		write(*,*) 'ca----', ca1(:,5)
		write(*,*) 'ca----', ca1(:,6)
		write(*,*)  '----'

		write(*,*) 'E1-L - E1-H', full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa12,Pb12)*hartree_in_invCm - full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa11,Pb11)*hartree_in_invCm,'cm^-1'
		write(*,*) 'E2-L - E2-H', full_energy_Huckel(vec2,n2,Z2,HHuckel,Pa22,Pb22)*hartree_in_invCm - full_energy_Huckel(vec2,n2,Z2,HHuckel,Pa21,Pb21)*hartree_in_invCm,'cm^-1'
		write(*,*) 'RE1-L - RE1-H', full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa12,Pb12)*hartree_in_invCm - full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)*hartree_in_invCm,'cm^-1'
		write(*,*) 'RE2-L - RE2-H', full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa22,Pb22)*hartree_in_invCm - full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa21,Pb21)*hartree_in_invCm,'cm^-1'

		call sleep(1)

open(unit=112,file='/home/olsij4am/prace/nose-carotenoidJ-fromfile.dat')
DO phi=0,0
DO x = 0,0

phi_D = phi/11.0_dp*PI
x_D = x*2.0_dp*Angstrom_in_Bohr_radii

!call prepare_carotenoid(vec1,n1,Z1,Dim)
!call prepare_carotenoid(vec2,n2,Z2,Dim)

!call rotate_and_move_molecule(vec2,n2,Z2, phi_D,0.0_dp,0.0_dp,x_D,0.0_dp,4.1_dp*Angstrom_in_Bohr_radii)

call fill_intermolecular_four_center_tensor_and_H_slater_directed(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12)

call fill_dipole_operator_elements_slater_directed(vec1,n1,Z1,S1,DD1)
call fill_dipole_operator_elements_slater_directed(vec2,n2,Z2,S2,DD2)

write(*,*) 'VEC1=',vec1
write(*,*) 'VEC2=',vec2

		tmpNeDipE = interaction_energy(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,excMol2_1,excMol2_2, &
							four_center_tensor1, H1, four_center_tensor2, H2, &
							ca1, cb1, ca1, cb1, ca2, cb2, ca2, cb2,&
							four_center_tensor_int_12, Hint1_2, Hint2_1)
		write(*,*) 'interaction energy:', tmpNeDipE,'=', tmpNeDipE*hartree_in_invCm,'cm^-1'

		tmpRR = 0.0_dp
		do i=1,Dim
			tmpRR = tmpRR + (vec1(:,i)-vec2(:,i))/Dim
		end do

		write(*,*)
		write(*,*) 'R=', sqrt(dot_product(tmpRR,tmpRR))
		write(*,*) 'v1=', vec1(:,1)
		write(*,*) "v1'=", vec1(:,2)
		write(*,*) 'v2=', vec2(:,1)
		write(*,*) "v2'=", vec2(:,2)

		do dip_ind=1,3
			tmpDD1(dip_ind) = dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMol1_1,excMol1_2,&
							ca1, cb1, ca1, cb1, DD1(dip_ind,:,:))
			tmpDD2(dip_ind) = dipole_moment_between_states(vec2,n2,Z2,vec2,n2,Z2, excMol2_1,excMol2_2,&
							ca2, cb2, ca2, cb2, DD2(dip_ind,:,:))
		end do


		write(*,*) 'd1',tmpDD1,sqrt(dot_product(tmpDD1,tmpDD1))
		write(*,*) 'd2',tmpDD2,sqrt(dot_product(tmpDD2,tmpDD2))

		tmpRabs = sqrt(dot_product(tmpRR,tmpRR))
		tmpDipE = dot_product(tmpDD1,tmpDD2)/tmpRabs**3 - 3*dot_product(tmpRR,tmpDD1)*dot_product(tmpRR,tmpDD2)/tmpRabs**5
		write(*,*)
		write(*,*) 'DIPOLE INTERACTION ENERGY', tmpDipE,'=',tmpDipE*hartree_in_invCm,'cm^-1'

		if(abs(tmpDipE-tmpNeDipE)/abs(tmpDipE) > 0.2) then
			write(*,*) 'ENERGIES DIFFER MORE THAN 20%'
		end if

		write(*,*)
		write(*,*)
		write(112,*) x_D/Angstrom_in_Bohr_radii, phi_D, tmpDipE*hartree_in_invCm, tmpNeDipE*hartree_in_invCm, sign(1.0_dp,dot_product(tmpDD1,tmpDD2))
		call flush(112)
END DO
END DO
close(112)

		END DO
		END DO

	end subroutine pokus3

	subroutine pokus4(Dim)
		integer(i4b), intent(in)	:: Dim !, outer_integer_parameter_A

		real(dp), dimension(3,Dim)		:: vec1, n1
		real(dp), dimension(Dim)		:: Z1, Ea1, Eb1
		real(dp), dimension(Dim,Dim) 	:: S1, H1, ca1Roothan, cb1Roothan, ca1, cb1, Fa1Roothan, Fb1Roothan, Pa11, Pb11, Pa12, Pb12, &
										   Hint2_1, Hint1_2

		real(dp), dimension(3,Dim)		:: vec2, n2
		real(dp), dimension(Dim)		:: Z2, Ea2, Eb2
		real(dp), dimension(Dim,Dim) 	:: S2, H2, ca2Roothan, cb2Roothan, ca2, cb2, Fa2Roothan, Fb2Roothan, Pa21, Pb21, Pa22, Pb22, HHuckel

		real(dp), dimension(Dim,Dim,Dim,Dim)	:: four_center_tensor1, four_center_tensor2, four_center_tensor_int_12

		real(dp), dimension(3,Dim,Dim)	:: DD1, DD2
		real(dp)						:: tmp1, tmp2, tmpRabs, tmpDipE, tmpNeDipE
		real(dp), dimension(3)			:: tmpRR, tmpDD1, tmpDD2

		logical, dimension(2,size(Z1)) :: excMol1_1,excMol1_2,excMolHOMO,excMolLUMO
		logical, dimension(2,size(Z2)) :: excMol2_1,excMol2_2

		real(dp)	:: parA, parB
		real(dp)	:: targetHOMO_LUMO
		real(dp), dimension(3)	:: targetD

		integer(i4b) 	:: i,j,k,l, c,d,mu,nu, a, dip_ind, phi, x
		real(dp)		:: phi_D, x_D

		write(*,*) 'POKUS4 called'
		! veryfiing transfer to dipole-dipole energy calculation


		call init(1.65_dp,1.65_dp)

		call prepare_carotenoid(vec1,n1,Z1,Dim)
		call prepare_carotenoid(vec2,n2,Z2,Dim)
!		call prepare_alkene(vec1,n1,Z1,Dim)
!		call prepare_alkene(vec2,n2,Z2,Dim)

		call rotate_and_move_molecule(vec2,n2,Z2, 0.0_dp,0.0_dp,0.0_dp,0.0_dp*Angstrom_in_Bohr_radii,0.0_dp,4.1_dp*Angstrom_in_Bohr_radii)

		call fill_four_center_tensor_and_H_S_slater_directed(four_center_tensor1,H1,S1, vec1,n1,Z1,.true.)
		call fill_four_center_tensor_and_H_S_slater_directed(four_center_tensor2,H2,S2, vec2,n2,Z2,.true.)


!!!     ------------- attempt to fir Huckel by Roothan ----------------
!!!		call solve_Roothan_equations(vec1,n1,Z1,S1,H1,four_center_tensor1, ca1Roothan,cb1Roothan,Fa1Roothan,Fb1Roothan,Ea1,Eb1)
!!!		call solve_Roothan_equations(vec2,n2,Z2,S2,H2,four_center_tensor2, ca2Roothan,cb2Roothan,Fa2Roothan,Fb2Roothan,Ea2,Eb2)

		write(*,*) 'VEC1=',vec1
		write(*,*) 'VEC2=',vec2
		!write(*,*) 'N1=',n1
		!write(*,*) 'N2=',n2
		write(*,*) 'Dim', Dim

!		write(*,*) 'S1',S1
!		write(*,*) 'S2',S2
!		write(*,*) 'H1',H1
		write(*,*) 'Ea1',Ea1*hartree_in_invCm

		call fill_HOMO_state(excMolHOMO, Dim)
		call fill_LUMO_state(excMolLUMO, Dim)

		call fill_intermolecular_four_center_tensor_and_H_slater_directed(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12)

		call fill_dipole_operator_elements_slater_directed(vec1,n1,Z1,S1,DD1)
		call fill_dipole_operator_elements_slater_directed(vec2,n2,Z2,S2,DD2)

!!!		call fill_P_matrix(Pa11,Pb11,ca1Roothan,cb1Roothan,excMolHOMO)
!!!		call fill_P_matrix(Pa12,Pb12,ca1Roothan,cb1Roothan,excMolLUMO)
!!!		call fill_P_matrix(Pa21,Pb21,ca2Roothan,cb2Roothan,excMolHOMO)
!!!		call fill_P_matrix(Pa22,Pb22,ca2Roothan,cb2Roothan,excMolLUMO)


!!!		targetHOMO_LUMO = full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa12,Pb12) - full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)
!!!		write(*,*) 'targetHOMO_LUMO',targetHOMO_LUMO*Hartree_in_invCm,'cm^-1'

!!!		do dip_ind=1,3
!!!			targetD(dip_ind) = dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMolHOMO,excMolLUMO,&
!!!							ca1Roothan, cb1Roothan, ca1Roothan, cb1Roothan, DD1(dip_ind,:,:))
!!!		end do

		targetD = 0.0_dp !(/ -1.07, 0.0, 2.85 /)

		if(Dim == 4) then
			targetHOMO_LUMO = 5.7453_dp/Hartree_in_eV
		elseif(Dim == 6) then
			targetHOMO_LUMO = 4.7484_dp/Hartree_in_eV
		elseif(Dim == 8) then
			targetHOMO_LUMO = 4.0878_dp/Hartree_in_eV
		elseif(Dim == 10) then
			targetHOMO_LUMO = 3.3248_dp/Hartree_in_eV
		elseif(Dim == 12) then
			targetHOMO_LUMO = 3.2629_dp/Hartree_in_eV
		elseif(Dim == 14) then
			targetHOMO_LUMO = 2.9864_dp/Hartree_in_eV
		elseif(Dim == 16) then
			targetHOMO_LUMO = 2.7632_dp/Hartree_in_eV
		elseif(Dim == 18) then
			targetHOMO_LUMO = 2.5815_dp/Hartree_in_eV
		elseif(Dim == 24) then
			targetHOMO_LUMO = 20700_dp/hartree_in_invCm
		else
			targetHOMO_LUMO = 0.0_dp
		end if

		!call fit_Huckel_parameters(vec1,n1,Z1,S1, targetHOMO_LUMO,targetD,excMolHOMO,excMolLUMO, parA,parB)

		! old values
		!parA = -11.4_dp/Hartree_in_eV
		!parB = parA*1.75_dp*2.1

		! fit 18-polyenu
		parA = 0.0_dp
		parB = -1.8644567287281373_dp

		! fit 24-carotenoidu
		!parA = 0.0_dp
		!parB =-2.4430311_dp


		call calculate_Huckel_method(vec1,n1,Z1,S1, ca1,HHuckel,Ea1, parA,parB)
		call calculate_Huckel_method(vec1,n1,Z1,S1, cb1,HHuckel,Eb1, parA,parB)
		call calculate_Huckel_method(vec2,n2,Z2,S2, ca2,HHuckel,Ea2, parA,parB)
		call calculate_Huckel_method(vec2,n2,Z2,S2, cb2,HHuckel,Eb2, parA,parB)


!		DO k=1,Dim
!		DO a=1,Dim
		DO k=Dim/2,Dim/2
		DO a=1,Dim/2+1

		call fill_HOMO_state(excMol1_1, Dim)
		call fill_HOMO_state(excMol2_1, Dim)
!		call fill_LUMO_state(excMol1_2, Dim)
!		call fill_LUMO_state(excMol2_2, Dim)

		call fill_HOMO_state(excMol1_2, Dim)
		call fill_HOMO_state(excMol2_2, Dim)

		if(.not.((excMol1_2(1,k) .eqv. .true.) .and. (excMol1_2(1,a) .eqv. .false.))) then
			cycle
		end if

		excMol1_2(1,k) = .false.
		excMol1_2(1,a) = .true.
		excMol2_2(1,k) = .false.
		excMol2_2(1,a) = .true.

		write(*,*) excMol1_1,k,a
		write(*,*) excMol1_2,k,a
		write(*,*) excMol2_1,k,a
		write(*,*) excMol2_2,k,a

		call fill_P_matrix(Pa11,Pb11,ca1,cb1,excMol1_1)
		call fill_P_matrix(Pa12,Pb12,ca1,cb1,excMol1_2)
		call fill_P_matrix(Pa21,Pb21,ca2,cb2,excMol2_1)
		call fill_P_matrix(Pa22,Pb22,ca2,cb2,excMol2_2)

		write(*,*) 'trace P----', trace(matmul(Pa11+Pb11,S1)),trace(matmul(Pa12+Pb12,S1)),trace(matmul(Pa21+Pb21,S2)),trace(matmul(Pa22+Pb22,S2))
		write(*,*) 'HR-----------------', H1
		write(*,*) 'Hü-----------------', HHuckel
		write(*,*) '----'
		write(*,*) 'ca1----', ca1(:,k)
		write(*,*) 'ca2----', ca1(:,a)
		write(*,*)  '----'
		write(*,*) 'ca----', ca1(:,1)
		write(*,*) 'ca----', ca1(:,2)
		write(*,*) 'ca----', ca1(:,3)
		write(*,*) 'ca----', ca1(:,4)
		write(*,*) 'ca----', ca1(:,5)
		write(*,*) 'ca----', ca1(:,6)
		write(*,*)  '----'

		write(*,*) 'parA----', parA
		write(*,*) 'parB----', parB
		write(*,*) 'Huckel energiesA----', Ea1
		write(*,*) 'Huckel energiesB----', Eb1

		write(*,*) 'E1-L - E1-H', full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa12,Pb12)*hartree_in_invCm - full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa11,Pb11)*hartree_in_invCm,'cm^-1  ', (full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa12,Pb12) - full_energy_Huckel(vec1,n1,Z1,HHuckel,Pa11,Pb11))*hartree_in_eV,' eV'
		write(*,*) 'E2-L - E2-H', full_energy_Huckel(vec2,n2,Z2,HHuckel,Pa22,Pb22)*hartree_in_invCm - full_energy_Huckel(vec2,n2,Z2,HHuckel,Pa21,Pb21)*hartree_in_invCm,'cm^-1'
		write(*,*) 'RE1-L - RE1-H', full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa12,Pb12)*hartree_in_invCm - full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)*hartree_in_invCm,'cm^-1'
		write(*,*) 'RE2-L - RE2-H', full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa22,Pb22)*hartree_in_invCm - full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa21,Pb21)*hartree_in_invCm,'cm^-1'

!		! test only
!		do dip_ind = 1,3
!		do a=1, size(Z1)
!			DDnuc1(dip_ind) = DDnuc1(dip_ind) + vec1(dip_ind,a)/size(Z1)
!		end do
!		do a=1, size(Z2)
!			DDnuc2(dip_ind) = DDnuc2(dip_ind) + vec2(dip_ind,a)/size(Z2)
!		end do
!		end do

		call sleep(1)
		return

open(unit=112,file='/home/olsij4am/prace/nose-carotenoidJ-Xphi.dat')
DO phi=0,0
DO x = 0,32

phi_D = phi/11.0_dp*PI
x_D = x*2.0_dp*Angstrom_in_Bohr_radii

call prepare_carotenoid(vec1,n1,Z1,Dim)
call prepare_carotenoid(vec2,n2,Z2,Dim)

call rotate_and_move_molecule(vec2,n2,Z2, phi_D,0.0_dp,0.0_dp,x_D,0.0_dp,4.1_dp*Angstrom_in_Bohr_radii)

call fill_intermolecular_four_center_tensor_and_H_slater_directed(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12)

call fill_dipole_operator_elements_slater_directed(vec1,n1,Z1,S1,DD1)
call fill_dipole_operator_elements_slater_directed(vec2,n2,Z2,S2,DD2)

write(*,*) 'VEC1=',vec1
write(*,*) 'VEC2=',vec2

		tmpNeDipE = interaction_energy(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,excMol2_1,excMol2_2, &
							four_center_tensor1, H1, four_center_tensor2, H2, &
							ca1, cb1, ca1, cb1, ca2, cb2, ca2, cb2,&
							four_center_tensor_int_12, Hint1_2, Hint2_1)
		write(*,*) 'interaction energy:', tmpNeDipE,'=', tmpNeDipE*hartree_in_invCm,'cm^-1'

		tmpRR = 0.0_dp
		do i=1,Dim
			tmpRR = tmpRR + (vec1(:,i)-vec2(:,i))/Dim
		end do

		write(*,*)
		write(*,*) 'R=', sqrt(dot_product(tmpRR,tmpRR))
		write(*,*) 'v1=', vec1(:,1)
		write(*,*) "v1'=", vec1(:,2)
		write(*,*) 'v2=', vec2(:,1)
		write(*,*) "v2'=", vec2(:,2)

		do dip_ind=1,3
			tmpDD1(dip_ind) = dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMol1_1,excMol1_2,&
							ca1, cb1, ca1, cb1, DD1(dip_ind,:,:))
			tmpDD2(dip_ind) = dipole_moment_between_states(vec2,n2,Z2,vec2,n2,Z2, excMol2_1,excMol2_2,&
							ca2, cb2, ca2, cb2, DD2(dip_ind,:,:))
		end do


		write(*,*) 'd1',tmpDD1,sqrt(dot_product(tmpDD1,tmpDD1))
		write(*,*) 'd2',tmpDD2,sqrt(dot_product(tmpDD2,tmpDD2))

		tmpRabs = sqrt(dot_product(tmpRR,tmpRR))
		tmpDipE = dot_product(tmpDD1,tmpDD2)/tmpRabs**3 - 3*dot_product(tmpRR,tmpDD1)*dot_product(tmpRR,tmpDD2)/tmpRabs**5
		write(*,*)
		write(*,*) 'DIPOLE INTERACTION ENERGY', tmpDipE,'=',tmpDipE*hartree_in_invCm,'cm^-1'

		if(abs(tmpDipE-tmpNeDipE)/abs(tmpDipE) > 0.2) then
			write(*,*) 'ENERGIES DIFFER MORE THAN 20%'
		end if

		write(*,*)
		write(*,*)
		write(112,*) x_D/Angstrom_in_Bohr_radii, phi_D, tmpDipE*hartree_in_invCm, tmpNeDipE*hartree_in_invCm, sign(1.0_dp,dot_product(tmpDD1,tmpDD2))
		call flush(112)
END DO
END DO
close(112)

		END DO
		END DO


	end subroutine pokus4

!	subroutine pokus5()
!		integer(i4b), parameter	:: Dim = 14
!
!		real(dp), dimension(3,Dim)		:: vec1, n1
!		real(dp), dimension(Dim)		:: Z1, Ea1, Eb1
!		real(dp), dimension(Dim,Dim) 	:: S1, H1, ca1, cb1, Fa1, Fb1, Pa11, Pb11, Pa12, Pb12, V1, V2, T,&
!										   Hint2_1, Hint1_2
!
!		real(dp), dimension(3,Dim)		:: vec2, n2
!		real(dp), dimension(Dim)		:: Z2, Ea2, Eb2
!		real(dp), dimension(Dim,Dim) 	:: S2, H2, ca2, cb2, Fa2, Fb2, Pa21, Pb21, Pa22, Pb22
!
!		real(dp), dimension(Dim,Dim,Dim,Dim)	:: four_center_tensor1, four_center_tensor2, four_center_tensor_int_12
!
!		real(dp), dimension(3,Dim,Dim)	:: DD1, DD2
!		real(dp)						:: tmp1, tmp2, tmpRabs, tmpDipE, tmpNeDipE
!		real(dp), dimension(3)			:: tmpRR, tmpDD1, tmpDD2
!
!		logical, dimension(2,size(Z1)) :: excMol1_1,excMol1_2
!		logical, dimension(2,size(Z2)) :: excMol2_1,excMol2_2
!
!		integer(i4b), parameter :: mol1ORB = 3, mol2ORB = 3
!
!		integer(i4b) :: i,j,k,l, c,d,mu,nu, a, dip_ind
!
!		write(*,*) 'POKUS4 called'
!		! veryfiing transfer to dipole-dipole energy calculation
!
!		DO i=1, 4
!
!		call prepare_carotenoid(vec1,n1,Z1,Dim)
!		call prepare_carotenoid(vec2,n2,Z2,Dim)
!		call rotate_and_move_molecule(vec2,n2,Z2, i*PI/4.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp)
!
!		call fill_HOMO_state(excMol1_1, Dim)
!		call fill_HOMO_state(excMol2_1, Dim)
!		call fill_LUMO_state(excMol1_2, Dim)
!		call fill_LUMO_state(excMol2_2, Dim)
!
!		write(*,*) excMol1_1
!		write(*,*) excMol1_2
!		write(*,*) excMol2_1
!		write(*,*) excMol2_2
!		write(*,*) 'VEC1=',vec1
!		write(*,*) 'VEC2=',vec2
!		write(*,*) 'N1=',n1
!		write(*,*) 'N2=',n2
!
!		call init(1.625_dp,1.625_dp)
!
!!		call random_number(tmp1)
!!		call random_number(tmp2)
!!
!!		call init_random_seed()
!!
!!		vec1(1,2) = tmp1/sqrt(tmp1**2 + tmp2**2)*1.4_dp
!!		vec1(2,2) = tmp2/sqrt(tmp1**2 + tmp2**2)*1.4_dp
!!
!!		call random_number(tmp1)
!!		call random_number(tmp2)
!!
!!		vec2(1,2) = tmp1/sqrt(tmp1**2 + tmp2**2)*1.4_dp
!		vec2(2,2) = tmp2/sqrt(tmp1**2 + tmp2**2)*1.4_dp
!
!
!		call fill_four_center_tensor_and_H_S_slater(four_center_tensor1,H1,S1, vec1,n1,Z1,T,V1,V2,mol1ORB)
!		call fill_four_center_tensor_and_H_S_slater(four_center_tensor2,H2,S2, vec2,n2,Z2,T,V1,V2,mol2ORB)
!		call solve_Roothan_equations(vec1,n1,Z1,S1,H1,four_center_tensor1, ca1,cb1,Fa1,Fb1,Ea1,Eb1)
!		call solve_Roothan_equations(vec2,n2,Z2,S2,H2,four_center_tensor2, ca2,cb2,Fa2,Fb2,Ea2,Eb2)
!		call fill_P_matrix(Pa11,Pb11,ca1,cb1,excMol1_1)
!		call fill_P_matrix(Pa12,Pb12,ca1,cb1,excMol1_2)
!		call fill_P_matrix(Pa21,Pb21,ca2,cb2,excMol2_1)
!		call fill_P_matrix(Pa22,Pb22,ca2,cb2,excMol2_2)
!
!		write(*,*)
!		write(*,*) '4centr', four_center_tensor1
!		write(*,*) '4centr', H1
!		write(*,*) '4centr', S1
!		write(*,*)
!		write(*,*) '4centr', four_center_tensor2
!		write(*,*) '4centr', H2
!		write(*,*) '4centr', S2
!
!!		write(*,*)
!!		write(*,*) '   T > ', T
!!		write(*,*) '   V1> ', V1
!!		write(*,*) '   V2> ', V2
!!		write(*,*)
!
!		call fill_intermolecular_four_center_tensor_and_H_slater(vec1,n1,Z1,vec2,n2,Z2,Hint1_2,Hint2_1,four_center_tensor_int_12,&
!			mol1ORB,mol2ORB)
!
!		call fill_dipole_operator_elements_slater(vec1,n1,Z1,S1,DD1,mol1ORB)
!		call fill_dipole_operator_elements_slater(vec2,n2,Z2,S2,DD2,mol2ORB)
!
!		write(*,*) 'ca1', ca1
!		write(*,*) 'cb1', cb1
!		write(*,*) 'ca2', ca2
!		write(*,*) 'cb2', cb2
!
!		write(*,*)
!		write(*,*) 'molecule 1 solved:'
!		write(*,*) '----------E1>',full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa11,Pb11)
!		write(*,*) '----------E2>',full_energy(vec1,n1,Z1,four_center_tensor1,H1,Pa12,Pb12)
!!		write(*,*) '----------tD>',dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMol1_1,excMol1_2,&
!!							ca1, cb1, ca1, cb1, DD1)
!!		write(*,*) '----------pD>',dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMol1_1,excMol1_1,&
!!							ca1, cb1, ca1, cb1, DD1)
!		write(*,*) 'molecule 2 solved:'
!		write(*,*) '----------E1>',full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa21,Pb21)
!		write(*,*) '----------E2>',full_energy(vec2,n2,Z2,four_center_tensor2,H2,Pa22,Pb22)
!!		write(*,*) '----------tD>',dipole_moment_between_states(vec2,n2,Z2,vec2,n2,Z2, excMol2_1,excMol2_2,&
!!							ca2, cb2, ca2, cb2, DD2)
!!		write(*,*) '----------pD>',dipole_moment_between_states(vec2,n2,Z2,vec2,n2,Z2, excMol2_1,excMol2_1,&
!!							ca2, cb2, ca2, cb2, DD2)
!
!
!		tmpNeDipE = interaction_energy(vec1,n1,Z1,vec2,n2,Z2, excMol1_1,excMol1_2,excMol2_1,excMol2_2, &
!							four_center_tensor1, H1, four_center_tensor2, H2, &
!							ca1, cb1, ca1, cb1, ca2, cb2, ca2, cb2,&
!							four_center_tensor_int_12, Hint1_2, Hint2_1)
!		write(*,*) 'interaction energy:'
!		write(*,*) '---------->',1.4,tmpNeDipE
!
!
!
!
!		tmpRR = (vec1(:,1) + vec1(:,2) - vec2(:,1) - vec2(:,2))/2
!		write(*,*)
!		write(*,*) 'R=', sqrt(dot_product(tmpRR,tmpRR))
!		write(*,*) 'v1=', vec1(:,1)
!		write(*,*) "v1'=", vec1(:,2)
!		write(*,*) 'v2=', vec2(:,1)
!		write(*,*) "v2'=", vec2(:,2)
!
!		write(*,*) 'DD1',DD1
!		write(*,*) 'DD2',DD2
!
!		do dip_ind=1,3
!			tmpDD1(dip_ind) = dipole_moment_between_states(vec1,n1,Z1,vec1,n1,Z1, excMol1_1,excMol1_2,&
!							ca1, cb1, ca1, cb1, DD1(dip_ind,:,:))
!			tmpDD2(dip_ind) = dipole_moment_between_states(vec2,n2,Z2,vec2,n2,Z2, excMol2_1,excMol2_2,&
!							ca2, cb2, ca2, cb2, DD2(dip_ind,:,:))
!		end do
!
!
!		write(*,*) 'd1',tmpDD1,sqrt(dot_product(tmpDD1,tmpDD1))
!		write(*,*) 'd2',tmpDD2,sqrt(dot_product(tmpDD2,tmpDD2))
!
!		tmpRabs = sqrt(dot_product(tmpRR,tmpRR))
!		tmpDipE = dot_product(tmpDD1,tmpDD2)/tmpRabs**3 - 3*dot_product(tmpRR,tmpDD1)*dot_product(tmpRR,tmpDD2)/tmpRabs**5
!		write(*,*)
!		write(*,*) 'DIPOLE INTERACTION ENERGY', tmpDipE
!
!		if(abs(tmpDipE-tmpNeDipE)/abs(tmpDipE) > 0.2) then
!			write(*,*) 'ENERGIES DIFFER MORE THAN 20%'
!		end if
!
!		write(*,*)
!		write(*,*)
!
!		END DO
!
!
!	end subroutine pokus5


 end module qch_lib
