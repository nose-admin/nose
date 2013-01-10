!
!
!
!
#include "util_allocation.h"

#define SUPERINDEX_FROM_K_L(k, l, lmax) ((l) + (lmax)*((k)-1))
#define K_FROM_SUPERINDEX(superindex, lmax) (((superindex) - 1) / (lmax) + 1)
#define L_FROM_SUPERINDEX(superindex, lmax) (mod(((superindex)-1) , (lmax)) + 1)

#define PALLOCATABLE pointer
!allocatable
!pointer
#define PALLOCATED   associated
!allocated
!associated
#define PNULLIFY(x) nullify(x)
! leave empty if not using pointers
!
!
!
!
module qme_excitonic_networks

	use prepare
	use twod

	use resources_qme	! includes temp, dt
	use std_types
	use nakajima_zwanzig_shared

	use numer_ode
	use numer_fft
	use numer_matrix
	use sci_misc

	use util_allocation

	implicit none

	! declaration
	real, parameter :: kB    = 8.6173e-5_dp		! Boltzmann constant [eV K-1]
	real, parameter :: cm2eV = 1.d0/8065.541_dp	! conversion from cm-1 to eV = hc [eV cm]

	! to be read as fct argument later on
	real :: j0in  =   0.0003_dp
	real :: wDin  = 100._dp			! Debye frequency [cm-1]
!	real :: Kfin  =   1.0d0/1000._dp	! Feeding rate
!	real :: Kdin  =   1.0d0/1000._dp  ! Draining rate

	logical :: reset_in = .false.	! reset feeding if true
! logical :: reset_in = .true.	! reset feeding if true

!	integer, parameter 	:: Nfin=1, Nexc=4, Nsrc=2, Nacp=2	! positions of the final (RC: 1), excited (from the wire, BChl 3: 4), source (Chlorosome: 9) and accepting (to the wire, BChl 1,6: 2,7) sites
!	integer, parameter 	:: Nfin=1, Nexc=4, Nsrc=9, Nacp=2	! positions of the final (RC: 1), excited (from the wire, BChl 3: 4), source (Chlorosome: 9) and accepting (to the wire, BChl 1,6: 2,7) sites
	integer 	:: Nfin, Nexc, Nsrc, Nacp	! positions of the final (RC: 1), excited (from the wire, BChl 3: 4), source (Chlorosome: 9) and accepting (to the wire, BChl 1,6: 2,7) sites
	integer 		::ind									! loop indice

	! test variables
	integer ::ind1, ind2
	complex(dpc) :: maxi

	! function declaration
	private :: hamiltonian			! Prepare the hamiltonian in the site basis
	private :: relaxation_dynamics	! Set the relaxation, feeding matrix
	private :: reversible_dynamics	! Set the reversible dynamics matrix
	private :: set_Superop			! Set superoperator to define the feeding matrix
	private :: init_pops			! Set the initial population matrix
	private :: delta				! Kronecker's delta function
	private :: m2v					! Transform matrix into vector 
	private :: v2m					! Transform vector into matrix 
	private :: SpD					! spectral density function
	private :: nthermal				! Bose-Einstein distribution function
	private :: thermal				! Thermalization function
	private :: w2J					! w^2.J(w) model

	contains



	subroutine main_excitonic_networks()

!*** apo form ***
		integer, parameter	:: Nlev   = 9		! number of levels
!*** holo form ***
!		integer, parameter	:: Nlev   = 10		! number of levels
		integer, parameter	:: Tini   = 0		! initial time step
		integer, parameter	:: Tsteps = 1000	! number of time step

		integer			:: Ntr1, Ntr2    		! dummy variables
		integer			:: tdt
		integer			:: aa, bb,kk, ii		! indices
		integer			:: amin, amax			! border range of indices

		complex(dpc), dimension(Nlev,Nlev)				:: HH			! Hamiltonian matrix
		complex(dpc), dimension(Nlev,Nlev)				:: En, SS, S1	! eigenvalue, eigenvectorsMatrix and inversed eigenvectorMatrix of the hamiltonian
		complex(dpc), dimension(Nlev,Nlev)				:: cc			! correlation between sites
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: RR			! Relaxation matrix
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: Rfeed		! Feeding matrix
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: LL_0, LL		! Reversible dynamics matrix
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: LD, UU, U1	! eigenvalue, eigenvectorsMatrix and inversed eigenvectorMatrix of the reversible dynamics matrix
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: ZZ, Z1		! superoperator U1xUU and inverse matrix
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: Udt			! evolution superoperator
		complex(dpc), dimension(Nlev,Nlev)				:: rho_m_0		! Initial population matrix
		complex(dpc), dimension(Nlev,Nlev)				:: rho_m		! Population matrix at one time step
		complex(dpc), dimension(Nlev,Nlev, Tsteps)		:: rho_m_t		! Population matrix at every time step
		complex(dpc), dimension(Nlev*Nlev)				:: rho_v		! Population matrix at one time step written in the vectorial form

		complex(dpc), dimension(Tsteps)			:: xx, st	! store the evolution of the amount of coherence and of populations
		complex(dpc), dimension(Tsteps+1)			:: std		! store time evolution of the populations, including initial condition
		complex(dpc), dimension(Tsteps+1,Nlev)	:: pop		! Store time evolution of populations, including intinial condition
		complex(dpc), dimension(Tsteps+1,Nlev*(Nlev-1)/2)	:: coh		! Store time evolution of coherence (off-diagonal term of rho_m_t)
		
		! temp declaration for testing
		complex(dpc), dimension(Nlev*Nlev,Nlev*Nlev)	:: TT

		! Variables to compute the von Neumann entropy and entanglement according to Sarovar_NP_6_2010_462
		complex(dpc), dimension(Nlev,Nlev)				:: rho_d		! Diagonal rho and then ln(rho_d), rho ln (rho_d)
		complex(dpc), dimension(Nlev,Nlev)				:: VV, V1		! eigenvector matrix of rho and inverse matrix
		complex(dpc), dimension(Tsteps)					:: entgl		! measure of entanglement according to Sarovar_NP_6_2010_462
		

		real(dp) , dimension(Tsteps)	:: tmem				! time evolution
		
		real(dp) :: j0
		real(dp) :: wD	! Debye frequency [cm-1]
		real(dp) :: Kf	! Feeding rate
		real(dp) :: Kd	! Draining rate

		logical :: reset
		logical :: loc	! true : local basis

		! Variable for non-secular effects
		real(dp) :: rho_osc, omega	! parameters of the oscillatory part of the source term
		rho_osc = 0.0
		omega   = 0.03


		!Initialization of vectors and matrices to avoid problem of memory allocation
		HH = (0.,0.)  	 
		En = (0.,0.)
		SS = (0.,0.)
		S1 = (0.,0.)
		cc = (0.,0.)  	 
		RR = (0.,0.) 	 
		Rfeed = (0.,0.)	 
		LL_0 = (0.,0.)
		LL = (0.,0.)
		LD = (0.,0.)
		UU = (0.,0.)
		U1 = (0.,0.)
		ZZ = (0.,0.)
		Z1 = (0.,0.)
		Udt = (0.,0.) 	 
		rho_m_0 = (0.,0.)  
		rho_m = (0.,0.)	 
		rho_m_t = (0.,0.)  
		rho_v = (0.,0.)	  
		rho_d = (0.,0.)	 
		VV = (0.,0.)
		V1 = (0.,0.)
		entgl = 0.
		
		Nfin = Npos(1)
		Nexc = Npos(2)
		Nsrc = Npos(3)
		Nacp = Npos(4)
						
		if (locBasis == 'yes') then
			loc = .true.
		else
			loc = .false.
		end if
		
		write(6,*) 'Nfin ', 'Nexc ', 'Nsrc ', 'Nacp ', Nfin, Nexc, Nsrc, Nacp, ' loc = ', loc

		reset = reset_in

		! Feeding and draining rates
		Kf = Kfin
		Kd = Kdin

		! bath parameters
		wD = wDin*cm2eV
		j0 = j0in

!===
!---	Define the Hamiltonian, and solve the eigenenergies and eigenstates problem
!===
		call hamiltonian(HH,Nlev)
		
		write(6,*) '********* real(HH) **********'
		write(6,20) transpose(real(HH))

		! Solve eigenvalues and eigenstates problem of the Hamiltonian
		call spec2(HH,SS,En)
		call inv(SS, S1)

		write(6,*) '********* Re eigenvalueMatrix**********'
		write(6,20) transpose(real(En))
		write(6,*) '********* eigenvectors SS**********'
		write(6,20) transpose(real(SS))
		write(6,*) '********* InvEigenvectors S1**********'
		write(6,20) transpose(real(S1))

		write(6,*) '********* SS * S1**********'
		write(6,10) transpose(matmul(SS,S1))
		write(6,*) '*******SS * HH*******'
		write(6,10) transpose(matmul(SS,HH))
		write(6,*) '*******HH* SS*******'
		write(6,10) transpose(matmul(HH,SS))
		write(6,*) '*******SS* En*******'
		write(6,10) transpose(matmul(SS,En))
		write(6,*) '*******En* SS*******'
		write(6,10) transpose(matmul(En,SS))
!===
!---	Define the relaxation matrix in the excitonic basis
!===

!		Correlation between sites
		cc = (0., 0.)
		forall(aa=1:Nlev) cc(aa,aa)=(1., 0.)
		
!		Define the relaxation matrix (Redfield tensor with secular approximation)
		if (doRelax == 'yes') then
			call relaxation_dynamics(RR, Nlev, En, SS, cc)
		else
			RR = (0. ,0.)
		endif


!		write(6,*) '********* Re(RR)  **********'
!		write(6,200) transpose(real(RR))
!		write(6,*) '********* Im(RR)  **********'
!		write(6,200) transpose(imag(RR))
  
!===
!---	Define the feeding matrix
!===

!		Definition in the local basis
		Rfeed = (0., 0.)

!		Feeding from the source pigment
! Warning: if the source pigment is within the molecure complex, feeding of the coherence terms should be added
      	call set_SuperOp(Rfeed,Nlev,Nacp,Nacp,Nsrc,Nsrc,(-1._dp*Kf))
		call set_SuperOp(Rfeed,Nlev,Nsrc,Nsrc,Nsrc,Nsrc,Kf)
		
!		Draining from the wire to the RC
		call set_SuperOp(Rfeed,Nlev,Nfin,Nfin,Nexc,Nexc,(-1._dp*Kd))
		call set_SuperOp(Rfeed,Nlev,Nexc,Nexc,Nexc,Nexc,Kd)
	
		if (Nexc.ne.1 .and. Nexc.ne.Nlev) then ! Nexc is one site of the FMO
			do aa=2,Nlev-1	! aa in molecular complex (e.g. FMO), local basis
				if (aa .ne. Nexc) then
					call set_SuperOp(Rfeed,Nlev,Nexc,aa,Nexc,aa,0.5*Kd)
					call set_SuperOp(Rfeed,Nlev,aa,Nexc,aa,Nexc,0.5*Kd)
				end if
			end do
		end if
		
		!Rfeed = (0., 0.)

!		write(6,*) '********* Re(Rfeed) loc basis **********'
!		write(6,200) transpose(real(Rfeed))

!		Prepare transformation into the excitonic basis
		ZZ = (0.,0.)
		forall(aa=1:Nlev, bb=1:Nlev, ii=1:Nlev, kk=1:Nlev)
			ZZ(Nlev*(bb-1)+aa, Nlev*(kk-1)+ii) = SS(aa, ii) * S1(kk, bb)	!ZZ(ai,kb) = SS(ai)S1(kb) = S1(kb)SS(ai)
		end forall

		call inv(ZZ, Z1)
		
!		write(6,*) '********* Re(Z1) loc basis **********'
!		write(6,200) transpose(real(Z1))
!		write(6,*) '********* ZZ **********'
!		write(6,200) transpose(real(ZZ))

!		Transformation of feeding matrix into excitonic basis
		Rfeed = matmul(Z1, matmul(Rfeed, ZZ))
!		write(6,*) '********* Re(Rfeed) exc basis **********'
!		write(6,200) transpose(real(Rfeed))

		RR = RR + Rfeed

		write(16,*) '**** RR R****'
		write(16,200)	transpose(real(RR))
		write(16,*) '**** RR I****'
		write(16,200)	transpose(imag(RR))
!===
!---	Define the coherent dynamics matrix in the excitonic basis
!===

		call reversible_dynamics(LL_0,En, Nlev)		! = En = S1*HH*SS = matmul(S1,matmul(HH,SS))

		write(16,*) '********* LL_0  R**********'
		write(16,200) transpose(real(LL_0))
		write(16,*) '********* LL_0  I**********'
		write(16,200) transpose(imag(LL_0))

!		LL_0 = (0.,0.)

		LL = LL_0 - (0., 1.)*RR

!		write(6,*) '********* Im(LL)  **********'
!		write(6,200) transpose(imag(LL))

!		Diagonalization to calculate dynamics

		UU= 0.
		LD= 0.
		call spec2(LL,UU,LD)
		call inv(UU, U1)
		write(16,*) '**** UU R****'
		write(16,200)	transpose(real(UU))
		write(16,*) '**** UU I****'
		write(16,200)	transpose(imag(UU))
		write(16,*) '**** LD R****'
		write(16,200)	transpose(real(LD))
		write(16,*) '**** LD I****'
		write(16,200)	transpose(imag(LD))

!---	Testing the error in diagonalization (temporary test, only on the real values)
		TT = (matmul(U1,matmul(LL,UU))-LD)
		maxi = (0.,0.)
		do ind1=1,Nlev*Nlev
			do ind2=1, Nlev*Nlev
				if (Abs(Real(TT(ind1, ind2))) > Real(maxi))  maxi = Abs(TT(ind1, ind2))
!				if (Abs(Imag(TT(ind1, ind2))) > Imag(maxi))  Imag(maxi) = Abs(Imag(TT(ind1, ind2)))
			end do
		end do
		write(6,*) 'Maxi error during diagonalization of LL into LD = ', maxi

!===
!---	Define the evolution superoperator for step dt
!===

!		Construction in the diagonal basis
		Udt = (0., 0.)
		forall(ind=1:Nlev*Nlev) Udt(ind, ind)=exp(LD(ind,ind) * (0.,-1.)*dt)
		write(16,*) '**** Udt diaR****'
		write(16,200)	transpose(real(Udt))
		write(16,*) '**** Udt diaI****'
		write(16,200)	transpose(imag(Udt))

!		Transforming back

		Udt = matmul(UU, matmul(Udt, U1))
		write(16,*) '**** Udt****'
		write(16,200)	transpose(real(Udt))

!===
!---	Set initial populations in the local basis
!===
		call init_pops(rho_m, Nsrc, Nlev)

		! WARNING: temp, non-secular effect
		!rho_m(Nsrc, Nsrc) = (1+ rho_osc)


		write(6,*) '********* real(rho_m_0) in loc. basis **********'
		write(6,20) transpose(real(rho_m))

!		Save initial conditions in rho_m_0
		if (loc) then			
		! Stay in the local basis
			rho_m_0 = rho_m
		else
		! Change to the excitonic basis
			rho_m_0 = matmul(S1,matmul(rho_m,SS))	! S1*rho_m*SS
		end if
		
!		Change to the excitonic basis to prepare rho for the evolution
		rho_m = matmul(S1,matmul(rho_m,SS))

!		write(6,*) '********* real(rho_m_0) in exc. basis **********'
!		write(6,20) transpose(real(rho_m))

!===
!---	Time loop
!===

		Ntr1  = 1
		Ntr2  = Nlev

		rho_m_t = (0., 0.)

		open(unit=222,file=trim(file_join(out_dir,"acoher.dat")))

		do tdt=1,Tsteps
		
			tmem(tdt)= Tini + dt * tdt

!--- Non-secular effect: oscillating inlet population rho(src, src)
			!	Go back to the local basis
			!rho_m = matmul(SS,matmul(rho_m,S1))
			!rho_m(Nsrc, Nsrc) = (1.0+ rho_osc*cos(omega*tdt))*exp(-Kf*tdt)
			
			! Change to the excitonic basis to prepare for the evolution
			!rho_m = matmul(S1,matmul(rho_m,SS))
!---			
		


			rho_v = m2v(rho_m, Nlev)	! Transform into vector 
			!write(6,*) '******** rho_v: before evolution*********'
			!write(6,1) rho_v
			rho_v = matmul(Udt,rho_v)	! Evolve populations by a step dt (in the excitonic basis)
			!write(6,*) '******** rho_v: after evolution*********'
			!write(6,1) rho_v
			rho_m = v2m(rho_v,Nlev)		! Transforming back into matrix

			! Save rho at current time step in rho_m_t according to the basis
			if (loc) then
				rho_m_t(1:Nlev,1:Nlev,tdt) = matmul(SS,matmul(rho_m,S1))
				! Exciton localized in the central matrix
				amin = 2
				amax = Nlev-1
			else
				rho_m_t(1:Nlev,1:Nlev,tdt) = rho_m
				! Excitons localized in the first terms of the matrix
				amin = 1
				amax = Nlev-2
			end if


! checking: verify that rho(a,b)^2 <= rho(a,a) * rho(b,b)
			do aa = amin, amax
				do bb = amin, amax
					if (aa .ne. bb .and. ( abs(rho_m(aa,bb))**2 /(abs(rho_m(aa,aa)*rho_m(bb,bb))) > 1.00001 ) ) then
						write (6,*) 'warning rho2>rhoxrho', tdt, aa, bb, abs(rho_m(aa,bb))**2/ abs(rho_m(aa,aa)*rho_m(bb,bb)), abs(rho_m(aa,bb)),abs(rho_m(aa,aa)), abs(rho_m(bb,bb))
					end if
				end do
			end do

!			Save total population
			st(tdt) = trace(rho_m_t(:,:,tdt))

!			Compute the amount of coherence in the basis defined by the user, using rho_m_t
!			(rho_m_t is in the basis defined by the user)

			xx(tdt) = (0.0, 0.0)

			do aa = amin, amax
				do bb = amin, amax
					if (aa > bb) then
						if ((abs(rho_m_t(aa,aa,tdt)) .ne. 0.).and.(abs(rho_m_t(bb,bb,tdt)).ne.0.)) then
							xx(tdt) = xx(tdt) + (rho_m_t(aa,aa,tdt)+rho_m_t(bb,bb,tdt)) &
					&			*(abs(rho_m_t(aa,bb,tdt))**2)/(rho_m_t(aa,aa,tdt)*rho_m_t(bb,bb,tdt))
						end if
					end if
				end do
			end do

!			if (mod(tdt,100) == 0) write(6,*) xx(tdt), trace(rho_m_t(amin:amax,amin:amax,tdt))
			if (trace(rho_m_t(amin:amax,amin:amax,tdt)) .ne. 0.) then
				xx(tdt) = xx(tdt)/(trace(rho_m_t(amin:amax,amin:amax,tdt))*(amax-amin))
			endif

!	Compute the amount of entanglement according to Sarovar_NP_6_2010_462
			call spec2(rho_m_t(:,:,tdt), VV, rho_d)
			call inv(VV, V1)
			
!			compute	log2(rho_d) and return to the orignal basis
			do aa=1,Nlev
				if (rho_d(aa,aa).ne.0.) rho_d(aa,aa) = log(rho_d(aa,aa))/log(2.)
			end do
			
			rho_d = matmul(VV,matmul(rho_d, V1))	! ln(rho)
			
!			if (tdt == 1) then
!				write(6,*) '********* ln(rho) in the basis defined by the user **********'
!				write(6,20) transpose(real(rho_d))
!			end if
			
			rho_d = matmul(rho_m_t(:,:,tdt),rho_d)	! rho x ln(rho)
			
!			if (tdt == 1) then
!				write(6,*) '********* rho x ln(rho) **********'
!				write(6,20) transpose(real(rho_d))
!			end if
			
			entgl(tdt) = 0.		
			do aa=amin,amax ! loop only within the FMO complex
				if (rho_m_t(aa,aa,tdt) .ne. 0.) entgl(tdt) = entgl(tdt) + rho_d(aa,aa) - rho_m_t(aa,aa,tdt)*log(rho_m_t(aa,aa,tdt))/log(2.)
!				if (tdt == 1) write (*,*) 'entgl', entgl(tdt)
			end do
			
			if (tdt == 1) write (*,*) 'entgl', entgl(tdt)
			
!	output amount of coherence
!			write (222,22) tmem(tdt), real(xx(tdt)), imag(xx(tdt))
			write (222,22) tmem(tdt), max(real(xx(tdt)),1.e-20), imag(xx(tdt)), real(entgl(tdt)), imag(entgl(tdt))

! 			reset feeding
! WARNING, rho_m_0 is basis dependent, and we want rho_m in the exc basis whatever the value of loc is.
!			in future, we should distinguish if loc or not.
			if (reset) rho_m(Nsrc,Nsrc) = rho_m_0(Nsrc,Nsrc)

		end do

		close(unit=222)

!		write(6,*) '********* rho_m at the end of the time loop **********', tdt
!		write(6,10) transpose((rho_m))
	
!===
!---	Output
!===

!		Output populations

		pop = (0., 0.)
		do aa = 1,Ntr2
			pop(1,aa) = rho_m_0(aa,aa);
			pop(2:(Tsteps+1),aa) = rho_m_t(aa,aa,:);
		end do
		
		std(1) = trace(rho_m_0(:,:))
		std(2:(Tsteps+1)) = st(:)
		
		open(unit=111,file=trim(file_join(out_dir,"populations.dat")))
		write(111, '(13a16)') 'Time_step','Real','Real','Real','Real','Real','Imag','Imag','Imag','Imag','Imag','Re(std)','Im(std)'
		do tdt = 1,Tsteps
!			write (111,11) tmem(tdt), real(pop(tdt,:)), imag(pop(tdt,:)), real(std(tdt)), imag(std(tdt))
			write (111,11) tmem(tdt), max(real(pop(tdt,:)),1.e-20), imag(pop(tdt,:)), real(std(tdt)), imag(std(tdt))
		end do
		close(unit=111)


!   	Output coherence

		coh = (0.,0.)
		kk = 0
		do aa = 1,Ntr2
			do bb = 1,Ntr2
				if (aa>bb) then
					kk = kk + 1
					if ( (aa<Ntr2) .and. (bb>=1) ) then
!					if ( (aa<Ntr2) .and. (bb>1) ) then
						coh(1,kk) = rho_m_0(aa,bb)
						coh(2:(Tsteps+1),kk) = rho_m_t(aa,bb,:)
!						write(6,*) aa, bb, kk, rho_m_0(aa,bb), coh(1,kk)
					end if
				end if
			end do
		end do

		open(unit=333,file=trim(file_join(out_dir,"coherences.dat")))
		do tdt = 1,Tsteps
			write (333,33) tmem(tdt), real(coh(tdt,:)), imag(coh(tdt,:)), abs(coh(tdt,:))
!			write (333,33) tmem(tdt), max(1.e-20,real(coh(tdt,:))), imag(coh(tdt,:)), abs(coh(tdt,:))
		end do
		close(unit=333)
		

!*** apo form ***
! for Nlev = 9
   10		format(9(f8.4," +",f6.2,"i,"))	! print complex format for dim(matrix)=N
   20		format(9f10.4)			! print real/imag format for dim(matrix)=N
  200		format(81e9.1)			! print real/imag format for dim(matrix)=N*N

   11		format(f16.2, 9e16.4, 9e16.4, 9e16.4)		! output file format
   22		format(f10.2, 4e20.6)						! output file format
   33		format(f10.2,3(36e16.4))					! output file format
   44		format(4i6,2e16.4)							! output file format

!*** holo form ***
! for Nlev = 10
!    1		format(1(f9.2," +",f9.2,"i,"))	! print complex format for vectors
!    2		format(1f9.2)			! print real/imag format for vectors
!   10		format(10(f8.4," +",f6.2,"i,"))	! print complex format for dim(matrix)=N
!   20		format(10f8.4)			! print real/imag format for dim(matrix)=N
!  100		format(100(f5.2," +",f5.2,"i,"))	! print complex format for dim(matrix)=N*N
!  200		format(100e9.1)			! print real/imag format for dim(matrix)=N*N

!   11		format(f16.2, 10e16.4, 10e16.4, 10e16.4)		! output file format
!   22		format(f10.2, 4e20.6)		! output file format
!   33		format(f10.2,3(45e16.4))		! output file format

	write(6,*) Kf, Kd, loc, doRelax, Npos


	end subroutine main_excitonic_networks


	subroutine hamiltonian(H,N)
		
		complex(dpc), dimension(N,N), intent(out)	:: H
		integer							:: N
		integer							:: i, j, k
		complex(dpc)					:: jc	! max coupling between sites
!*** apo form ***
		complex(dpc), dimension(9)	:: ee	! WARNING, manually fixed dimension. Allocatable instead of dimension(N) can lead to segmentation default during compilation
		real(dp), dimension(21)		:: je	! WARNING, manually fixed dimension. coupling distribution between the exciton dim= ((N-2)^2 - (N-2))/2
!*** holo form ***
!		complex(dpc), dimension(10)	:: ee	! WARNING, manually fixed dimension. Allocatable instead of dimension(N) can lead to segmentation default during compilation
!		real(dp), dimension(28)		:: je	! WARNING, manually fixed dimension. coupling distribution between the exciton dim= ((N-2)^2 - (N-2))/2


		!N = size(H,1)
		! initialization of the hamiltonian
		H = (0.,0.)

		! energies of the sites
		! Qy absorption of the chlorosome is at 12 500 cm-1

!*** apo form ***
! FMO with 7 BChls: Energies and interactions from Adolphs_BJ_2006_2778
!		! Energie shift of E3 = 12 210 cm -1
		ee = (/ (-290.,0.), (200.,0.), (320.,0.), (0.0,0.), (110.0,0.), (270.0,0.), (420.0,0.), (230.0,0.), (290.,0.)/)*cm2eV

 		je = (/ (-87.80,0.), &
 		&		(5.5,0.), (30.8,0.), &
 		&		(-5.9,0.), (8.2,0.), (-53.5,0.), &
 		&		(6.7,0.), (0.7,0.), (-2.2,0.), (-70.7,0.), &
 		&		(-13.7,0.), (11.8,0.), (-9.6,0.), (-17.0,0.), (81.1,0.), &
 		&		(-9.9,0.), (4.3,0.), (6.0,0.), (-63.3,0.), (-1.3,0.), (39.7,0.)/)*cm2eV
!		je = (/ (-00.00,0.), &
!		&		(0.0,0.), (30.8,0.), &
!		&		(-0.0,0.), (8.2,0.), (-53.5,0.), &
!		&		(0.0,0.), (0.7,0.), (-2.2,0.), (-70.7,0.), &
!		&		(-00.0,0.), (11.8,0.), (-9.6,0.), (-17.0,0.), (81.1,0.), &
!		&		(-0.0,0.), (4.3,0.), (6.0,0.), (-63.3,0.), (-1.3,0.), (39.7,0.)/)*cm2eV

!*** holo form ***
! FMO with 8 BChls: Energies and interactions from Olbrich_JPC_2011_8609
! Average values for the trimer
! in eV (Table 1)
! Energie shift of E3 = 1 493 meV
!		ee = (/ (-56.8,0.), (23.0,0.), (10.0,0.), (0.0,0.), (14.0,0.), (8.0,0.), (11.0,0.), (61.0,0.), (27.0,0.), (56.8,0.)/)*1.e-3 ! from meV to eV
! in cm-1
!		ee = (/ (-290.,0.), (186.,0.), (081.,0.), (0.0,0.), (113.0,0.), (065.0,0.), (089.0,0.), (492.0,0.), (218.0,0.), (290.,0.)/)*cm2eV

! Energie shift of E3 = 12 042 cm-1
!		je = (/ (-80.3,0.), &
!		&		(3.5,0.), (23.5,0.), &
!		&		(-4.0,0.), (6.7,0.), (-49.8,0.), &
!		&		(4.5,0.), (0.5,0.), (-1.5,0.), (-63.4,0.), &
!		&		(-10.2,0.), (07.5,0.), (-6.5,0.), (-13.3,0.), (55.8,0.), &
!		&		(-4.9,0.), (1.5,0.), (1.2,0.), (-42.2,0.), (4.7,0.), (33.0,0.), &
!		&		(21.0,0.), (3.3,0.), (0.7,0.), (-1.2,0.), (2.8,0.), (-7.3,0.), (-8.7,0.) /)*cm2ev


		forall(i = 1:N)   H(i,i) = ee(i)	! initialization of the diagonal elements

		jc = (1., 0.)

		k = 0
		do i = 2,(N-1)
			do j = 2,(N-2)
				if (i > j) then
					k = k + 1
					H(i,j) = jc * je(k)
					H(j,i) = H(i,j)
				end if
			end do
		end do

	end subroutine hamiltonian



	subroutine relaxation_dynamics(R, N, E, SS, cc)

		complex (dpc), dimension(N*N,N*N), intent(inout)	:: R
		complex(dpc), dimension(N,N), intent(in)			:: SS
		complex(dpc), dimension(N,N), intent(in)			:: E, cc
		integer, intent(in)								:: N
		
		complex(dpc), dimension(N,N)		:: S1
		complex(dpc), dimension(N,N,N)	:: KK
		complex(dpc), dimension(N,N)		:: kl
		complex(dpc), dimension(N)		:: gg
		integer								:: aa,bb,mm,nn
		complex(dpc)						:: sum


		call inv(SS, S1)
		
		KK = (0.,0.)
		do aa = 1,N
			do bb = 1,N
				do mm = 1,N
					KK(aa,bb,mm) = S1(aa,mm)*SS(mm,bb)
!					write(6,'(3i5,2f5.2)') aa, bb, mm, KK(aa,bb,mm)
  				end do
			end do
		end do

!		population transfer rates (kab := a->b)
		kl = (0.,0.)
	    do aa=1,N
			do bb=1,N
				if (aa.ne.bb) then
           			do mm = 1,N 
               			do nn = 1,N
                   			kl(aa,bb) = kl(aa,bb) + SpD(E(aa,aa)-E(bb,bb))*cc(mm,nn)*KK(aa,bb,mm)*KK(bb,aa,nn)
               			end do
           			end do
				end if
			end do
        end do


!		Verification of kl(a,b): check that kl(a,b) = kl(b,a) * exp(h wab / kb T)  (Eq. 3.290)
		do aa = 1,N
			do bb = 1,N
				if ( abs(kl(aa,bb)) .ne. abs( kl(bb,aa)*thermal(E(aa,aa)-E(bb,bb)) ) ) then
					write(*,'(a30,5e15.3)') 'check detailed balance (=1)',abs(kl(aa,bb))/abs(kl(bb,aa)*thermal(E(aa,aa)-E(bb,bb))), real(E(aa,aa)-E(bb,bb)), real(kl(aa,bb)), real(kl(bb,aa)), real(kl(bb,aa)*thermal(E(aa,aa)-E(bb,bb)))
				end if
			end do
		end do

		!write(6,*) KK
		!write(6,*) '***********kl**************'
		!write(6,*) kl
   
!		depopulation rates
		do aa = 1,N
			sum = (0.0, 0.0)
			do bb = 1,N
!				sum = sum + kl(bb,aa)	! original matlab code
				sum = sum + kl(aa,bb)	! cf. Eq. 3.286 [MayKuehnBook] corrected on 2012.05.16
			end do
			kl(aa,aa) = sum
		end do
		


!		population dynamics

		forall( aa=1:N ) R(N*(aa-1)+aa,N*(aa-1)+aa) = kl(aa,aa)
!		forall( aa=1:N, bb=1:N, aa .ne. bb) R(N*(aa-1)+aa,N*(bb-1)+bb) = -kl(aa,bb) ! original matlab code
		forall( aa=1:N, bb=1:N, aa .ne. bb) R(N*(aa-1)+aa,N*(bb-1)+bb) = -kl(bb,aa) ! cf. Eq. 3.286 [MayKuehnBook] corrected on 2012.05.16

!		dephasing rates

		gg = (0., 0.)
		forall (aa = 1:N) 
			gg(aa) = 0.5*kl(aa,aa)
		end forall
    
!		dephasing dynamics

		forall( aa=1:N, bb=1:N, aa.ne.bb)  
                R(N*(aa-1)+bb,N*(aa-1)+bb) = gg(aa) + gg(bb)
        end forall


!		write(6,*) '-------Test RR-------'
!		do bb = 1,N
!			sum = 0.
!			do aa = 1,N
!				sum = sum + R(N*(aa-1)+aa,N*(bb-1)+bb)
!			end do
!			write(6,*) sum
!		end do


!   	Output relaxation rates
		open(unit=444, file= trim(file_join(out_dir,"relaxation.dat")))		
		write(444,'(4a5, 2a15)') 'aa', 'bb', 'cc', 'dd', 'R', 'kl'
		do aa=1,N
			do bb=1, N
 				do mm=1, N
 					do nn=1, N
						write(444,'(4i5, 2e15.3)') aa, bb, mm, nn, real(R(N*(aa-1)+bb, N*(mm-1)+nn)), real(kl(aa,bb))
 					end do
 				end do
!			write(444,'(4i5, 2e15.3)') aa, bb, mm, nn, real(R(N*(aa-1)+bb, N*(aa-1)+bb)), real(kl(aa,bb))
			end do
		end do
		close(unit=444)
		
	end subroutine relaxation_dynamics


	subroutine set_Superop(Rf, N, Na, Nb, Nc, Nd, val)
		complex(dpc), dimension(N*N, N*N), intent(inout) :: Rf
		integer	, intent(in)							:: Na, Nb, Nc, Nd
		integer	, intent(in)							:: N
		!real(dp)										:: N2
		real(dp), intent(in)							:: val

		!N2 = size(Rf,1)*1.
		!N2 = sqrt(N2)
		!N = int(N2)

		Rf(N*(Na-1)+Nb,N*(Nc-1)+Nd) = val
		write(6,'(a15, 2i4,2e15.2)') 'set_Superop', (N*(Na-1)+Nb), (N*(Nc-1)+Nd), Rf(N*(Na-1)+Nb,N*(Nc-1)+Nd)

	end subroutine


	subroutine reversible_dynamics(LL,En, N)

		complex (dpc), dimension(N*N,N*N), intent(out)	:: LL
		complex (dpc), dimension(N,N), intent(in)		:: En
		integer, intent(in)								:: N
		integer											:: i, j, k, l

		LL = (0.,0.)

		forall(i=1:N, j=1:N, k=1:N, l=1:N) &
		&		LL(N*(i-1)+j,N*(k-1)+l) = En(i,k)*delta(j,l) - delta(k,i)*En(l,j)


	end subroutine reversible_dynamics



	subroutine init_pops(rho, src, N)

		complex (dpc), dimension(N,N), intent(out)	:: rho
		integer	, intent(in)			:: N	! number of levels
		integer	, intent(in)			:: src	! position of the source site
		integer										:: ii, jj

		rho = (0.,0.)

		rho(src, src) = (1,0.)

! 		Completely coherent system
!		forall (ii=2:(N-1), jj=2:(N-1))	 rho(ii, jj) = (0.142857,0.)

	end subroutine init_pops


	pure function m2v(m,N)
!	transform the matrix lines into a vector element
	
		complex(dpc), dimension(N*N)				:: m2v
		complex(dpc), dimension(N,N), intent(in)	:: m
		integer, intent(in)						:: N	! number of levels
		integer										:: i
	
		forall(i=1:N) m2v(N*(i-1)+1:N*i) = m(i,1:N)
	
	end function


	pure function v2m(v,N)
!	transform a vector into lines of a matrix
		
		complex(dpc), dimension(N,N)				:: v2m
		complex(dpc), dimension(N*N), intent(in)	:: v
		integer, intent(in)						:: N	! number of levels
		integer										:: i

		forall(i=1:N) v2m(i,1:N) = v(N*(i-1)+1:N*i)
	
	end function


	pure function delta(i,j)
	
		real					:: delta
		integer, intent(in)	:: i,j

		delta = 0.
		if (i==j) delta = 1.
	
	end function
	
	
	pure function SpD(w)
	
		complex(dpc)				:: SpD, wD
		complex(dpc), intent(in)	:: w
		
		wD = wDin * cm2ev
		if (abs(w) == 0.) then
			SpD = j0in/wD**2
		else
			SpD = ( 1. + nthermal(w))* (w2J(w) - w2J(-1_dpc*w))
		end if
			
	end function


	pure function nthermal(w)
	! Bose-Einstein distribution function
	
		complex(dpc)				:: nthermal
		complex(dpc), intent(in)	:: w

		nthermal = 1._dp / (exp(w/(kB*temp))-1._dp)

	end function


	pure function thermal(w)
	! Bose-Einstein distribution function
	
		complex(dpc)				:: thermal
		complex(dpc), intent(in)	:: w

		thermal = exp(w/(kB*temp))

	end function


	pure function w2J(w)
	! w**2 . J(w) using Debye spectral density
		complex(dpc)				:: w2J, wD
		complex(dpc), intent(in)	:: w

		wD = wDin*cm2ev
		if (real(w) >= 0.) then
			w2J = j0in*w/(w*w + wD*wD)
		else
			w2J = (0., 0.)
		endif

	end function

end module qme_excitonic_networks

