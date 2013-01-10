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
module qme_vibronic_networks

	use prepare
	use twod

	use resources_qme	! includes temp, dt
	use std_types
	use nakajima_zwanzig_shared

	use numer_ode
	use numer_fft
	use numer_matrix
	use sci_misc
	use numer_gamma

	use util_allocation

	implicit none
	
	logical :: oneParticle=.true.	! if true, one-particle approximation: molecules in their electronic ground state are also in their vibrational ground state								
!	Define if we want to compute the whole SPECTRUM or only the diagonal of the spectrum
	logical, parameter:: spectrum=.true.
!-- declaration of parameters
!	real(dp) :: j0in  = 0.0003_dp
	real(dp) :: lambdain= 35._dp			! Bath reorganization energy [cm-1]
	real(dp) :: wDin  = 130._dp			! Debye frequency [cm-1]
 	real(dp) :: HR = 0.05_dp			! Huang-Rys factor (lambdain = HR * wHin)


!--- generate dgap randomly
!    real(dp), dimension(:), allocatable :: ran    ! random number with normal distribution
    real(dp),dimension(1) :: ran    ! random number with normal distribution
    real(dp) :: dgap    ! shift of the gap energy due to static disorder
	
	integer :: iterate	! counter for the disorder loop   
	real(dp):: dgapi    ! shift of the gap energy due to static disorder 
	real(dp):: PDF		! Probability density function of the normal distribution
	real(dp):: PDFs		! sum of the probability density function of the normal distribution (for verification only)

	
!    real(dp) :: ran    ! random number with normal distribution
!    real(dp) :: den    ! energy shift due to static disorder



	real, parameter :: kB    = 8.6173e-5_dp			! Boltzmann constant [eV K-1]
	real, parameter :: cm2eV = 1.d0/8065.541_dp		! conversion from cm-1 to eV = hc [eV cm]
! 	real, parameter :: cm2eV = 1._dp		! conversion from cm-1 to eV = hc [eV cm]

	integer 	:: Nfin, Nexc, Nsrc, Nacp	! positions of the final (RC: 1), excited (from the wire, BChl 3: 4), source (Chlorosome: 9) and accepting (to the wire, BChl 1,6: 2,7) sites
	integer		::ind									! loop indice

!	integer	:: beyond_wH			! 1 if Abs(gap)>wHin. Used to track the vibrational coherences after redistribution of the eigenenergies, when computing the enhancement relative to J=0
!	integer	:: neg_gap				! 1 if gap <0.
	
	logical	:: signal=.true.		! if true: calculate evolution of total signal (4 interactions with light, only non-rephasing pathway)
	                                ! else: calculate coherences (2 interactions with light)
	
	type aggregate									! Characteristics of the molecular aggregate
		
		integer, pointer 					:: N	  !number of molecules
		integer, pointer					:: nu	  !number of vibrational levels
!		real(dp), pointer					:: w0	  !frequency of the vibrational mode (in eV)
		real(dp),pointer					:: theta  !angle between the different molecules
		real(dp),pointer					:: gap	  !initial energy gap in cm-1 (without energetic disorder)
		real(dp),pointer					:: dwidth !standard deviation of the normal distribution for energy disorder (in cm-1)
		real(dp), dimension(:), pointer		:: d	  !coordinate shift from the ground state position
		real(dp), dimension(:), pointer		:: E	  !site energy (in eV)
		real(dp), dimension(:,:), pointer	:: J	  !coupling energy (in eV)
		real(dp), dimension(:,:), pointer	:: ori	  !orientation of the molecules (x,y,z)
	
	end type aggregate
	
	type(aggregate), private, pointer 		:: agg	
	
	integer, private						:: Ntot, Nnu	! Total number of levels, and vib. levels (including source and optical levels)
		
	private :: Hamiltonian		
	private :: matrice_init					! matrix/vector allocation
	private :: matrice_deinit				! matrix/vector deallocation
	private :: aggregate_init
	private :: aggregate_deinit
	private :: SpD0							! Normalized low-frequency Spectral Density
	private :: SpD							! Total Spectral Density (includes the high-frequency mode)
	private :: Crl							! Correlation function
	private :: Gamma_init					! Initialize coefficient matrix for pure dephasing relaxation
	private :: Relaxation					! Compute Redfield matrix
	private :: reversible_dynamics			! compute the matrix of coherence transport
	private :: init_pops					! prepare initial population
	private :: dynamics						! dynamical time loop
	private :: signal2D						! calculate the 2D signal from the results of dynamics
	private :: m2v							! transform a matrix N,N into vector N*N
	private :: v2m							! transform a vector N*N into a matrice N,N
	private :: delta						! delta function
	private :: CrD						! Debye spectral density
	private :: w2J							! 
	private :: FC							! compute Franck-condon factor 
	private :: relax_time					! compute the relaxation time

    private :: gap_disorder

	private :: goftHT						! Compute the line shape function in the high temperature limit (approx. OK for 77K)
	private :: GFT							! Fourier transform of G(w), including the goft (in high temp. limit)

contains


	subroutine main_vibronic_networks()
	

		character(len=15)	:: formN, formNtot, formNN, formNtotNtot, formSpec	! output format
		character(len=20)	:: formPop, formCoh, formSig	!output formats
		character(len=4)	:: dim             	!output formats
		character(len=20)	:: out_file			!name of the output files for the coherences 
		character(len=20)	:: out_file0			!name of the output files for the coherences 

		integer										:: nn, mm, aa, bb, tdt, xi, aam, bbm
		integer										:: nu, mu
		integer										:: opt		! for output only
		integer										:: ww, ww3	! loop over the frequency of observation
		complex(dpc)								:: w1, w3	! frequencies at which we oberve the 2D signal
		complex(dpc)								:: w 				! contains the projectors from loc to vib basis
		complex(dpc), dimension(2000)				:: OD		! linear absorption spectrum (arbitrary dimension)
		complex(dpc)								:: dw		! frequency difference
		type(aggregate)								:: agg
		complex(dpc), dimension(:,:), allocatable	:: HH, En, SS, S1
		complex(dpc), dimension(:,:), allocatable	:: SSavg, Enavg, DMavg	! composition of the eigenvectors, eigenenergies and transient dipole moment averaged over energetic disorder
		complex(dpc), dimension(:,:), allocatable	:: KK				! contains the projectors from loc to vib basis
		complex(dpc), dimension(:,:), allocatable	:: cc				! correlation between site
		complex(dpc), dimension(:,:), allocatable	:: Gamma			! matrix coefficient to transformation the correlation function in the excitonic basis
		complex(dpc), dimension(:,:), allocatable	:: RR, RRD			! relaxation matrix

		complex(dpc), dimension(:,:), allocatable	:: pop		! Store time evolution of populations, including intinial condition
		complex(dpc), dimension(:,:), allocatable	:: coh		! Store time evolution of coherence (off-diagonal term of rho_m_t)
		complex(dpc), dimension(:)  , allocatable	:: std		! Store total population
		real(dp), dimension(:), allocatable			:: tmem			! save time vector

		complex(dpc), dimension(:,:), allocatable	:: Xsi		! Vibrational character of the coherences, integrated over the energy disorder

!		complex(dpc), parameter	:: wini=12000._dpc	! inital frequency for sampling of the signal (cm-1)
!		complex(dpc), parameter :: wdw =020._dpc		! discretization step for w1 (cm-1)
!		integer, parameter		:: wsteps= 20 		! number of iteration on w

		real(dp) :: wHin           							! High-frequency mode treated explicitely [cm-1]
		real(dp), parameter :: wHini   =117._dp			! initial mode frequency [cm-1]
		real(dp), parameter :: wHdw    =50._dp				! incremental step for the mode frequency
		integer , parameter :: wHsteps =  1  				! number of iteration of wHin
		integer				:: wwH			   				! incremental variable used in the loop on wH
!		real(dp) :: lambdain  = 185.*0.05_dp	! bath reorganization energy [cm-1]
		
!----   Prepare dynamics
		integer, parameter	:: Tsteps = 1500		! number of time step
		integer		:: Nfin, Nexc, Nsrc, Nacp	! positions of the final (RC: 1), excited (from the wire, BChl 3: 4), source (Chlorosome: 9) and accepting (to the wire, BChl 1,6: 2,7) sites
		logical		:: loc	! true : local basis

  		complex(dpc), dimension(:,:), allocatable	:: Rfeed	   ! Feeding matrix
		complex(dpc), dimension(:,:), allocatable	:: LL,LL_0 		! Reversible dynamics matrix
		complex(dpc), dimension(:,:), allocatable	:: LD, UU, U1	! eigenvalue, eigenvectorsMatrix and inversed eigenvectorMatrix of the reversible dynamics matrix
		complex(dpc), dimension(:,:), allocatable	:: TT			! used for testing diagonalization of LL
!		complex(dpc), dimension(:,:), allocatable	:: ZZ, Z1	   ! superoperator U1xUU and inverse matrix
		complex(dpc), dimension(:,:), allocatable	:: Udt  	   ! evolution superoperator
		complex(dpc), dimension(:,:), allocatable	:: rho_m_0  	! Initial population matrix
!		complex(dpc), dimension(:,:), allocatable	:: rho_m		! Population matrix at one time step
		complex(dpc), dimension(:,:,:), allocatable :: rho_m_t		! Population matrix at every time step
!		complex(dpc), dimension(:), allocatable 	:: rho_v		! Population matrix at one time step written in the vectorial form

		complex(dpc), dimension(:,:,:), allocatable	:: Amp			! Amplitude of the signal at initial time, after interaction with the pulses, before any evolution
		complex(dpc), dimension(:,:,:,:), allocatable :: Amp_t		! Amplitude of the signal at every time step

		complex(dpc), dimension(:,:), allocatable	:: sig			! For output of Amp_t: signal of all difference pathways
		complex(dpc), dimension(:), allocatable :: totsig			! Sum of the signal over all pathways
		complex(dpc), dimension(:,:,:), allocatable :: spec_tot     ! Total spectrum (w1, w3, t)

!		real(dp) , dimension(Tsteps)	:: tmem				! time evolution

		complex(dpc) :: maxi	! using in testing LL diagonalization
		real(dp), dimension(:), allocatable	:: maxAmp		! Maximum amplitude of the coherence after decay of the electronic one. 
		real(dp)	:: maxSig	! Maximum tot signal, output in spectrum.dat


!--- For output / verification purpose
		complex(dpc), dimension(:), allocatable	:: gft2	! GFT(w1)*GFT(w3) averaged over the energetic disorder
		complex(dpc), dimension(:,:), allocatable	:: sigdE			! For output of Amp_t: signal of all difference pathways
		complex(dpc), dimension(:,:,:), allocatable	:: dummy		! transitoire variable

		real(dp) :: j0
		real(dp) :: wD	! Debye frequency [cm-1]
		real(dp) :: Kf	! Feeding rate
		real(dp) :: Kd	! Draining rate

!---    Define variable from the input file

		Nfin = Npos(1)
		Nexc = Npos(2)
		Nsrc = Npos(3)
		Nacp = Npos(4)

		if (locBasis == 'yes') then
			loc = .true.
		else
			loc = .false.
		end if
!		write(6,*) 'vib: Nfin, Nexc, Nsrc, loc ', Nfin, Nexc, Nsrc, loc
		
		write(6,*) 'Starting calculation for: ', molSystem

		! Feeding and draining rates
		Kf = Kfin
		Kd = Kdin

!---	Test gamma function for complex numbers
!		open(unit=101, file= trim(file_join(out_dir,"cgamma.dat")))
!		write(101,'(7a15)') 'Re(w)', 'Im(w)', 'Re(Gamma(w))', 'Im(Gamma(w))', 'Abs(Gamma(w))',  'Re(IncGamma(w))', 'Im(IncGamma(w))'
!		do aa = -50, 50 	! [cm-1]
!			w = aa*(0.1,0.2)		! complex number
!			call cgamma(0,w,w1)	! w1=gamma(w)
!			write(101,'(7e15.3)') real(w), imag(w), real(w1), imag(w1), abs(w1), real(cincgamma(w,w)), imag(cincgamma(w,w))
!		end do 
!		close (unit=101)
		
!-- Allocatable matrix for the total spectrum
		allocate(spec_tot(spec_fq(3),spec_fq(3),Tsteps+1)) ! spec_fq(3) is the number of steps of w1, w3 in 2D spectrum
		spec_tot=0.0
!--


!---

!  prepare file were the maximum amplitude of the vibrational coherences will be saved 
   open(unit=10, file=trim(file_join(out_dir,"max.dat")))
   write(10,'(4a10)') "#        w","wHin","Max amp.","Xsi(3,4)"
   
!   open(unit=555, file=trim(file_join(out_dir,"gft2.dat")))
!   write(555,'(2a15, 2a25)') '#         w1', 'w3', 'Re(gft(w1)*gft(w3))', 'imag'

   ! open file only for calculation of spectrum
   if (spectrum) then
   		open(unit=11,file=trim(file_join(out_dir,"spectrum.dat")))
   		write(11,'(2a15,22a25)') '#         w1', 'w3', 'Re(Max(Totsig(t>1ps)))', &
			& 'Re(totsig(1))'   , 'Re(totsig(100))' , 'Re(totsig(200))' , 'Re(totsig(300))' , 'Re(totsig(400))' , 'Re(totsig(500))' ,	&	
			& 'Re(totsig(600))' , 'Re(totsig(700))' , 'Re(totsig(800))' , 'Re(totsig(900))' , 'Re(totsig(1000))', 'Re(totsig(1100))',	&	
			& 'Re(totsig(1200))', 'Re(totsig(1300))', 'Re(totsig(1400))', 'Re(totsig(1500))', 'Re(totsig(1600))', 'Re(totsig(1700))',	&	
			& 'Re(totsig(1800))', 'Re(totsig(1900))', 'Re(totsig(2000))'
	end if
	
   open(unit=22, file=trim(file_join(out_dir,"sigdE.dat")))
   if (pathway=='R1') write(22,'(2a15, 7a15)') '#         w1', 'dE', 'Max(sigdE(0))', 't>500', 't>1000', 't>1500', 't>2000', 'Re(gft2)', 'Im(gft2)'
   if (pathway=='R4') write(22,'(2a15, 7a15)') '#         w1', 'dE', 'sigdE (dimer)', 'sigdE (monomer)', 'Re(gft2)', 'Im(gft2)'


   open(unit=33, file=trim(file_join(out_dir,"Ew.dat")))
   write(33,'(2a15, 12a15)') '#         w0', 'w1', 'En1', 'En2', 'En3', 'En4', 'En5', 'En6', &
     &                      'Xsi34', 'Xsi35','SS44', 'SS45','SS54', 'SS55'
   write(33,'(2a15, 16a15)') '#         w0', 'w1', 'Gd,nu0', 'Gd,nu1', 'Exc1', 'Exc2', 'Exc3', 'Exc4', &
     &                      'Xsi12', 'Xsi13','SS exc2,n1nu1', 'SS exc2,n2nu0','SS exc3,n1nu1', 'SS exc3,n2nu0', &
     &                      'SS exc1,n1nu0', 'SS exc1,n1nu1','SS exc2,n1nu0','SS exc3,n1nu0'

!   open(unit=33, file=trim(file_join(out_dir,"Amp.dat")))
!   if (pathway=='R1') write(33,'(10a15)') 'Amp(4,3,1)', 'Amp(3,3,2)', 'DM(3,1)', 'DM(3,2)', 'DM(4,1)', 'DM(4,2)'


!### Loop over the frequency of observation
!	SPECTRUM
  do ww3=1, spec_fq(3) 
  w3 = spec_fq(1) + spec_fq(2)*(ww3-1)
  write(6,'(a20,f10.0)') 	'###### w3 ######', real(w3)

   do ww=1, spec_fq(3) 
    w1 = spec_fq(1) + spec_fq(2)*(ww-1)  ! effective frequency at wich we sample the signal by call of signal2D
	write(6,'(a20,f10.0)') 	'###### w1 ######', real(w1)

   do wwH=1, wHsteps
    wHin = wHini + wHdw*(wwH-1)
	write(6,'(a20,f10.0)') 	'###### wHin ####', real(wHin)
	
	if (.not.spectrum) w3 =w1
 	if (molSystem=='jonas') w3=12185. !E2-Eopt
	
!=== Loop over gap energy disorder. 
	PDF = 0.
	PDFs = 0.
	do iterate=1,Nruns	!run over a representative width of the gaussian
		if (Nruns==1) PDF=1	! 1 single only
		dgapi = dgapini+iterate	 ! minimum should be low enough to cover the normal distribution
!---

		call aggregate_init(agg, wHin)

		!  dimension of the Hamiltonian
		if (oneParticle) then
			Ntot = agg%N * agg%nu		
		else 
			Ntot = agg%N * agg%nu**(agg%N-2)	 ! N-2 excited states
		end if 
		Nnu	 = agg%nu					! save the number of vib. level 

!---	Compute the PDF of the normal distribution to weight results with the corresponding probability
!		First used in init_pops: should be computed before call to this subroutine
		PDF = exp(-(dgapi-agg%gap)**2/(2*agg%dwidth**2))/(Sqrt(2*pi)*agg%dwidth) ! dwidth=sigma, gaussian standard deviation
		if (Nruns==1) PDF=1 !doesn't affect results for only 1 run, with the original energy gap
		PDFs = PDFs + PDF
		if (iterate==Nruns)	write(6,*) 'sum of PDF = ', PDFs , ' after ', Nruns, ' iterations (should be 1)'
!---

		! Prepare output format 
		! for matrices of dimension N (real or imag part)
		write(dim,'(i2)') agg%N
		formN = '('//trim(adjustl(dim))//"f10.4"//')'
		! for matrices of dimension N*nu (real or imag part)
		write(dim,'(i2)') Ntot
		formNtot = '('//trim(adjustl(dim))//"f10.4"//')'
		! for matrices of dimension N*N (real or imag part)
		write(dim,'(i2)') agg%N*agg%N
		formNN = '('//trim(adjustl(dim))//"f10.4"//')'
		! for matrices of dimension (N*nu)*(N*nu) (real or imag part)
		write(dim,'(i3)') Ntot*Ntot
		formNtotNtot = '('//trim(adjustl(dim))//"e9.1"//')'
		write(6,*) 'Format = formN, formNtot, formNN, formNtotNtot ', formN, formNtot, formNN, formNtotNtot

!		if (iterate==1) then	
!---    Allocate vectors and matrices	
			call matrice_init(HH,Ntot)
			call matrice_init(SS,Ntot)
			call matrice_init(S1,Ntot)
			call matrice_init(En,Ntot)

			call matrice_init(cc,Ntot)
			call matrice_init(Gamma,Ntot)
	
			call matrice_init(KK,Ntot*Ntot)
			call matrice_init(RR,Ntot*Ntot)
			call matrice_init(RRD,Ntot*Ntot)
		
			call matrice_init(Rfeed,Ntot*Ntot)
			call matrice_init(LL,Ntot*Ntot)
			call matrice_init(LL_0,Ntot*Ntot)
			call matrice_init(LD,Ntot*Ntot)
			call matrice_init(UU,Ntot*Ntot)
			call matrice_init(U1,Ntot*Ntot)
			call matrice_init(TT,Ntot*Ntot)
!			call matrice_init(ZZ,Ntot*Ntot)
!			call matrice_init(Z1,Ntot*Ntot)
			call matrice_init(Udt,Ntot*Ntot)
		
			call matrice_init(rho_m_0,Ntot)
!			call matrice_init(rho_m,Ntot)

			allocate(Amp(Ntot, Ntot, agg%nu))
			Amp = 0.

			allocate(sigdE(Tsteps+1,(agg%nu*Ntot*(Ntot-1)/2)))
			sigdE= 0.
			
			allocate(dummy(Ntot,agg%nu,Tsteps+1))
			dummy=0.

			allocate(gft2(Ntot))
			gft2=0.0
			
		if (iterate==1) then	
	
			OD = 0.0
	 
	 		call matrice_init(SSavg,Ntot)
			call matrice_init(Enavg,Ntot)
			call matrice_init(DMavg,Ntot)

			allocate(rho_m_t(Ntot,Ntot,Tsteps+1))
!			allocate(rho_v(Ntot*Ntot))
			allocate(pop(Tsteps+1,Ntot*agg%nu))
			allocate(coh(Tsteps+1,Ntot*(Ntot-1)/2))
			allocate(tmem(Tsteps+1))
			allocate(std(Tsteps+1))
			rho_m_t = 0.0
			pop = 0.0
			coh = 0.0
			tmem = 0.0
!			rho_v	= 0.0

			allocate(maxAmp(agg%nu*Ntot*(Ntot-1)/2))
			maxAmp = 0.

			call matrice_init(Xsi,Ntot)

			allocate(Amp_t(Ntot,Ntot,agg%nu,Tsteps+1))
			allocate(sig(Tsteps+1,(agg%nu*Ntot*(Ntot-1)/2)))
			allocate(totsig(Tsteps+1))
			Amp_t  = 0.
			sig    = 0.
			totsig = 0.
			
!			allocate(gft2(Ntot))
!			gft2=0.0

		end if

!---

		
!---    Prepare Hamiltonian and solve eigenstate problem
		call Hamiltonian(HH, agg, wHin)
		
		! Solve eigenvalues and eigenstates problem of the Hamiltonian
		call spec(HH,SS,En)
		call inv(SS, S1)


!---    Set initial populations in the excitonic basis
		select case (pathway)
		case ('R1','R2') ! coherences in the excited-states. Will dephase. 
			call init_pops(rho_m_0, Amp, Nsrc, agg, SS, S1, DMavg)
		case ('R4','R3')	! coherences in the ground state
							! initial everything at 1 for the dynamics. Real amplitudes are taken cared of after call to dynamics
							! Should be ok because R=1 : no transfer between different coherences.
			rho_m_0=1.
			Amp=1.
		end select

!		write(6,*) '*******SS*******'
!		write(6,formNtot) real(transpose(SS))

		!-- Diverse Outputs
		if (iterate<2) then
			write(6,*) 'Ntot = ', Ntot
			write(6,'(a8, f6.1, a9, f8.3, a3)') 'wH = ', wHin, ' cm-1 = ', wHin*cm2eV, ' eV'
			write(6,*) '*******E*******'
			write(6,formN) agg%E
			write(6,*) '*******J*******'
			write(6,formN) agg%J
			write(6,*) '*******HH*******'
			write(6,formNtot) real(transpose(HH))
			do aa=1, Ntot
				write(6,*) aa, real(En(aa,aa))
			end do

			write(6,*) '*******En*******'
			write(6,formNtot) real(transpose(En))
!			write(6,*) '*******S1 * HH* SS*******'
!			write(6,formNtot) real(transpose(matmul(S1,matmul(HH,SS))))
			write(6,*) '*******SS*******'
			write(6,formNtot) real(transpose(SS))
!			write(6,*) '*******S1*******'
!			write(6,formNtot) real(transpose(S1))
!			write(6,*) '*******S1 * SS*******'
!			write(6,formNtot) real(transpose(matmul(S1,SS)))
!			write(6,*) '*******SS * HH*******'
!			write(6,formNtot) real(transpose(matmul(SS,HH)))
!			write(6,*) '*******HH* SS*******'
!			write(6,formNtot) real(transpose(matmul(HH,SS)))
!			write(6,*) '*******SS* En*******'
!			write(6,formNtot) real(transpose(matmul(SS,En)))
!			write(6,*) '*******En* SS*******'
!			write(6,formNtot) real(transpose(matmul(En,SS)))
			write(6,*) '********* real(rho_m_0) in exc. basis **********'
			write(6,formNtot) transpose(real(rho_m_0))		
			write(6,*) '********* real(Amp_0) in exc. basis **********'
			write(6,formNtot) transpose(real(Amp(:,:,1)))		
			
!			write(6,*) '********* real(rho_m_0) in loc. basis **********'
!			write(6,formNtot) transpose(real(matmul(SS,matmul(rho_m_0,S1))))

!---  		Output Spectral density and compare different models
			open(unit=100, file= trim(file_join(out_dir,"SpD.dat")))
			write(100,'(8a15)') 'w', 'SpD(w)', 'CrD(w)', 'Crl(w)', 'nthermal(w)', 'Crl(-w)', 'SpD(-w)', 'CrD(-w)'
			do aa = 0, 600 	! [cm-1]
				w = aa*1._dp * cm2eV
			    write(100,'(8e15.3)') real(w), real(SpD(w)), real(CrD(w)), real(Crl(w)), real(nthermal(w)), real(Crl(-1.*w)), real(SpD(-1.*w)), real(CrD(-1.*w)) !CrD should be compared with Crl
			end do 
			close(unit=100)
!---		
	
		end if
!---
		
		write(6,*) '********* dE = ', dgapi

!---

		
		
!---	Define correlation function in the site basis (sites are assumed to be totally correlated in the first instance)
		cc = 0.
		if (oneParticle) then
			do nn=2,(agg%N-1)
				do mm=2,(agg%N-1)
					do aa=1,agg%nu
						do bb=1,agg%nu
							cc(agg%nu*(nn-1)+aa,agg%nu*(mm-1)+bb)= delta(nn,mm)!*delta(aa,bb)
						end do
					end do
				end do
			end do
		else
!			do nn=2,(agg%N-1)
!				do mm=2,(agg%N-1)
!					do aa=1,agg%nu
!						do bb=1,agg%nu
!							do aam=1,agg%nu
!								do bbm=1,agg%nu
							forall (nn=2:(agg%N-1),mm=2:(agg%N-1),aa=1:agg%nu,bb=1:agg%nu,aam=1:agg%nu,bbm=1:agg%nu) &
							& cc((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(aa-1)+bb,(agg%nu**(agg%N-2))*(mm-1)+agg%nu*(aam-1)+bbm)= delta(nn,mm)
!								end do
!							end do
!						end do
!					end do
!				end do
!			end do

		end if 
		write(6,*) '*******cc*loc***'
		write(6,formNtot) transpose(real(cc))
!---

!---	Prepare the coefficients to compute the correlation function in the exciton basis
!       Later used to compute the pure dephasing rates
!		if (pureDephasing=='yes') then
			call Gamma_init(Gamma, SS, agg)
!		end if
!---		


!---	Compute the relaxation tensor
		if (doRelax == 'yes') then
			call Relaxation(RR, En, SS, S1, cc, RRD)
		end if
!---		
		
		
		
!---    Prepare transformation into the vibronic basis
		! Calculate the matrix with projection of the local basis onto the eigenvectors		
!		forall (nn=1:Ntot, mm=1:Ntot, aa=1:Ntot, bb=1:Ntot) KK(Ntot*(nn-1)+mm, Ntot*(aa-1)+bb)=SS(mm,bb)*S1(aa,nn)

!		forall(aa=1:Ntot, bb=1:Ntot, nn=1:Ntot, mm=1:Ntot)
!			ZZ(Ntot*(bb-1)+aa, Ntot*(mm-1)+nn) = SS(aa, nn) * S1(mm, bb)	!ZZ(an,mb) = SS(an)S1(mb) = S1(mb)SS(an)
!		end forall

!		call inv(ZZ, Z1)
!---		
	
		!call set_SuperOp and define Rfeed
		
		
!		RR = RR + Rfeed
!		write(116,*) '****vib: RRD R****'
!		write(116,formNtotNtot) transpose(real(RRD))
!		write(116,*) '****vib: RRD I****'
!		write(116,formNtotNtot) transpose(imag(RRD))
		
!---    Define the evolution operator

		call reversible_dynamics(LL_0,En)	

!		write(116,*) '****vib: LL_0 R****'
!		write(116,formNtotNtot) transpose(real(LL_0))
!		write(116,*) '****vib: LL_0 I****'
!		write(116,formNtotNtot) transpose(imag(LL_0))
		
		select case (pathway)
			case ('R1','R2')	! pathways through electronic excited state, dephasing of excited states
				LL = LL_0 - (0., 1.)*RRD
			case ('R4','R3')	! pathways through electronic ground state, no dephasing of vib. coherence
				LL = LL_0 
		end select

		call spec(LL,UU,LD)
		call inv(UU, U1)

!	   write(116,*) '****vib: LL R****'
!	   write(116,formNtotNtot) transpose(real(LL))
!	   write(116,*) '****vib: UU R****'
!	   write(116,formNtotNtot) transpose(real(UU))
!	   write(116,*) '****vib: UU I****'
!	   write(116,formNtotNtot) transpose(imag(UU))
!	   write(116,*) '****vib: LD R****'
!	   write(116,formNtotNtot) transpose(real(LD))
!	   write(116,*) '****vib: LD I****'
!	   write(116,formNtotNtot) transpose(imag(LD))

!       Testing error in the diagonalization (test only on the real values)
		TT = (matmul(U1,matmul(LL,UU))-LD)
		maxi = (0.,0.)
		do aa=1,Ntot*Ntot
			do bb=1, Ntot*Ntot
				if (Abs(Real(TT(aa, bb))) > Real(maxi))  maxi = Abs(TT(aa, bb))
!				if (Abs(Imag(TT(aa, bb))) > Imag(maxi))  maxi = Abs(Imag(TT(aa, bb)))
			end do
		end do
		if (abs(maxi)>0.01) write(6,*) '*********WARNING*****  \n	vib Maximum error in LL diagonalization = ', abs(maxi)

!		Construct the evolution superoperator for step dt in the diagonal basis
		Udt = (0., 0.)
		forall(aa=1:Ntot*Ntot) Udt(aa, aa)=exp(LD(aa,aa) * (0.,-1.)*dt)

!		write(116,*) '****vib: Udt diaR****'
!		write(116,formNtotNtot) transpose(real(Udt))
!		write(116,*) '****vib: Udt diaI****'
!		write(116,formNtotNtot) transpose(imag(Udt))

!		Transforming back in the vibronic basis
		Udt = matmul(UU, matmul(Udt, U1))

!		write(116,*) '****vib: Udt****'
!		write(116,formNtotNtot) transpose(real(Udt))
!---



!---    Dynamical loop
!		write(6,*) 'Start dynamical loop'
		call dynamics(rho_m_0, Amp, Udt, Gamma, rho_m_t, Amp_t, Tsteps, loc, SS, S1, agg, tmem)
!---

		if (pathway=='R4' .or. pathway=='R3') then	! multiply the evolved signal by the amplitude of the dipole moments for the ground-state pathways (account for the 4 interactions with laser pulses)
			call init_pops(rho_m_0, Amp, Nsrc, agg, SS, S1, DMavg)
			dummy = Amp_t(1,:,:,:)			
			do aa=1,Ntot	! Evolved signal is for state g_k g_0, i.e. in Amp_t(1,1,xi,t)
				do bb=1, Ntot
					do xi=1,agg%nu
						Amp_t(aa,bb,xi,:) = dummy(xi,xi,:) * Amp(aa,bb,xi)
!						if (aa==bb) write(6,*) aa, xi,dummy(xi,xi,1000),Amp_t(1,2,xi,1000), Amp(aa,bb,xi)
					end do
				end do 
			end do
		end if

!		write(6,*) '*********Amp_t_0*********' 
!		write(6,'(8f10.4)') transpose(real(Amp_t(:,:,1,1000)))
!		write(6,*) '*********Amp_t_1*********'
!		write(6,'(8f10.4)') transpose(real(Amp_t(:,:,2,1000)))

!--- 	Compute the signal for w1, w3 (after FT of the response function over t1 and t3)

! 		non-rephasing pathway: the same state (aa) in involved in the optical coherence at t1 and t3
!		WARNING: in the two-particule model, we should decay to a different ground state, i.e. to n=g, nu.ne.0, than during the excitation
!		aa  = agg%nu+1				! lowest exciton (pos_g0+1)
!		opt = (agg%N-2)*agg%nu+1	!WARNING: to generalize for FMO or N levels 
!		opt = 1						!WARNING: to generalize for FMO or N levels 
!		One particule model
!		aa  = agg%nu+1				! lowest exciton (pos_g0+1)
!		opt = 1						!WARNING: to generalize for FMO or N levels 
!		w   = w1*cm2eV - (En(aa,aa)-En(opt,opt))
!		call signal2D(w, -w, aa, aa, opt, opt, rho_m_t, RRD, Gamma)

!		aa  = agg%nu+1				! lowest exciton (pos_g0+1)
!		opt = 1						!WARNING: to generalize for FMO or N levels 

!---    compute linear spectrum 
		if (pathway == 'OD') then
			open(unit=100, file= trim(file_join(out_dir,"OD.dat")))
			call init_pops(rho_m_0, Amp, Nsrc, agg, SS, S1, DMavg)
			do nn=1,2000	! frequency in [cm-1]
				w = (11500. + nn)
				do aa=agg%nu+1,Ntot-agg%nu		! eigenvalues in the excited manifold
					dw  = (En(aa,aa)-En(1,1)) - w*cm2eV
					OD(nn) = OD(nn) + Amp(aa,1,1)*GFT(dw,aa,1,RR,Gamma)
!					write (6,*) aa, OD(1)
				end do		
!				OD(nn) = OD(nn) + w * cm2eV * OD(nn)*PDF
				OD(nn) = w * cm2eV * OD(nn)
				if (iterate == Nruns) then
					write (6,*) 'Print OD with PDF = ', PDF
					write(100,'(f15.5,e15.3)') real(w), real(OD(nn))
				end if
			end do
			! averaging over energy disorder
			close(unit=100)
		end if
!---

		
		call signal2D(w1, w3, Amp_t, En, RRD, Gamma, agg, gft2)

!	   write(6,*) '*********Amp_t_nu0 at t=1000*********' 
!	   write(6,formNtot) transpose(real(Amp_t(:,:,1,1000)))
!	   write(6,*) '*********Amp_t_nu1 at t=1000*********'
!	   write(6,formNtot) transpose(real(Amp_t(:,:,2,1000)))
		
		
!		output Fourier transform of the goft		
		if (iterate==1) then		
			open(unit=102, file= trim(file_join(out_dir,"GFT.dat")))
			write(102,'(7a15)') 'Re(dw)', 'Re(GFT(w))', 'Im(GFT(w))', 'Re(GFT(w)GFT(w))'
			do nn = -200, 200 	! [cm-1]
				w = nn*(1,0.) * cm2eV
				opt = 1
				aa  = agg%nu+1
				write(102,'(4e15.3)') real(w), real(GFT(w, aa, opt, RRD, Gamma)), imag(GFT(w, aa, opt, RRD, Gamma)), real(GFT(w, aa, opt, RRD, Gamma)*GFT(w, aa, opt, RRD, Gamma))
			end do 
			close (unit=102)
		end if
		
!---


!---	Energetic characterisic averaged over the energetic disorder
!		Composition of the eigenvectors
		do aa=1,Ntot
			do bb=1,Ntot
				SSavg(aa,bb) = SSavg(aa,bb) + PDF*(SS(aa,bb))**2
			end do
		end do 
!		Averaged eigenvalues		
		Enavg = Enavg + PDF*En
!---


!---	Evaluate the vibrational character of each coherences
		do aa=1, Ntot		!loop over the eigenstates
			do bb=1, Ntot	
				nn = 2		! evaluates the vibrational coherence on excited state nn. (nn=2 for exciton 1)
				Xsi(aa,bb) = SS(agg%nu*(nn-1)+1, aa)**2 * SS(agg%nu*(nn-1)+2, bb)**2 + SS(agg%nu*(nn-1)+2, aa)**2 * SS(agg%nu*(nn-1)+1, bb)**2 
				Xsi(aa,bb) = Xsi(aa,bb) + Xsi(aa,bb) * PDF	! Average over energy disorder
			end do
		end do
!---


!---	Output populations and coherences

!		Populations (diagonal terms)
!		do aa = agg%nu**(agg%N-2)+1,(agg%N-1)*agg%nu**(agg%N-2)
		do aa = 1, Ntot
			do xi=1, agg%nu
				pop(:,(xi-1)*Ntot + aa) = pop(:,(xi-1)*Ntot + aa) + PDF * Amp_t(aa,aa,xi,:)
				if (pathway=='R4'.or.pathway=='R3') then
					sigdE(:,(xi-1)*Ntot + aa) = Amp_t(aa,aa,xi,:)
					totsig(:)    = totsig(:)    + PDF * Amp_t(aa,aa,xi,:)
				end if
			end do
		end do
		
!		Total of populations
		do tdt=1,(Tsteps+1)
			std(tdt) = trace(rho_m_t(:,:,tdt))
		end do
		
!		Coherence
		nn = 0
		do aa = 1,Ntot
			do bb = 1,Ntot
				if (aa>bb) then
					nn = nn + 1
					if ( (aa<Ntot) .and. (bb>=1) ) then
						coh(:,nn) = coh(:,nn) + PDF * rho_m_t(aa,bb,:)
						do xi=1, agg%nu
							! Signal for each pathway
							sig(:,(xi-1)*(Ntot*(Ntot-1)/2)+nn) = sig(:,(xi-1)*(Ntot*(Ntot-1)/2)+nn) + PDF * Amp_t(aa,bb,xi,:)
							! Signal for each pathway, not averaged over disorder
							if (pathway=='R1' .or. pathway=='R2') then
								sigdE(:,(xi-1)*(Ntot*(Ntot-1)/2)+nn) = Amp_t(aa,bb,xi,:)
								! Sum over all pathways
								totsig(:)    = totsig(:)    + PDF * Amp_t(aa,bb,xi,:)
							end if
						end do
					end if
				end if
			end do
		end do
		
		
		
!       Total signal of the calculated pathway: sum the contribution of all coherences aa, bb for all vibrational level xi.
!		do aa = agg%nu**(agg%N-2)+1,(agg%N-1)*agg%nu**(agg%N-2)
!			do bb = agg%nu**(agg%N-2)+1,(agg%N-1)*agg%nu**(agg%N-2)
!				do xi=1, agg%nu
!					totsig(:)    = totsig(:)    + PDF * Amp_t(aa,bb,xi,:)
!				end do
!			end do
!		end do

!       Save the total signal in spec_tot
		spec_tot(ww,ww3,:) = totsig(:)
		

!---


!---	Compute and output relaxation time
!		call relax_time(RRD, Gamma, Tsteps, wHin)
!---


!---	Save the maximum amplitude at different times 
		if (iterate==1) write(22,'(/a6, i4)') '# i = ',ww-1
		maxAmp = 0.
		if (pathway=='R1') then
			do tdt= 1, 300	! Max at initial time, coherence 12
				if (Abs(Real(sigdE(tdt     ,6))) > maxAmp(1))  maxAmp(1) = Abs(Real(sigdE(tdt     ,6)))
				if (Abs(Real(sigdE(tdt+500 ,6))) > maxAmp(2))  maxAmp(2) = Abs(Real(sigdE(tdt+500 ,6)))
				if (Abs(Real(sigdE(tdt+1000,6))) > maxAmp(3))  maxAmp(3) = Abs(Real(sigdE(tdt+1000,6)))
				if (Abs(Real(sigdE(tdt+1500,6))) > maxAmp(4))  maxAmp(4) = Abs(Real(sigdE(tdt+1500,6)))
				if (Abs(Real(sigdE(tdt+2000,6))) > maxAmp(5))  maxAmp(5) = Abs(Real(sigdE(tdt+2000,6)))
			end do
			write(22,'(2f15.2,9e15.3)') Real(w1), dgapi, Real(maxAmp(1:5)), Real(gft2(3)), Imag(gft2(3))
		else if (pathway=='R4') then
			do tdt= 1,1000	
				if (Abs(Real(sigdE(tdt     ,12))) > maxAmp(1))  maxAmp(1) = Abs(Real(sigdE(tdt     ,12)))   ! dimer
				if (Abs(Real(sigdE(tdt     ,10))) > maxAmp(2))  maxAmp(2) = Abs(Real(sigdE(tdt     ,10)))    ! monomer
			end do
			write(22,'(2f15.2,7e15.3)') Real(w1), dgapi, Real(maxAmp(1:2)), Real(gft2(4)), Imag(gft2(4))
		end if
!---


!---    Clean memory			
		call aggregate_deinit(agg)	

		call matrice_deinit(HH)
		call matrice_deinit(SS)
		call matrice_deinit(S1)
		call matrice_deinit(En)

		call matrice_deinit(cc)
		call matrice_deinit(Gamma)

		call matrice_deinit(KK)
		call matrice_deinit(RR)
		call matrice_deinit(RRD)
		
		call matrice_deinit(Rfeed)
		call matrice_deinit(LL)
		call matrice_deinit(LL_0)
		call matrice_deinit(LD)
		call matrice_deinit(UU)
		call matrice_deinit(U1)
		call matrice_deinit(TT)
!		call matrice_deinit(ZZ)
!		call matrice_deinit(Z1)
		call matrice_deinit(Udt)
		
		call matrice_deinit(rho_m_0)
		
		deallocate(Amp)
		deallocate(sigdE)
		deallocate(dummy)

		deallocate(gft2)
		
!=== end the loop over energy disorder
	end do
		
!=== Output results of the dynamics, averaged of the runs (energetic disorder of the gap)		

!---	Prepare output format
!		Populations (real part) for all nu
		write(dim,'(i4)') Nnu*Ntot
		formPop = '(f15.2,'//trim(adjustl(dim))//'e15.3)'
!		Coherences: N (N-1)/2 , real, abs
		write(dim,'(i4)') 3*Ntot*(Ntot-1)/2
		formCoh = '(f15.2,'//trim(adjustl(dim))//'e15.3)'
		write(dim,'(i4)') 2*Nnu*Ntot*(Ntot-1)/2
		formSig = '(f15.2,'//trim(adjustl(dim))//'e15.3)'
!		write(6,*) 'formPop, formCoh, formSig: ', formPop, formCoh, formSig
		

!---	Write into files
!		write(out_file,'(f7.0,a3)') real(w1),"dat"
		write(out_file0,'(f7.0)') real(w1)
		write(out_file,'(f7.0,a6,a3)')  real(wHin),trim(adjustl(out_file0)),"dat"
		open(unit=111,file=trim(file_join(out_dir,"sigdi"//trim(adjustl(out_file)))))
!		write(111, '(13a16)') 'Time_step','Real','Real','Real','Real','Real','Imag','Imag','Imag','Imag','Imag','Re(std)','Im(std)'
		do tdt = 1,Tsteps
!			write (111,11) tmem(tdt), real(pop(tdt,:)), imag(pop(tdt,:)), real(std(tdt)), imag(std(tdt))
			write (111,formPop) tmem(tdt), real(pop(tdt,:))
		end do
		close(unit=111)
		
		
!		write(out_file,'(f7.0,a3)') real(w1),"dat"
!		open(unit=222,file=trim(file_join(out_dir,"coh"//trim(adjustl(out_file)))))
!		do tdt = 1,Tsteps
!			write (222,formCoh) tmem(tdt), real(coh(tdt,:)), imag(coh(tdt,:)), abs(coh(tdt,:))
!		end do
!		close(unit=222)

		
!		write(out_file,'(f7.0,a3)') real(w1),"dat"
!		open(unit=333,file=trim(file_join(out_dir,"sig"//trim(adjustl(out_file)))))
		write(out_file0,'(f7.0)') real(w1)
		write(out_file,'(f7.0,a6,a3)')  real(wHin),trim(adjustl(out_file0)),"dat"
		open(unit=333,file=trim(file_join(out_dir,"sig"//trim(adjustl(out_file)))))
		do tdt = 1,Tsteps
			write (333,formSig) tmem(tdt), real(sig(tdt,:)), imag(sig(tdt,:))
		end do
		close(unit=333)

		write(out_file,'(f7.0,a3)') real(w1),"dat"
		open(unit=444,file=trim(file_join(out_dir,"totsig"//trim(adjustl(out_file)))))
		do tdt = 1,Tsteps
			write (444,'(f15.2,3e15.3)') tmem(tdt), real(totsig(tdt)), imag(totsig(tdt)), abs(totsig(tdt))
		end do
		close(unit=444)
		

!		write (555,'(2f15.2,100e25.4)') real(w1), real(w3), real(gft2(:)), imag(gft2(:))

		
!
		if (spectrum) then
			if (ww==1) write(11,'(/a6, i4)') '# i = ',ww3-1
			maxSig = 0.
			do tdt=  1000, (Tsteps+1)	! Max amplitude after 1ps
				if (Abs(Real(totsig(tdt))) > maxSig)  maxSig = Abs(Real(totsig(tdt)))
			end do
			write(11,'(2f15.2,22e25.4)') real(w1), real(w3), maxSig, &
			& Real(totsig(1))   , Real(totsig(100)) , Real(totsig(200)) , Real(totsig(300)) , Real(totsig(400)) , Real(totsig(500)) ,	&	
			& Real(totsig(600)) , Real(totsig(700)) , Real(totsig(800)) , Real(totsig(900)) , Real(totsig(1000)), Real(totsig(1100)),	&	
			& Real(totsig(1200)), Real(totsig(1300)), Real(totsig(1400)), Real(totsig(1500)), Real(totsig(1600)), Real(totsig(1700)),	&	
			& Real(totsig(1800)), Real(totsig(1900)), Real(totsig(2000))
		end if 

!		Save the maximum amplitude after decay of the electronic coherence (for time>700s)
		maxAmp = 0.
		do aa=1, (Nnu*Ntot*(Ntot-1)/2)
			do tdt= 1000, (Tsteps+1)	! 700 is arbritrarily chosen, but should be big enough to be beyond the decay of any electronic coherence
				if (Abs(Real(sig(tdt,aa))) > maxAmp(aa))  maxAmp(aa) = Abs(Real(sig(tdt,aa)))
			end do
		end do
		write(10,formCoh) Real(w1), real(wHin), Real(maxAmp), real(Xsi(3,4))
		
!		Save the vibrational character of the coherence
		open(unit=222,file=trim(file_join(out_dir,"Xsi.dat")))
		write(222,'(2a6,2a20)') '#   aa', 'bb', 'Real(Xsi(aa,bb))',  'Imag(Xsi(aa,bb))'
		write(222,formNtot) real(transpose(Xsi))
		write(222,'(//a15)') '#=== SSavg === (<SS^2>)'
		write(222,formNtot) real(transpose(SSavg))
		write(222,'(//a15)') '#=== Enavg ==='
		write(222,formNtot) real(transpose(Enavg))
		write(222,'(//a15)') '#=== DMavg ==='
		write(222,formNtot) real(transpose(DMavg))
		close(unit=222)
		
! Save eigenenergies and vibrational character of the coherences of interest
		write(33,'(2f15.0,6f15.5,10e15.3)') real(wHin), real(w1), real(Enavg(1,1)), real(Enavg(2,2)), real(Enavg(3,3)), &
     &             real(Enavg(4,4)), real(Enavg(5,5)), real(Enavg(6,6)), real(Xsi(3,4)), real(Xsi(3,5)), &
     &             real(SSavg(4,4)),real(SSavg(4,5)),real(SSavg(5,4)),real(SSavg(5,5)), &
     &             real(SSavg(3,3)),real(SSavg(3,4)),real(SSavg(4,3)),real(SSavg(5,3))   
		
!		
		
!===
		
		
!---    Clean memory	

		call matrice_deinit(SSavg)
		call matrice_deinit(Enavg)
		call matrice_deinit(DMavg)
		
		deallocate(rho_m_t)
!		call matrice_deinit(rho_m)
!		deallocate(rho_v)
		deallocate(pop)
		deallocate(coh)
		deallocate(tmem)
		deallocate(std)
		
		deallocate(maxAmp)
		call matrice_deinit(Xsi)

		deallocate(Amp_t)
		deallocate(sig)
		deallocate(totsig)
		
!		deallocate(gft2)
!---

!###    End loop over wH
   end do 

!###    End loop of w1, ww
   end do 
!###    End loop of w3, ww3
! SPECTRUM
  end do 

!--- Output the spectrum in matrix form  
		open(unit=11, file=trim(file_join(out_dir,"spec_tot.dat")))
		write(dim, '(i4)') spec_fq(3)
		formSpec='('//trim(adjustl(dim))//"e13.3"//')'
		do nn=1, 10	! number of output in time
			write(11, formSpec) real(transpose(spec_tot(:,:,1+(100*nn-1))))
			write(11,'(//)') 
			write(11,'(//)') 
		end do
!---

		deallocate(spec_tot)
   
   close(unit=10)
   close(unit=555)
   close(unit=11)
   close(unit=33)
		

	end subroutine main_vibronic_networks


	subroutine matrice_init(MM,N)
	
		integer, intent(in)		:: N
		complex(dpc), dimension(:,:),allocatable		:: MM
		
		if (.not. allocated(MM)) then
			allocate(MM(N,N))
			MM = (0.,0.)
		end if
		
	end subroutine matrice_init
	
	subroutine matrice_deinit(H)
		complex(dpc), dimension(:,:),allocatable	:: H
		
		if (allocated(H)) then
			deallocate(H)
		end if
		
	end subroutine matrice_deinit
	
	subroutine aggregate_init(agg, wHin)
	
		type(aggregate), intent(out)	:: agg
		real(dp), intent(in)			:: wHin
		real(dp)						:: dd	! dimensionless coordinate shift
	
		integer	:: nn ! used for loops
		
		! Allocation of memory for all arguments of the aggregate
		allocate(agg%N)
		allocate(agg%nu)
		allocate(agg%theta)
		allocate(agg%gap)
		allocate(agg%dwidth)

		select case (molSystem)
			case ('FMO')
				agg%N  = 9
				agg%nu = 2	! >=1 (zero-ground state)
			case ('dimer')
				agg%N  = 4	!(+2 states for the optical transitions)
				agg%nu = 2	! >=1 (zero-ground state)
				! WARNING: in future change the dimension of theta for FMO
				agg%theta=1.870193147	! angle between the dipole moment of BChl 3 and 4
			case ('jonas')	! parameters used in Tiwari2012a paper
				agg%N  = 4	!(+2 states for the optical transitions)
				agg%nu = 2	! >=1 (zero-ground state)
				agg%theta=1.870193147	! angle between the dipole moment of BChl 3 and 4
			case ('monomer')
				agg%N  = 3	!(+2 states for the optical transitions)
				agg%nu = 2	! >=1 (zero-ground state)		
		end select
			
		allocate(agg%d(agg%N))
		allocate(agg%E(agg%N))
		allocate(agg%J(agg%N,agg%N))            
		allocate(agg%ori(agg%N*agg%nu,3)) 	!(x,y,z)  
		agg%d  =0.
		agg%E  =0.		  
		agg%J  =0.
		agg%ori=0.
		
!		Define the orientation of the molecules
		do nn=1,(agg%N*agg%nu)
			agg%ori(nn,1) = (1.,0.)	! oriente all molecules along x
		end do
!		rotate the second site by theta !WARNING: to be completed to account for the full FMO. Only valid for dimer		
		agg%ori(2*agg%nu+1,1)  = cos(agg%theta)
		agg%ori(2*agg%nu+2,1)  = cos(agg%theta)
		agg%ori(2*agg%nu+1,2)  = sin(agg%theta)
		agg%ori(2*agg%nu+2,2)  = sin(agg%theta)


  		dd = sqrt(2*HR)	! Huang-Rys factor S=0.05 [Christensson_2012]
!		dd = 0.
!  		dd = sqrt(2*0.05)	! Huang-Rys factor S=0.05 [Christensson_2012]

!		Define the site energy and respective coupling
		select case (molSystem)

    	   case ('FMO')

				agg%d = (/ dd, dd, dd, dd, dd, dd, dd, dd, dd /)		   ! Shift of every oscillator relative to the ground state

!				FMO with 7 BChls: Energies E and coupling J from Adolphs_BJ_2006_2778
!				Energie shift of E3 = 12 210 cm -1
!				Table 4, Trimer C. Tepidum
				agg%E = (/ -12210., 200., 320., 0.0, 110.0, 270.0, 420.0, 230.0, 3500./)*cm2eV
!				Table 4, Monomer C. Tepidum
!				Energie shift of E3 = 12 205 cm -1
!				agg%E = (/ -290., 240., 315., 0.0, 130.0, 285.0, 435.0, 245.0, 290./)*cm2eV
!				Table 5, C.tepidum
!				agg%E = (/ -290., 220.,  55., 0.0, 225.0,  35.0, 300.0, 105.0, 290./)*cm2eV
				agg%J = transpose(reshape((/ 						 &
				&	  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,&
				&	  0.0,  0.0,-87.7,  5.5,  -5.9,  6.7,-13.7, -9.9, 0.0,&
				&	  0.0,-87.7,  0.0, 30.8,   8.2,  0.7, 11.8,  4.3, 0.0,&
				&	  0.0,  5.5, 30.8,  0.0, -53.5, -2.2, -9.6,  6.0, 0.0,&
				&	  0.0, -5.9,  8.2,-53.5,   0.0,-70.7,-17.0,-63.3, 0.0,&
				&	  0.0,  6.7,  0.7, -2.2, -70.7,  0.0, 81.1, -1.3, 0.0,&
				&	  0.0,-13.7, 11.8, -9.6, -17.0, 81.1,  0.0, 39.7, 0.0,&
				&	  0.0, -9.9,  4.3,  6.0, -63.3, -1.3, 39.7,  0.0, 0.0,&
				&	  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0/)*cm2eV, &
				&						(/agg%N,agg%N/)))			   !Size of the allocated matrix

			case ('dimer')

				agg%d = (/ dd, dd,  dd, dd/)	! Shift of every oscillator relative to the ground state
!				agg%gap    = wHin			! initial gap energy
!				agg%gap    = 181.3_dp			! initial gap energy
 				agg%gap    = 110._dp			! initial gap energy
				agg%dwidth = 80._dp/2.*sqrt(log(2.))	!width_gap = sqrt(2)*width_individual_pigment
				
				if (Nruns==1) dgapi = agg%gap 	! start with the original gap energy for the one run
				
				agg%E = (/ -12210._dp, 0._dp, dgapi, 3500._dp /)*cm2eV		! Site Energy	
				agg%J = transpose(reshape((/ 0.,   0. ,   0. , 0., & 
											 0.,   0. , -53.5, 0.,  &
											 0., -53.5,   0. , 0.,  &
											 0.,   0.,    0. , 0.  /)*cm2eV,  &
											 (/agg%N,agg%N/)))				!Size of the matrix to allocate
				
				
			case ('jonas')

				agg%d = (/ dd, dd,  dd, dd/)	! Shift of every oscillator relative to the ground state
 				agg%gap    = 150._dp			! initial gap energy
				agg%dwidth = 80._dp/2.*sqrt(log(2.))	!width_gap = sqrt(2)*width_individual_pigment
				
				if (Nruns==1) dgapi = agg%gap 	! start with the original gap energy for the one run
				
				agg%E = (/ -12210._dp, 0._dp, dgapi, 3500._dp /)*cm2eV		! Site Energy	
				agg%J = transpose(reshape((/ 0.,   0. ,   0. , 0., & 
											 0.,   0. ,  66.0, 0.,  &
											 0.,  66.0,   0. , 0.,  &
											 0.,   0.,    0. , 0.  /)*cm2eV,  &
											 (/agg%N,agg%N/)))				!Size of the matrix to allocate
				


			case ('monomer')

				agg%d = (/ dd, dd, dd/)	! Shift of every oscillator relative to the ground state
!				agg%gap    = wHin			! initial gap energy
!				agg%gap    = 110._dp			! initial gap energy
 				agg%gap    = 0._dp			! initial site energy
				agg%dwidth = 80._dp/2.*sqrt(2.*log(2.))
				
				if (Nruns==1) dgapi = agg%gap 	! start with the original gap energy for the one run
				
				agg%E = (/ -12210._dp, dgapi, 3500._dp /)*cm2eV		! Site Energy	
				agg%J = transpose(reshape((/ 0.,   0. ,  0., & 
											 0.,   0. ,  0.,  &
											 0.,   0.,   0.  /)*cm2eV,  &
											 (/agg%N,agg%N/)))				!Size of the matrix to allocate
			
		end select

		
				
!		Flag to know if the energy gap is below the vibrational frequency
!		if (Abs(dgapi)>wHin) then
!			beyond_wH=1
!		else
!			beyond_wH=0
!		end if
!		if (dgapi<0.) then
!			neg_gap=1
!		else
!			neg_gap=0
!		end if

!		Calculate disorder except for the first run		
!		if (iterate>1) then
!			randomly calculate dgap with a gaussian distribution
!			call gap_disorder(agg)
!			agg%gap    = agg%gap + dgap	!gap energy with disorder

!			sum over a discretized normal distribution and control energy gap
!			agg%gap = dgapi	!control the energy gap directly in the loop
!		end if


	



	end subroutine aggregate_init


	subroutine aggregate_deinit(agg)
	
		type(aggregate), intent(out)		:: agg	
		
		! Deallocation of the aggregate type
		deallocate(agg%N)
		deallocate(agg%nu)
		deallocate(agg%theta)
		deallocate(agg%gap)
		deallocate(agg%dwidth)
!		deallocate(agg%w0)
		deallocate(agg%d)
		deallocate(agg%E)
		deallocate(agg%J)
		deallocate(agg%ori)

	
	end subroutine aggregate_deinit


	
	subroutine Hamiltonian(HH, agg, wHin)
		
		complex(dpc), dimension(:,:), intent(out)	:: HH		! Total Hamiltonian 
		type(aggregate), intent(in)					:: agg		! molecular aggregate
		real(dp), intent(in)						:: wHin		! frequency of the vibrational mode treated explicitely
		real(dp)									:: wH		! frequency of the vibrational mode treated explicitely
		
		integer			:: n, m 	! index of the excitonic (molecular) level
		integer			:: nuN, muN	! index of the vibrational level when molecule N is excited
		integer			:: nuM, muM	! index of the vibrational level when molecule M is excited
									! nuN, muM refers to the vib. level of the excited molecule
		
		real(dp)		:: dd,d		! for output only
		
		
		wH = wHin * cm2eV
		
		HH  = 0._dp
	
!		write(6,*) 'Prepare Hamiltonian, dim = ', size(HH,1)
	
		if (oneParticle) then ! H dimension (N*nu)*(N*nu)
			do n=1, agg%N
				do m=1, agg%N
					do nuN=1, agg%nu
						do muM=1, agg%nu
						
							if (n==m .and. nuN==muM) then 		! diagonal terms
							
								HH(agg%nu*(n-1)+nuN, agg%nu*(m-1)+muM) = agg%E(n) + (nuN-1)*wH
						
							else if (n.ne.m) then				! coupling terms
						
								HH(agg%nu*(n-1)+nuN, agg%nu*(m-1)+muM) = agg%J(n,m) * &	! WARNING: here with still excite from the nu=0 for the ground state
							      	& FC(1, nuN, agg%d(n)) * FC(1, muM, agg%d(m))
!							      	& franc_condon_factor((nu-1),wH,0,wH,agg%d(n)/sqrt(wH)) * franc_condon_factor((mu-1),wH,1,wH,agg%d(m)/sqrt(wH))
										
							end if
						
						end do					
!						write(6,'(a10,i4, 3e10.2)') 'nu, FC = ', nu, agg%d(n), franc_condon_factor((nu-1),wH,0,wH,agg%d(n)/sqrt(wH)), FC((nu-1), agg%d(n))
					end do
				end do
			end do
		else	! H dimension (N*nu^N)*(N*nu^N)
			! Molecule n
			do n=1, agg%N
				do nuN=1, agg%nu
					do muN=1, agg%nu
						! Molecule m
						do m=1, agg%N
							do nuM=1, agg%nu
								do muM=1, agg%nu
						
									if (n==m .and. nuN==muM .and. muN==nuM) then 		! diagonal terms
																						! nuN = nu(e)N, muM = mu(e)M
							
										HH((agg%nu**(agg%N-2))*(n-1)+agg%nu*(nuN-1)+muN, (agg%nu**(agg%N-2))*(m-1)+agg%nu*(muM-1)+nuM) = agg%E(n) + (nuN-1)*wH
						
									else if (n.ne.m) then				! coupling terms
						
										HH((agg%nu**(agg%N-2))*(n-1)+agg%nu*(nuN-1)+muN, (agg%nu**(agg%N-2))*(m-1)+agg%nu*(muM-1)+nuM) = &
										&   agg%J(n,m) * FC(nuM, nuN, agg%d(n)) * FC(muN, muM, agg%d(m)) ! FC(nu(g), nu(e))
!							      			& franc_condon_factor((nu-1),wH,0,wH,agg%d(n)/sqrt(wH)) * franc_condon_factor((mu-1),wH,1,wH,agg%d(m)/sqrt(wH))
										
									end if
								end do
							end do
						end do					
!						write(6,'(a10,i4, 3e10.2)') 'nu, FC = ', nu, agg%d(n), franc_condon_factor((nu-1),wH,0,wH,agg%d(n)/sqrt(wH)), FC((nu-1), agg%d(n))
					end do
				end do
			end do
		end if 
		


		if (iterate==1) then
			open(unit=101, file= trim(file_join(out_dir,"FC.dat")))
			write(101,'(11a20)') 'd', 'FC(0,0,d) wave-packet', '[Chang2005]', 'FC(0,1,d)', '[Chang2005]', 'FC(0,2,d)', '[Chang2005]', 'FC(1,0,d)', '[Chang2005]', 'FC(1,1,d)', '[Chang2005]' 
			do n = 0, 100	
			    dd = n* 0.1			! dd is dimensionless
				d =  dd / sqrt(wH)	! d dimension [m]
				write(101,'(11e20.3)') dd,  FC(1, 1, dd), franc_condon_factor(0,wH,0,wH,d),  &
					&	FC(1, 2, dd), franc_condon_factor(0,wH,1,wH,d),	&
					&	FC(1, 3, dd), franc_condon_factor(0,wH,2,wH,d), &	
					&	FC(2, 1, dd), franc_condon_factor(1,wH,0,wH,d),	&
					&	FC(2, 2, dd), franc_condon_factor(1,wH,1,wH,d)	
			end do 
			close(unit=101)
		end if
		
	end subroutine Hamiltonian

	
	
	subroutine Gamma_init(Gamma, SS, agg)
		
		integer										:: aa, bb, nn, mm		! indices used in the different loops
		integer										:: nuN, muN, nuM, muM	! indices used in the different loops
		type(aggregate), intent(in)					:: agg
		complex(dpc), dimension(:,:), intent(in)	:: SS
		complex(dpc), dimension(:,:), intent(out)	:: Gamma					! result matrix


!   	Prepare the coefficients to compute the correlation function in the exciton basis
		if (oneParticle) then
			do aa=1,Ntot
				do bb=1,Ntot
					do nn=1,agg%N
						do nuN=1,agg%nu
							do muN=1,agg%nu
									Gamma(aa,bb)= Gamma(aa,bb) + SS((agg%nu*(nn-1)+nuN),aa)**2 * SS((agg%nu*(nn-1)+muN),bb)**2
! 								Equivalent calculation
!								if (nuN==muN) then
!									Gamma(aa,bb)= Gamma(aa,bb) + SS((agg%nu*(nn-1)+nuN),aa)**2 * SS((agg%nu*(nn-1)+muN),bb)**2
!								else
!									Gamma(aa,bb)= Gamma(aa,bb) + 0.5*(SS((agg%nu*(nn-1)+nuN),aa)**2 * SS((agg%nu*(nn-1)+muN),bb)**2 + SS((agg%nu*(nn-1)+muN),aa)**2 * SS((agg%nu*(nn-1)+nuN),bb)**2)
!								end if 
							end do
						end do
					end do
				end do
			end do
		else 
			do aa=1,Ntot
				do bb=1,Ntot
					do nn=1,agg%N
						do nuN=1,agg%nu
							do muN=1,agg%nu
								do nuM=1,agg%nu
									do muM=1,agg%nu
										Gamma(aa,bb)= Gamma(aa,bb) +  &
									&                  SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nuN-1)+muN),aa)**2 &
									&                * SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(muM-1)+nuM),bb)**2
!										if (nuN==muM .and. muN==nuM) then ! same vibrational levels on the excited and ground state molecules
!											Gamma(aa,bb)= Gamma(aa,bb) + SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nuN-1)+muN),aa)**2 * SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(muM-1)+nuM),bb)**2
!										else
!											Gamma(aa,bb)= Gamma(aa,bb) +  &
!									&           0.5 *( SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nuN-1)+muN),aa)**2 * SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(muM-1)+muN),bb)**2 &
!									&              +   SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(muM-1)+muN),aa)**2 * SS(((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nuN-1)+muN),bb)**2 )
!										end if
									end do
								end do 
							end do
						end do
					end do
				end do
			end do

		endif
		
		
		! output gamma matrix
		if (iterate==1) then
			open(unit=111, file=trim(file_join(out_dir,"Gamma.dat")))
			if (oneParticle) then 
				write(111,'(8f5.2)') real(transpose(Gamma))
			else
				write(111,'(16f5.2)') real(transpose(Gamma))
			end if 
			close(111)
		end if
		
		
!       Prepare the matrix to compute the pure dephasing terms 
!		Gamma_ab = gam_aa + gam_bb - 2 gam_ab for coherent terms rho_ab
!		forall (nn=1:Ntot, mm=1:Ntot, aa=1:Ntot, bb=1:Ntot) &
!     &         Gamma(Ntot*(mm-1)+nn,Ntot*(bb-1)+aa) = (gam(nn,nn)+gam(mm,mm)-2*gam(nn,mm))*delta(aa,nn)*delta(bb,mm)
!	
!		call spec(Gamma,GD, GammaD)
!		call inv(GD, GD1)


	end subroutine Gamma_init

	
	subroutine Relaxation(RR, E, SS, S1, cor, RRD)
	
		complex(dpc), dimension(:,:), intent(out)		:: RR				! Relaxation matrix 
		complex(dpc), dimension(:,:), intent(in)		:: E, SS, S1, cor 	! Eigenenergies, Eigenvectors, site correlation
		complex(dpc), dimension(:,:), allocatable		:: G				! contains the transfer rates 
		complex(dpc), dimension(:,:), allocatable		:: GG				! sum of the rates over some states
		integer 										:: Ntot				! Dimension of total system (Nlev*nu vib)
		integer								:: aa, bb, cc, dd, ee, nn, mm	! indices used in the different loops
		complex(dpc)						:: dummy, dummy2				! dummy variable
		
		complex(dpc), dimension(:,:), allocatable		:: GD, GGD			! same as G and GG but using Debye spectral density
		complex(dpc), dimension(:,:), intent(out)		:: RRD				! Relaxation matrix 

		real(dp)	:: eps	! criteria used for relative error
		eps = 1.e-5
		
		
		Ntot = size(E,1)
!		write(6,*) 'Calculate relaxation tensor. Ntot = ', Ntot


!		write(6,*) '*****cor(exc)**'
!		write(6,200) transpose(real(cor))
		
		allocate(GG(Ntot*Ntot, Ntot*Ntot))
		allocate(G(Ntot, Ntot))
		
		GG = 0.
		G  = 0.

		allocate(GGD(Ntot*Ntot, Ntot*Ntot))
		allocate(GD(Ntot, Ntot))
		
		GGD = 0.
		GD  = 0.
		
		! Loop over the superindices of G
		do aa=1,Ntot
			do bb=1, Ntot 
 				do cc=1, Ntot
 					do dd=1, Ntot
						if (cc.ne.dd .and. pureDephasing=='yes') then	! remove pure dephasing terms from the redfield tensor and calculate them separately
							! Loop over the site vectors
 							do nn=1, Ntot
 								do mm=1, Ntot
									GG(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) = GG(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) + &
										S1(aa,nn)*SS(nn,bb)*S1(dd,mm)*SS(mm,cc) * cor(nn,mm)*Crl(E(cc,cc)-E(dd,dd))
									GGD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) = GGD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) + &
										S1(aa,nn)*SS(nn,bb)*S1(dd,mm)*SS(mm,cc) * cor(nn,mm)*CrD(E(cc,cc)-E(dd,dd))
								end do
	 						end do
						else ! include the pure dephasing terms in the redfield tensor
							! Loop over the site vectors
 							do nn=1, Ntot
 								do mm=1, Ntot
									GG(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) = GG(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) + &
										S1(aa,nn)*SS(nn,bb)*S1(dd,mm)*SS(mm,cc) * cor(nn,mm)*Crl(E(cc,cc)-E(dd,dd))
									GGD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) = GGD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) + &
										S1(aa,nn)*SS(nn,bb)*S1(dd,mm)*SS(mm,cc) * cor(nn,mm)*CrD(E(cc,cc)-E(dd,dd))
								end do
 							end do
						end if
!						if ( GG(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd).ne. GGD(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd)) &
!						& write(6,*) aa, bb, cc, dd, GG(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd), GGD(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd)
 					end do
 				end do
			end do
		end do

		
!--- Verification of detailed balance : GGDabab = GGDbaba * thermal(Ea-Eb)		
		do aa=1,Ntot
			do bb=1, Ntot 
				dummy =  GG(Ntot*(aa-1)+bb, Ntot*(aa-1)+bb) * thermal(E(aa,aa)-E(bb,bb)) 
				if (  (abs(GG(Ntot*(bb-1)+aa, Ntot*(bb-1)+aa) - dummy )/abs(dummy))> eps  ) then
					write(6,*) 'Check detailed balance of GG (should be 1:)', aa, bb, GG(Ntot*(bb-1)+aa, Ntot*(bb-1)+aa) / dummy
				end if
			end do
		end do
!---		
		
		
		do aa=1,Ntot
			do cc=1, Ntot
 				do ee=1, Ntot
					G(aa, cc) = G(aa, cc) + GG(Ntot*(ee-1)+aa, Ntot*(ee-1)+cc)
					GD(aa, cc) = GD(aa, cc) + GGD(Ntot*(ee-1)+aa, Ntot*(ee-1)+cc)
				end do
!				if (G(aa, cc).ne. GD(aa, cc)) write(6,*) aa, cc, G(aa, cc), GD(aa, cc)
			end do
		end do
		
		
!--- Calculate RR + output
		open(unit=200, file= trim(file_join(out_dir,"relax_vib.dat")))
		write(200,'(4a5, 10a15)') 'aa', 'bb', 'cc', 'dd', 'Ea', 'Eb','Ec','Ed','RRD', '1/RRD (tau)','GGD', 'RR', 'GD', 'GG'
		do aa=1,Ntot
			do bb=1, Ntot
 				do cc=1, Ntot
 					do dd=1, Ntot
!					if ((bb==aa.and.dd==cc) .or. (cc==aa.and.dd==bb)) then ! compute only the secular terms
						RR(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) =  0.5 * ( &						
							+ delta(aa,cc) * G(bb, dd) & 
							+ delta(bb,dd) * G(aa, cc) &
							- GG(Ntot*(aa-1)+cc, Ntot*(bb-1)+dd) &
							- GG(Ntot*(bb-1)+dd, Ntot*(aa-1)+cc) &
							)
						RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) =  0.5 * ( &						
							+ delta(aa,cc) * GD(bb, dd) & 
							+ delta(bb,dd) * GD(aa, cc) &
							- GGD(Ntot*(aa-1)+cc, Ntot*(bb-1)+dd) &
							- GGD(Ntot*(bb-1)+dd, Ntot*(aa-1)+cc) &
							)
						write(200,'(4i5, 20e15.3)') aa, bb, cc, dd, &
							real(E(aa,aa)), real(E(bb,bb)), real(E(cc,cc)), real(E(dd,dd)),&
							real(RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)),  &
							real(1/RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)),  &
							real(GGD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)),  &
							real(RR(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)),   &
							real(GD(aa,cc)),  real(GG(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)) , &
							imag(RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc))
!						if (abs(RR(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd)) .ne. abs(RRD(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd))) &
!						& write(6,*) RR(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd), RRD(Ntot*(aa-1)+bb, Ntot*(cc-1)+dd)
!					end if
					if ((doSecAp == 'yes') .and. ((bb.ne.aa.or.dd.ne.cc) .and. (cc.ne.aa.or.dd.ne.bb))) then
					! reset non-secular terms to zero to compute using the secular approximation
						RR(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) =   0.0
						RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc) =   0.0
					end if 
					end do
 				end do
			end do
		end do
		close(unit=200)
	

		! output in matrix format
		if (iterate==1) then
			open(unit=201, file= trim(file_join(out_dir,"relax_mat.dat")))
			write(201,*) '# real part '
			write(201,'(14e15.3)') transpose(real(RR))
			write(201,*) ''
			write(201,*) ''
			write(201,*) '# imaginary part'
			write(201,'(14e15.3)') transpose(imag(RR))
			close(unit=201)
		end if 

		
!--- Verify hermiticity of d rho_ab/dt (=> sum_cd abs(Rabcd)-abs(Rbacd) = 0 )
		do aa=1,Ntot
			do bb=1, Ntot 
				dummy = 0.
				do cc=1,Ntot
					do dd=1, Ntot 
 						dummy = dummy + abs(RRD(Ntot*(bb-1)+aa, Ntot*(dd-1)+cc)) - abs(RRD(Ntot*(aa-1)+bb, Ntot*(dd-1)+cc))
					end do
				end do
				if ( abs(dummy) > 1.e-7 ) then
					write(6,'(a55,2i4,1e15.3)') 'WARNING: a, b, sum_cd abs(Rabcd)-abs(Rbacd) .ne. 0',aa, bb, dummy
				endif
			end do
		end do
	

!---
		
		
!=== Verifications
		!-- Stationary solution: sum_c R_aacc exp(-Ec/kB T) = 0 (Eq 3.293b)
		do aa=1, Ntot
			dummy  = 0.
			do cc=1, Ntot
				dummy = dummy + RRD(Ntot*(aa-1)+aa, Ntot*(cc-1)+cc) * exp(-1.*E(cc,cc)/(kB*temp)) 
			end do 
			if (abs(dummy) > eps) write(6,*) 'Stationary solution (should be 0): ', dummy
		end do 
		
		!--- Raacc = Gac delta(ac) - GGcaca
		do aa=1,Ntot
			do cc=1,Ntot				
				if (aa.ne.cc) then
					dummy = RRD(Ntot*(aa-1)+aa, Ntot*(cc-1)+cc)+GGD(Ntot*(aa-1)+cc, Ntot*(aa-1)+cc)
					if (abs(dummy).ne.0.) write(6,*) 'WARNING!! aa = ', aa, '  cc = ', cc, 'Raacc+GGcaca (0)= ', dummy
				else
					dummy =  RRD(Ntot*(aa-1)+aa, Ntot*(cc-1)+cc)- (GD(aa,cc)-GGD(Ntot*(aa-1)+cc, Ntot*(aa-1)+cc))
					if (abs(dummy).ne.0.) write(6,*) 'WARNING!! aa = ', aa, '  cc = ', cc, 'Raaaa-(Gaa-GGcaca) (0)= ', dummy
				end if
			end do
		end do	
		
		!--- Rabab = 0.5 (Gbb + Gaa - GGaabb - GGbbaa) 
		do aa=1,Ntot
			do bb=1,Ntot				
				if (aa.ne.bb) then
					dummy = RRD(Ntot*(bb-1)+aa, Ntot*(bb-1)+aa)-0.5*(GD(bb,bb) + GD(aa,aa)-GGD(Ntot*(aa-1)+aa, Ntot*(bb-1)+bb)-GGD(Ntot*(bb-1)+bb, Ntot*(aa-1)+aa))
					if (abs(dummy).ne.0.) write(6,*) 'WARNING!! aa = ', aa, '  bb = ', bb, 'Rabab - 0.5 (Gbb + Gaa - GGaabb - GGbbaa) =', dummy
				end if
			end do
		end do	
		
!--- Verify that Raaaa = sum_ee.ne.aa k_aa,ee
		do aa=1,Ntot
			dummy = 0.
			do ee=1, Ntot
				if (ee.ne. aa) dummy = dummy + GGD(Ntot*(ee-1)+aa, Ntot*(ee-1)+aa)
			end do
			if (abs(dummy-RRD(Ntot*(aa-1)+aa, Ntot*(aa-1)+aa))/abs(dummy)>eps) write(6,*) 'WARNING: Raaaa', aa, dummy, ' should = ', RRD(Ntot*(aa-1)+aa, Ntot*(aa-1)+aa)
		end do
!---		
	
!===	
		
		deallocate(GG)
		deallocate(G)

		deallocate(GGD)
		deallocate(GD)
		
	!	allocate(G())
	
  200		format(18f10.4)			! print real/imag format for dim(matrix)=N*nu
	
	end subroutine Relaxation
	

	subroutine reversible_dynamics(LL,En)

		complex(dpc), dimension(:,:), intent(inout)	:: LL
		complex(dpc), dimension(:,:), intent(in)	:: En
		integer										:: Ntot
		integer										:: i, j, k, l

		Ntot= size(En,1)

		forall(i=1:Ntot, j=1:Ntot, k=1:Ntot, l=1:Ntot) &
		&		LL(Ntot*(i-1)+j,Ntot*(k-1)+l) = En(i,k)*delta(j,l) - delta(k,i)*En(l,j)


	end subroutine reversible_dynamics


	subroutine init_pops(rho, Amp, src, agg, SS, S1, DMavg)

		complex(dpc), dimension(:,:), intent(in)	:: SS, S1
		complex(dpc), dimension(:,:), intent(inout)	:: DMavg	! energy averaged transition dipole moments
		complex(dpc), dimension(:,:), intent(out)	:: rho
		complex(dpc), dimension(:,:,:), intent(inout)	:: Amp	! Amplitude of the signal after interaction with the 4 laser pulses
		complex(dpc), dimension(:,:), allocatable	:: DM	! transition dipole moment of the exciting light (matrix mu)
		complex(dpc), dimension(:,:), allocatable	:: ori	! orientation matrix in the excitonic basis
		complex(dpc), dimension(:,:), allocatable	:: ori_fac	! correcting factor for the average over an ensemble of 
																! randomly oriented molecules, [Hochstrasse2011a], excitonic basis
		integer, intent(in)							:: src	! position of the source site
		integer										:: pos	! position corrected with the number of vibrational levels
		type(aggregate), intent(in)					:: agg
		integer :: nn,nu,nun,num, mu,mum, mm, Ntot,ii
		integer	:: aa,bb,xi	! indices used in the excitonic basis
		integer	:: pos_g0		! position of the ground state eigenvector (n=g, nu=0)
		real :: val, norm
		real :: dd	! transition dipole moment
		complex(dpc)	:: dummy	! intermediate variable
	
		Ntot = size(rho,1)
		
		allocate(DM(Ntot,Ntot))
		allocate(ori_fac(Ntot,Ntot))
		allocate(ori(Ntot,3))
		DM      = 0.
		ori_fac = 0.
		ori     = 0.

!---	Compute the orientation factor to account for the averaging over an ensemble of randomly oriented molecules
!		Orientation matrix in the excitonic basis 
		do aa=1,Ntot
				do bb=1,Ntot
					ori(aa,:) = ori(aa,:) + SS(bb,aa) * agg%ori(bb,:)
				end do
!			end do
		end do
		! normalization
		do aa=1,Ntot
			norm=0.
			do ii=1, 3	!loop over the coordinate
				norm = norm+ori(aa,ii)**2
			end do
			if (norm.ne.0.) ori(aa,:)= ori(aa,:)/sqrt(norm)
		end do
!		
!		write(6,*) '***ori exc****'
!		write(6,'(3f10.4)') real(transpose(ori))

!		1/15 (2cos^2(theta_ab)+1)
		do aa=1,Ntot
			do bb=1,Ntot
				val=0.
				do ii=1,3
					val=val+ori(aa,ii)*ori(bb,ii)	! Scalar product between orientation of eigenvector a and b
				end do
				ori_fac(aa,bb)=1./15.*(1.+2.*val**2)
			end do
		end do

!		write(6,*) '***ori_fac****'
!		write(6,'(8f10.4)') real(transpose(ori_fac))
		
		
!		Modulus of the dipole moment for each BChl
		dd = 1.
		
		! define the matrix of transition dipole moment in the site basis
		if (oneParticle) then
			do nn=2, (agg%N-1)			! remove first and last site (optical coherences)
				do nu=1, agg%nu
					do xi=1, agg%nu
!						DM(1, agg%nu*(nn-1)+nu) = dd * FC((nu-1), agg%d(nn)) 
!						DM(agg%nu*(nn-1)+nu,1)  = dd * FC((nu-1), agg%d(nn)) 
						DM(xi, agg%nu*(nn-1)+nu) = dd * FC(xi, nu, agg%d(nn)) 
						DM(agg%nu*(nn-1)+nu,xi)  = dd * FC(xi, nu, agg%d(nn)) 
					end do
				end do
			end do
		else
			do nn=2, (agg%N-1)			! remove first and last site (optical coherences)
				do nun=1, agg%nu				! vib. level of the excited mol in its excited state
					do xi=1, agg%nu				! vib. level of the ground state molecule which will be excited
						do num=1, agg%nu 		! vibrational level of the coupled molecule which remains on its electronic ground state. 
							DM(agg%nu*(xi-1)+num, (agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nun-1)+num) = dd * FC(xi, nun, agg%d(nn)) 
							DM((agg%nu**(agg%N-2))*(nn-1)+agg%nu*(nun-1)+num,agg%nu*(xi-1)+num)  = dd * FC(xi, nun, agg%d(nn)) 
						end do 
					end do
				end do
			end do		
		endif

!		Transform the dipole moment matrix into the excitonic basis
		write(6,*) '***DM loc*****'
		write(6,'(16f10.4)') real(transpose(DM))
		DM=matmul(S1, matmul(DM,SS))
! 		write(6,*) '***DM exc*****'
! 		write(6,'(8f10.4)') real(transpose(DM))

		
		rho(1,1) = 1.
!		Transform the initial density matrix into the excitonic basis
		rho=matmul(S1, matmul(rho,SS))

!		write(6,*) '***rho exc*before excitation by 4pulses****'
!		write(6,'(8f10.4)') real(transpose(rho))

		if (signal) then
!---	Amplitude of the total response function after 4 interactions with light: ori_fac da da rho0 db db		
			!pos_g0=	Ntot-2*agg%nu+1	! WARNING: check that the position in correctly defined for more general (FMO) case (position depends on the energy of the eigenstates)
			pos_g0=	1	! always true if the ground state has much lower energy than the excited ones
			dummy = rho(pos_g0,pos_g0)	! save initial rho before rewritting over 
			select case (pathway)
			case ('R1','R2')
				do aa=1,Ntot
					do bb=1,Ntot
						rho(aa,bb)= ori_fac(aa,bb) * DM(aa,pos_g0)**2 * dummy * DM(bb,pos_g0)**2
						do xi=1,agg%nu	!loop over the ground state vibrational levels
							Amp(aa,bb,xi)= ori_fac(aa,bb) * DM(aa,pos_g0+(xi-1)) * DM(aa,pos_g0) * dummy * DM(bb,pos_g0) * DM(bb,pos_g0+(xi-1))
						end do
					end do
					do xi=1,agg%nu
						DMavg(aa, pos_g0+(xi-1)) = DMavg(aa, pos_g0+(xi-1)) + PDF*DM(aa,pos_g0+(xi-1))**2 
					end do
				end do
			case ('R4','R3')
				do aa=1,Ntot
					do bb=1,Ntot
						do xi=1,agg%nu
							Amp(aa,bb,xi) = ori_fac(aa,bb) * DM(bb,pos_g0) * DM(bb,pos_g0+(xi-1))* DM(aa,pos_g0+(xi-1)) * DM(aa,pos_g0)* dummy
							DMavg(aa, pos_g0+(xi-1)) = DMavg(aa,pos_g0+(xi-1)) + PDF*DM(aa, pos_g0+(xi-1))**2 
						end do
					end do 
				end do
			case ('OD')
				do aa=1,Ntot
					Amp(aa,pos_g0,:) = DM(aa,pos_g0)**2
				end do
			end select
			if (.false.)	then ! RFmono: monomer response function intensity in dimer calc.
				Amp = 0.
				Amp(3,3,1)= 1.   
				Amp(3,4,1)= 1.   
				Amp(4,3,1)= Amp(3,4,1)
				Amp(4,4,1)= 1.   
			end if	
			
			
! 			write(6,*) '========rho0=in init_pops========'
! 			write(6,'(8f10.4)') transpose(real(rho))

!			write(33,'(15f10.4)') Real(Amp(3,4,1)), Real(Amp(3,3,2)), Real(DM(3,1)), Real(DM(3,2)), Real(DM(4,1)), Real(DM(4,2))
!			write(6,*) '========Amp0,0====================='
!			write(6,'(8f10.4)') transpose(real(Amp(:,:,1)))
!			write(6,*) '========Amp0,1====================='
!			write(6,'(8f10.4)') transpose(real(Amp(:,:,2)))

		else
!---	Compute the coherent signal after two interaction with light. 
!			Initialization in the excitonic basis 
!       	WARNING: no correction from the orientational factor yet
			rho = matmul(DM, matmul(rho,DM))

! 		normalization
			do nn=1, Ntot
				norm=0
				do mm=1,Ntot
					norm = norm + rho(mm,nn)**2
				end do
				if (norm.ne.0.) rho(:,nn) = rho(:,nn)/sqrt(norm)
			end do

		end if
		
		

		! Excite all vib levels of the pigment
!		val=1./agg%nu
!		do nu=1, agg%nu
!			do mu=1, agg%nu
!				rho(agg%nu*(src-1)+nu, agg%nu*(src-1)+mu) = FC((nu-1), agg%d(src)) * FC((mu-1), agg%d(src))
!			end do
!		end do

		! Excite only ground level of the src pigment
!		pos = (src-1)*agg%nu+1
!		rho(pos, pos) = (1,0.)

!		Return result in the local basis
!		rho = matmul(SS,matmul(rho,S1))

!		write(6,*) '***rho loc****'
!		write(6,'(8f10.4)') real(transpose(rho))

		deallocate(DM)
		deallocate(ori_fac)
		deallocate(ori)

	end subroutine init_pops


	subroutine dynamics(rho_m_0, Amp_0, Udt, Gamma, rho_m_t, Amp_t, Tsteps, loc, SS, S1, agg, tmem)
		
		integer, parameter	:: Tini   = 0		! initial time step
 		integer, intent(in)	:: Tsteps			! number of time step

		logical , intent(in)	:: loc
		complex(dpc), dimension(:,:), intent(in)	:: SS, S1  	   ! transformation matrix and its inverse
!		complex(dpc), dimension(:,:), intent(inout) :: pop,coh 	   ! saved populations and coherences
		real(dp),     dimension(:),   intent(inout) :: tmem 	   ! saved time vector
		complex(dpc), dimension(:,:) 				:: Udt  	   ! evolution superoperator
		complex(dpc), dimension(:,:), intent(in)	:: Gamma 	   ! coeffiecent for the pure dephasing terms
		complex(dpc), dimension(:,:), intent(inout)	:: rho_m_0  	! Initial population matrix
		complex(dpc), dimension(:,:,:), intent(out) :: rho_m_t		! Population matrix at every time step
		complex(dpc), dimension(:,:), allocatable	:: rho_m		! Population matrix at one time step
		complex(dpc), dimension(:,:), allocatable	:: rho_pd		! Population matrix containing pure dephasing relaxation (rho_pd = rho_m * exp(-Gamma*g(t)))
		complex(dpc), dimension(:), allocatable 	:: rho_v		! Population matrix at one time step written in the vectorial form

		complex(dpc), dimension(:,:,:), intent(inout)	:: Amp_0  	! Amplitude of the signal at initial time
		complex(dpc), dimension(:,:,:), allocatable		:: Amp_m	! Amplitude at one time step 
		complex(dpc), dimension(:,:,:), allocatable		:: Amp_pd	! Amplitude at one time step containing the pure dephasing terms
		complex(dpc), dimension(:,:,:,:), intent(inout)	:: Amp_t  	! Amplitude for every time step

		complex(dpc), dimension(Tsteps)				:: xx			! store the evolution of the amount of coherence and of populations
		
		type(aggregate), intent(in)					:: agg
	
		integer										:: Ntot				! Dimension of total system (Ntot*nu vib)
		integer										:: aa, bb, xi 		! indices used in the different loops
		integer	:: tdt	! time indice

		Ntot = size(rho_m_0,1)

		call matrice_init(rho_m,Ntot)
		call matrice_init(rho_pd,Ntot)
		allocate(rho_v(Ntot*Ntot))
		allocate(Amp_m(Ntot,Ntot,agg%nu))
		allocate(Amp_pd(Ntot,Ntot,agg%nu))
		rho_v  = 0.
		Amp_m  = 0.
		Amp_pd = 0.

		rho_m = rho_m_0	! excitonic basis 
		Amp_m = Amp_0

!		Save initial conditions in rho_m_0
		if (loc) then
			rho_m_0 = matmul(SS, matmul(rho_m_0,S1))
		end if

!		Save initial conditions in rho_m_t(:,:,1)
		if (loc) then
			rho_m_t(1:Ntot,1:Ntot,1) = matmul(SS,matmul(rho_m_0,S1))
		else
			rho_m_t(1:Ntot,1:Ntot,1) = rho_m_0
			Amp_t(1:Ntot,1:Ntot,1:agg%nu,1) = Amp_0
		end if
		
!		write(6,*) '*********rho_m_0*********'
!		write(6,'(8f10.4)') transpose(real(rho_m_0))
!		write(6,*) '*********Amp-0*********'
!		write(6,'(8f10.4)') transpose(real(Amp_0(:,:,1)))
		
		tmem(1) = Tini
		
!		open(unit=222,file=trim(file_join(out_dir,"vib_aco.dat")))
!		open(unit=20, file=trim(file_join(out_dir, "goftHT.dat")))
!		write(20,'(3a15)') "#Time (fs)","Re[goftHT]","Im[goftHT]"

		do tdt=1,Tsteps
		
			tmem(tdt+1)= Tini + dt * tdt

!           Evolve rho_m (pure dephasing terms are added separately in rho_pd)
			rho_v = m2v(rho_m, Ntot)	! Transform into vector 
			rho_v = matmul(Udt,rho_v)	! Evolve populations by a step dt (in the excitonic basis)
			rho_m = v2m(rho_v,Ntot)		! Transforming back into matrix

!			Evolve the signal amplitude 
			do xi=1,agg%nu
				rho_v = m2v(Amp_m(:,:,xi),Ntot)
				rho_v = matmul(Udt,rho_v)
				Amp_m(:,:,xi) = v2m(rho_v,Ntot)
			end do
			
!			write(6,'(8f10.4)') Real(Udt(2,2))
		
!			write(6,*) '*********rho_m*********'
!			write(6,'(8f10.4)') transpose(real(rho_m))
!			write(6,*) '*********Amp_m_0********* tdt=', tdt 
!			write(6,'(8f10.4)') transpose(real(Amp_m(:,:,1)))
!			write(6,*) '*********Amp_m_1********* tdt=', tdt 
!			write(6,'(8f10.4)') transpose(real(Amp_m(:,:,2)))
		

!           Separately add pure dephasing terms, which are not accounted for in the redfield
			if (pureDephasing=='yes') then
				do aa=1,Ntot
					do bb=1,Ntot
						rho_pd(aa,bb)   = rho_m(aa,bb)   * exp( (-Gamma(aa,aa)-Gamma(bb,bb)+2.*Gamma(aa,bb)) * goftHT(tmem(tdt)) )
						if (pathway=='R1') then 
							Amp_pd(aa,bb,:) = Amp_m(aa,bb,:) * exp( (-Gamma(aa,aa)-Gamma(bb,bb)+2.*Gamma(aa,bb)) * goftHT(tmem(tdt)) )
						else
							Amp_pd(aa,bb,:) = Amp_m(aa,bb,:)
						end if
					end do
				end do
			else
				rho_pd = rho_m
				Amp_pd = Amp_m
			end if

!			output goftHT in dat file 
!			if (iterate<2) write(20,'(3e15.3)') tmem(tdt), goftHT(tmem(tdt))


!           Store results in rho_m_t at time (tdt+1) (tdt=1 corresponds to time=0)
			! WARNING: to be changed to account for vibrational levels
			! CHECK THE BOUNDARY AND IF IT IS CONSISTENT WITH THE ORDERING OF VECTORS. 
			! Save rho at current time step in rho_m_t according to the basis
			if (loc) then
				rho_m_t(1:Ntot,1:Ntot,tdt+1) = matmul(SS,matmul(rho_pd,S1))
			else
				rho_m_t(1:Ntot,1:Ntot,tdt+1) = rho_pd
				Amp_t(1:Ntot,1:Ntot,1:agg%nu,tdt+1) = Amp_pd
			end if

		end do

!		close(unit=222)
!		close(unit=20)
	
		call matrice_deinit(rho_m)
		call matrice_deinit(rho_pd)
		deallocate(rho_v)
		deallocate(Amp_m)
		deallocate(Amp_pd)
		
	end subroutine dynamics


	subroutine signal2D(w1, w3, Amp_t, En, RR, Gamma, agg, gft2)
	! compute the signal in the 2D spectrum for the frequencies w1 and w3, after fourier transform of the dynamical response function at the corresponding frequency. 
	
		complex(dpc), intent(in)					:: w1, w3		! frequency of the excited level after pulse 1 and 3, respectively
		complex(dpc)								:: dw1, dw3		! difference between the excited and optical frequency corresponding to pulse 1 and 3, respectively
!		integer, intent(in)							:: aa1, aa3		! excited and deexcited levels after interaction with pulses 1 and 3
		integer										:: opt1, opt3	! ground state level from/to (opt1/opt3) which the electronic state is excited/deexcites
		complex(dpc), dimension(:,:,:,:), intent(inout) :: Amp_t	! Time-dependent amplitude of the signal, complete model 
		complex(dpc), dimension(:,:), intent(in)	:: En			! Eigenvalue matrix
		complex(dpc), dimension(:,:), intent(in)	:: RR			! Relaxation tensor (=Redfield without the pure dephasing terms)
		complex(dpc), dimension(:,:), intent(in)	:: Gamma		! coefficient matrix
		complex(dpc), dimension(:), intent(inout)	:: gft2			! GFT(w1)*GFT(w3)
		type(aggregate), intent(in)	:: agg

		integer	:: aa, bb, xi, aa1	! used for internal loops

!		aa1  = agg%nu+1
		select case (pathway)
		case ('R1')
			do aa=1, Ntot	
				do bb=1, Ntot
					opt1 = 1
					dw1 = (En(aa,aa)-En(opt1,opt1)) - w1*cm2eV
					! Be careful with the indices between rho_ba and G_aa. Only the lower diagonal matrix in output in coh. 
!					if (bb==4 .and. aa==3) write(6,*) bb, aa, rho_m_t(bb,aa,1), GFT(dw1,aa,opt1,RR,Gamma)
!					rho_m_t(bb,aa,:) = rho_m_t(bb,aa,:) * GFT(dw1,aa,opt1,RR,Gamma) * GFT(dw3,aa,opt3,RR,Gamma)				
					do xi=1,agg%nu
						opt3 = xi
						dw3  = (En(aa,aa)-En(opt3,opt3))  - w3*cm2eV
						Amp_t(bb,aa,xi,:) = Amp_t(bb,aa,xi,:) * GFT(dw1,aa,opt1,RR,Gamma) * GFT(dw3,aa,opt3,RR,Gamma)
					end do 
				end do
				opt3 = 1
				dw3 = (En(aa,aa)-En(opt3,opt3)) - w3*cm2eV
				gft2(aa)= GFT(dw1,aa,opt1,RR,Gamma) * GFT(dw3,aa,opt3,RR,Gamma)
!				if (aa==3) write(6,*) 'Enaa, Enopt, RR33, GFT, dw1, dw3', En(aa,aa),En(opt1,opt1),RR(aa,aa), GFT(dw3,aa,opt3,RR,Gamma), dw1, dw3
			end do
		case ('R2','R3')
			do aa=1,Ntot
				do bb=1, Ntot
					opt1 = 1	! interval t1 always involves the vibrational ground state of the electronic ground state
					dw1  = -((En(aa,aa)-En(opt1,opt1)) - w1*cm2eV)
					do xi=1,agg%nu
						opt3 = xi
						dw3  = (En(bb,bb)-En(opt3,opt3)) - w3*cm2eV
						Amp_t(aa,bb,xi,:) = Amp_t(aa,bb,xi,:)*GFT(dw1,aa,opt1,RR,Gamma)*GFT(dw3,bb,opt3,RR,Gamma)
					end do
 					if (aa==bb) gft2(aa)= GFT(dw1,aa,opt1,RR,Gamma) * GFT(dw3,bb,opt3,RR,Gamma) ! output only diagonal terms
! 					if (aa==bb .and. aa==4) write(6,*) 'Enaa, Enopt, RR33, GFT, dw1, dw3', En(aa,aa),En(opt1,opt1),RR(aa,aa),GFT(dw3,bb,opt3,RR,Gamma), dw1, dw3
				end do
			end do		
		case ('R4')
			do aa=1,Ntot
				do bb=1, Ntot
					opt1 = 1	! interval t1 always involves the vibrational ground state of the electronic ground state
					opt3 = 1	! interval t3 always involves the vibrational ground state of the electronic ground state
					dw1  = (En(aa,aa)-En(opt1,opt1)) - w1*cm2eV
					dw3  = (En(bb,bb)-En(opt3,opt3)) - w3*cm2eV
					Amp_t(aa,bb,:,:) = Amp_t(aa,bb,:,:)*GFT(dw1,aa,opt1,RR,Gamma)*GFT(dw3,bb,opt3,RR,Gamma)
!					write(6,*) 'aa,Amp_t(1000)', aa, Amp_t(1,aa,:,1000)
 					if (aa==bb) gft2(aa)= GFT(dw1,aa,opt1,RR,Gamma) * GFT(dw3,bb,opt3,RR,Gamma) ! output only diagonal terms
!					if (aa==bb .and. aa==3) write(6,*) 'Enaa, Enopt, RR33, GFT, dw1, dw3', En(aa,aa),En(opt1,opt1),RR(aa,aa),GFT(dw3,bb,opt3,RR,Gamma), dw1, dw3
				end do
			end do		
		end select

		
!		do aa=1,Ntot
!			do bb=1,Ntot
!				opt3 = opt1
!				dw3  = w3*cm2eV - (En(aa3,aa3)-En(opt3,opt3))
!				rho_m_t(aa,bb,:) = rho_m_t(aa,bb,:) * GFT(dw1,aa1,opt1,RR,Gamma) * GFT(-dw3,aa3,opt3,RR,Gamma)
!				do xi=1,agg%nu
!				opt3 = opt1+xi-1
!					dw3  = w3*cm2eV - (En(aa3,aa3)-En(opt3,opt3))
!					Amp_t(aa,bb,xi,:) = Amp_t(aa,bb,xi,:) * GFT(dw1,aa1,opt1,RR,Gamma) * GFT(-dw3,aa3,opt3,RR,Gamma)
!				end do 
!			end do
!		end do

	end subroutine signal2D
	
	
	complex(dpc) function GFT(dw, aa, opt, RR, Gamma)
!	One-side Fourier transform of the goftHT
	
		integer,      intent(in)	:: aa						! eigenstate involved in the optical transition
		integer,      intent(in)	:: opt						! postition of the ground state from which the electronic state is excited
		complex(dpc), intent(in)	:: dw						! frequency of interest in the 2D spectrum - frequency of the coherence
		complex(dpc), dimension(:,:), intent(in)	:: RR		! Relaxation tensor (=Redfield without the pure dephasing terms)
		complex(dpc), dimension(:,:), intent(in)	:: Gamma	! coefficient matrix
		complex(dpc)	:: cst, dummy, dumGamma, dumIncGamma 	! intermediate variable 
		complex(dpc)	:: relax	! relaxation rate calculated from the tensor of the optical coherence of interest
		real(dp)		:: wD,lambda
		integer			:: Ntot
		
		Ntot   = INT(sqrt(size(RR,1)*1.))
		wD     = wDin*cm2ev
		lambda = lambdain * cm2eV
		relax  = RR(Ntot*(opt-1)+aa,Ntot*(opt-1)+aa)
		cst    = Gamma(aa,aa)*(2*lambda*kB*temp/wD**2 - (0,1)*lambda/wD)	!constant from the goftHT
		dummy  = (relax + (0,1)*dw)/wD + cst 
		

		call cgamma(0,dummy,dumGamma)	! dumGamma = Gamma(dummy), Gamma function for complex argument
		dumIncGamma = cincgamma(dummy, cst)	! Incomplete gamma function

		GFT = Exp(cst)*(dumGamma - dumIncGamma) / ((cst**dummy)*wD) ! -i factor comes from the one-sided fourier transform (see Eq 5.30g [Mukamelbook])
		
!		if (opt==2) write(6,*) 'Ntot,opt, aa, Rindex,relax, Gammaaa,cst,dw,dummy, dumIncGamma, GFT', Ntot,opt, aa, Ntot*(opt-1)+aa,relax,Gamma(aa,aa),cst, dw, dummy, dumIncGamma, GFT
!		if (iterate==1) write(6,*) 'aa, dw, relax, gamma, GFT, dumGamma, dumIncGamma = ', aa, dw, relax, Gamma(aa,aa), GFT , dumGamma, dumIncGamma
		
	end function


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
		integer, intent(in)							:: N	! number of levels
		integer										:: i

		forall(i=1:N) v2m(i,1:N) = v(N*(i-1)+1:N*i)
	
	end function


	
	complex(dpc) function SpD0(w)				! Normalized low-frequency Spectral Density [Eq. 26 Adolphs_BJ_2006_91_2778]
		
		complex(dpc), intent(in)	:: w
		real(dp)					:: wD, j0, lambda
		real(dp)					:: s1, s2, w1, w2
		
		s1 = 0.8_dp
		s2 = 0.5_dp
		w1 = 0.069E-3_dp	! in eV
		w2 = 0.24E-3_dp		! in eV
		
		wD = wDin * cm2eV
		lambda =  lambdain * cm2eV

		j0 = 2*wD*lambda/pi
		
		if (abs(w) == 0.) then
!			SpD0 = j0/wD**2
			SpD0 = 8*lambda/wD*(kB*temp)
		else
			SpD0 = w**3/((s1+s2)*2*factorial(7)) * ( s1/w1**4*exp(-sqrt(abs(w/w1))) + s2/w2**4*exp(-sqrt(abs(w/w2)))) 
		end if
			
	end function


	complex(dpc) function SpD(w)				! Total Spectral Density, include explicit treatment on 1 high-frequency model [Eq. 25 Adolphs_BJ_2006_91_2778]
		
		complex(dpc), intent(in)	:: w
		real(dp)					:: wD, wH,j0, lambda
		real(dp)					:: S0, SH	! Huang-Rys factors of pigment-protein coupling and high-frequency mode, resp.
		
		wD = wDin * cm2eV

		lambda =  lambdain * cm2eV

		j0 = 2*wD*lambda/pi
		
		S0 = 0.35_dp	! WARNING to be checked with Fig. 2
		SH = 0.22_dp
		
		if (abs(w) == 0.) then
!			SpD = j0/wD**2
			SpD = 8*lambda/wD*(kB*temp)
!		else if (w == wH) then		! explicit treament of one single mode
!			SpD = S0 * SpD0(w) + SH
!			write(6,*) 'you are in the high frequency mode and SpD = ', SpD
		else
			SpD = S0 * SpD0(w)  
		end if
			
	end function
	
	
	complex(dpc) function Crl(w)
		
		complex(dpc), intent(in)	:: w
		real(dp)            		:: wD, wH, j0, lambda
		
		wD = wDin * cm2eV

		lambda =  lambdain * cm2eV

		j0 = 2*wD*lambda/pi

		if (abs(w) == 0.) then
!			Crl = j0/wD**2
			Crl = 8*lambda/wD*(kB*temp)
		else
			Crl = 2.*pi * w**2*( 1. + nthermal(w))* (SpD(w) - SpD(-1._dp*w))
		end if
				
	end function
	
	
	complex(dpc) function nthermal(w)
	  ! Bose-Einstein distribution function

		complex(dpc), intent(in)	:: w
		
		nthermal = 1._dp / (exp(w/(kB*temp))-1._dp)

	end function	
	
	
	complex(dpc) function thermal(w)
	! Bose-Einstein distribution function
	
		complex(dpc), intent(in)	:: w

		thermal = exp(w/(kB*temp))

	end function


	pure integer function delta(aa,bb)        
		
		integer, intent(in)	:: aa, bb
		
		if (aa==bb) then
			delta = 1
		else
			delta = 0
		end if
		
	end function
	
	
! For verifications

	complex(dpc) function CrD(w)
	
		complex(dpc), intent(in)	:: w
		real(dp)					:: wD,j0, lambda
		
		wD = wDin * cm2ev
		lambda =  lambdain * cm2eV

!		j0 = 2*wD*lambda/pi

		if (abs(w) == 0.) then	! only used to compute pure dephasing terms (wcd=0)
!			CrD = j0in/wD**2
 			CrD = 4*lambda*kB*temp/wD
!			CrD =0
!			write(6,*) 'CrD(0)=', CrD 
		else
			CrD = 2*pi*( 1. + nthermal(w))* (w2J(w) - w2J(-1.*w))
		end if
			
	end function
	
	
	complex(dpc) function w2J(w)
	! w**2 . J(w) using Debye spectral density

		real(dp)					:: wD,j0, lambda
		complex(dpc), intent(in)	:: w

		wD = wDin*cm2ev
		lambda =  lambdain * cm2eV
		
		j0 = 2*wD*lambda/pi

		if (real(w) >= 0.) then
			w2J = 2*wD*lambda/pi*w/(w*w + wD*wD)
		else
			w2J = (0., 0.)
		endif

	end function
	
	
	real(dp) function FC(nu, mu, dd) !--- nu=1 represents the vibrational quantum zero
	! Calculate Franck-Condon factor with the wave function overlap
	! after application of the shift operator XX = -sqrt(S) (a^+ - a) 
	! Output: <g_nu|e_mu> between ground and excited states

	! Alternative to compute <g_0|e_mu> is to use: [Spano_JCP_116_2002_5877]
	! FC = sqrt(  dd**(2*nu) / (2.**nu * factorial(nu)) * exp(-1.*dd**2./2.)  )
		
		real(dp), intent(in)	:: dd
		integer, intent(in)		:: nu, mu
		complex(dpc), dimension(:,:), allocatable	:: Aan, Acr	!Annihilation and creation operators
		complex(dpc), dimension(:,:), allocatable	:: XX, Xn, XD, XD1	! solution matrix, eigenstates, transformation matrix and its inverse
		integer	:: nn, mm 	! loop indice
		integer :: dim	! dimension of the matrices
		
		dim = Ntot*3	!compute a matrix of dimension larger than what is needed to have precise results
		call matrice_init(Aan, dim)     
		call matrice_init(Acr, dim)     
		call matrice_init(XX, dim)	   
		call matrice_init(Xn, dim)	   
		call matrice_init(XD, dim)	   
		call matrice_init(XD1, dim) 
		
		! Definition of the annihilation Aan and creation Acr operator
		do nn=1, dim-1	
			Aan(nn, nn+1) = sqrt(1.*nn)	! Aan |nn> = sqrt(nn) |nn-1>
			Acr(nn+1, nn) = sqrt(1.*nn)	! Acr |nn> = sqrt(nn+1) |nn+1>
		end do
		
!		write(6,*) 'd= ', dd
!		write(6,*) ' =======Aan========= '
!		write(6,'(24f10.5)') transpose(real(Aan))
!		write(6,*) ' =======Acr========= '
!		write(6,'(24f10.5)') transpose(real(Acr))
!		write(6,*) '  Aan(3,4) = ', Aan(3,4)
!		write(6,*) '  Acr(2,1) = ', Acr(2,1)

		! Compute the shift operator <g_nu|XX|g_mu> = <g_nu|e_mu>
		forall(nn=1:dim, mm=1:dim) XX(nn,mm) = -dd/sqrt(2.)*(Acr(nn,mm)-Aan(nn,mm))


!		write(6,*) ' ====== Shift op. ========= ' ! identical to AA[S0] in mathematica
!		write(6,'(24f10.5)') transpose(real(XX))

		! Transform in the diagonal basis 
		call spec2(XX,XD,Xn)
		call inv(XD,XD1)
		
		! Solution of the eq. diff after application of the shift operator
		XX = 0.
		forall(nn=1:dim) XX(nn,nn) = exp(Xn(nn,nn) * (1.,0.))
		
		! Transformation back into the initial basis (eg)
		XX = matmul(XD, matmul(XX,XD1))

!		write(6,*) ' ====== XX before transposition========= '
!		write(6,'(24f10.5)') transpose(real(XX))

		
		! output the results for <g_nu|e_mu> in XX(nu+1, mu+1) if nu, mu start from zero
		!XX = conjg(transpose(XX))
		!WARNING: somehow, there is a transposition during matmul(XD, matmul(XX,XD1))!!!... 
		! We should not transpose again to have the results as is mathematica!... but where is the bug???? 
		
		FC = real(XX(nu, mu)) !<g_nu|e_mu>

!		write(16,*) ' ====== XX ========= '
!		write(16,'(24f10.5)') transpose(real(XX))


		call matrice_deinit(Aan)	   
		call matrice_deinit(Acr)	   
		call matrice_deinit(XX) 	   
		call matrice_deinit(Xn) 	   
		call matrice_deinit(XD) 	   
		call matrice_deinit(XD1)		
	
		
		
	
	end function


    !--- inspired from excitons.F90
	!
    ! Creates energy shifts due to the disorder
    !
    subroutine gap_disorder(agg)
        type(aggregate)			:: agg
		integer 				:: i,j

!        if (allocated(ran)) then
!             DEALLOCATE(ran)
!        end if
!        ALLOCATE(ran,((agg%N*agg%nu)**2))

        call random_Normal(ran)

!        if (allocated(den)) then
!            DEALLOCATE(den)
!        end if
!        ALLOCATE(den,((agg%N*agg%nu),(agg%N*agg%nu)))

!		dgap = agg%dwidth*ran(2)
		dgap = agg%dwidth*ran(1)
		
		write (600,*) iterate, dgap, ran


    end subroutine gap_disorder


	! Compute the line shape function goft in the limit of high temperature and using debye spectral density
	! (Eq. 8.49 from [MukamelBook])
	
	complex(dpc) function goftHT(time)

		real(dp)				:: wD,lambda
		real(dp), intent(in)	:: time

		wD = wDin*cm2ev
		lambda =  lambdain * cm2eV
		
		goftHT = (exp(-(time*wD)) + time*wD - 1.)*(2*lambda*kB*temp/wD**2 - ((0,1)*lambda/wD))

	end function
	

	subroutine relax_time(RR, Gamma, Tsteps)
	
		integer, intent(in)							:: Tsteps	! number of time step
		real(dp), dimension(Tsteps, Ntot*(Ntot-1)/2):: tau		! save relaxation time for elements of rho
		real(dp), dimension(Tsteps, Ntot, Ntot)		:: tau_m	! save relaxation time for elements of rho
		complex(dpc), dimension(:,:), intent(in)	:: RR		! Relaxation matrix
		complex(dpc), dimension(:,:), intent(in)	:: Gamma	! Coefficient for the correlation fct in excitonic basis
		real(dp)									:: wD, lambda
		
		integer		:: aa,bb,tdt,kk! used in iterative loops
		real(dp)	:: xsi		! intermediate variable: xsi_ab = Gamma_aa + Gamma_bb - 2 Gamma_ab
		real(dp)	:: cst		! constant, intermediate result
		real(dp)	:: rrab		! coefficient of the relaxation matrix 

		wD = wDin*cm2ev
		lambda =  lambdain * cm2eV
	
		do tdt=1,Tsteps
			do aa=1,Ntot
				do bb=1,Ntot	
					xsi = Gamma(aa,aa)+Gamma(bb,bb)-2.*Gamma(aa,bb)
					cst = (2*lambda*kB*temp/wD**2)
					rrab= real( RR( Ntot*(bb-1)+aa,Ntot*(bb-1)+aa ) )	! RRabab
					tau_m(tdt, aa, bb) = tdt/(  rrab*tdt+xsi*cst*(exp(-(tdt*wD)) + tdt*wD - 1.)  )
				end do
			end do
		end do

		tau=(0.,0.)
		kk = 0
		do aa = 1,Ntot
			do bb = 1,Ntot
				if (aa>bb) then
					kk = kk + 1
					if ( (aa<Ntot) .and. (bb>=1) ) then
						tau(1:Tsteps,kk) = tau_m(:,aa,bb)
					end if
				end if
			end do
		end do
			
		open(unit=21,file=trim(file_join(out_dir,"relax_time.dat")))	

		do tdt = 1,Tsteps
			write(21,*) tdt, tau(tdt,:)
		end do

		close(unit=21)

	end subroutine relax_time
	
	
end module qme_vibronic_networks

