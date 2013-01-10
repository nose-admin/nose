module resp3_rss

	use std_types

    real(dp), dimension(:,:), allocatable :: rhoL, rhoR          ! initial window density matrix
    real(dp), dimension(:,:), allocatable :: dge, deg,dfe        ! transition dipole moments
    real(dp), dimension(:,:), allocatable :: oeg,ofe             ! transition frequencies
    real(dp), dimension(:),   allocatable :: pp                  ! equilibrium populations

    complex(dpc), dimension(:,:,:,:,:), pointer :: Ugege => NULL()
    complex(dpc), dimension(:,:,:,:,:), pointer :: Uegeg => NULL()
    complex(dpc), dimension(:,:,:,:,:), pointer :: Ueeee => NULL()
    complex(dpc), dimension(:,:,:,:,:), pointer :: Ugggg => NULL()
    complex(dpc), dimension(:,:,:,:,:), pointer :: Ufefe => NULL()
    complex(dpc), dimension(:,:,:,:,:), pointer :: Uefef => NULL()


    complex(dpc), dimension(:,:), allocatable :: R1g0, R2g0, R3g0, R4g0
    complex(dpc), dimension(:,:), allocatable :: R1f0, R2f0, R3f0, R4f0

    complex(dpc), dimension(:,:,:), allocatable :: W1g,W1f,W2g,W2f,W3g,W3f,W4g,W4f
    complex(dpc), dimension(:,:,:), allocatable :: D1g,D1f,D2g,D2f,D3g,D3f,D4g,D4f

    complex(dpc), dimension(:,:,:,:,:), allocatable :: extD

    integer :: Ng, Ne, Nf
    integer :: Nt1, Nt3
    integer :: Nbe, Nbf, Np

    real(dp) :: ddt1

    logical :: window_doorway

end module resp3_rss
