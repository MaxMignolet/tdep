submodule(type_phonon_dispersions) type_phonon_dispersions_qpoint_generalized_operators

use konstanter, only: lo_freqtol, r8
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_blas_lapack_wrappers, only: lo_zheev
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix

implicit none
contains

!> calculate generalized group velocities taking degeneracies into account
subroutine return_generalized_group_velocity(ompoint,p,fc,qpoint,genvel,mem)
    !> harmonic properties at this q
    class(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> qpoint
    class(lo_qpoint), intent(in) :: qpoint
    !> generalized group velocities (xyz,mode,mode)
    complex(r8), dimension(:,:,:), intent(out) :: genvel
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    complex(r8), dimension(:,:,:), allocatable :: Dq
    complex(r8), dimension(:,:), allocatable :: m0,m2
    integer :: ialpha,b1,b2,nb
    logical, parameter :: use_syms = .false.

    nb=p%na*3

    allocate(m0(nb,nb))
    allocate(m2(nb,nb))
    allocate(Dq(nb,nb,3))
    m0=0.0_r8 ! dynamical matrix
    Dq=0.0_r8 ! q-gradient of the dynamical matrix
    call fc%dynamicalmatrix(p, qpoint, m0, mem, Dq, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])
    ! qdirection only relevant when at Gamma for non-analycities
    m0(:,:) = ompoint%egv(:,:)

    genvel(:,:,:) = 0.0_r8 ! zero initialisation otherwise acoustic modes are not set
    do ialpha=1,3
        m2=0.0_r8
        do b1=1,nb
        do b2=1,nb
            m2(b1,b2)=m2(b1,b2)+dot_product(m0(:,b1),matmul(Dq(:,:,ialpha),m0(:,b2)))
        enddo
        enddo
        ! m2 holds the component ialpha of the sandwiched gradient of the dyn mat

        ! Store as generalized group velocity
        do b1=1,nb
        do b2=1,nb
            if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
            if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
            genvel(ialpha,b1,b2) = m2(b1,b2) / (2*sqrt(ompoint%omega(b1)*ompoint%omega(b2)))
        enddo
        enddo

    enddo

    deallocate(m2)
    deallocate(m0)
    deallocate(Dq)

end subroutine

!> returns the \mathbb{L}_{\alpha} thingy I defined
subroutine return_generalized_angmom(ompoint,p,qpoint,Lalpha)
    !> harmonic properties at this q
    class(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> qpoint
    class(lo_qpoint), intent(in) :: qpoint
    !> generalized group velocities (xyz,mode,mode)
    complex(r8), dimension(:,:,:), intent(out) :: Lalpha


    integer, dimension(3,3,3) :: lcv
    complex(r8), dimension(:,:), allocatable :: m2
    complex(r8), dimension(3) :: cw0,cw1
    integer :: j,ialpha,ibeta,igamma,b1,b2,nb
    logical, parameter :: use_syms = .false.


    ! levi-cevita tensor
    lcv=0
    lcv(1,2,3)=1
    lcv(2,3,1)=1
    lcv(3,1,2)=1
    lcv(1,3,2)=-1
    lcv(2,1,3)=-1
    lcv(3,2,1)=-1

    nb=p%na*3

    allocate(m2(nb,nb))
    do ialpha=1,3
        m2=0.0_r8
        do b1=1,nb
        do b2=1,nb
            if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
            if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle

            do j=1,p%na
                cw0=conjg(ompoint%egv( (j-1)*3+1:j*3,b1 ))
                cw1=      ompoint%egv( (j-1)*3+1:j*3,b2 )
                do ibeta=1,3
                do igamma=1,3
                    m2(b1,b2)=m2(b1,b2) - lo_imag * lcv(ialpha,ibeta,igamma)*cw0(ibeta)*cw1(igamma)
                enddo
                enddo
            enddo
        enddo
        enddo

        Lalpha(ialpha,:,:)=m2
    enddo

    deallocate(m2)

end subroutine

end submodule
