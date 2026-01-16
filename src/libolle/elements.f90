!> calculate generalized group velocities taking degeneracies into account
subroutine lo_return_generalized_group_velocity(p,fc,qpoint,ompoint,genvel,mem)
        !> structure
        type(lo_crystalstructure), intent(in) :: p
        !> forceconstant
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        !> qpoint
        class(lo_qpoint), intent(in) :: qpoint
        !> harmonic properties at this q
        type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
        !> generalized group velocities (xyz,mode,mode)
        real(r8), dimension(:,:,:), intent(out) :: genvel
        !> memory helper
        type(lo_mem_helper), intent(inout) :: mem

        complex(r8), dimension(:,:,:), allocatable :: Dq
        complex(r8), dimension(:,:), allocatable :: m0,m1,m2,m3
        real(r8), dimension(:), allocatable :: v0,v1
        real(r8), dimension(3) :: w0,w1
        real(r8) :: f0
        logical, dimension(:), allocatable :: modefixed
        integer :: i,j,ii,jj,ialpha,b1,b2,nb,iop

        nb=p%na*3

        allocate(m0(nb,nb))
        allocate(Dq(nb,nb,3))
        m0=0.0_r8
        Dq=0.0_r8
        call fc%dynamicalmatrix(p, qpoint, m0, mem, Dq, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])

        allocate(modefixed(nb))
        allocate(v1(nb))
        v1=0.0_r8
        modefixed=.false.
        do ialpha=1,3
            ! Find the properly rotated eigenvectors
            modefixed=.false.
            v1=0.0_r8
            m0=0.0_r8
            do b1=1,nb
                if ( modefixed(b1) ) cycle

                if ( ompoint%omega(b1) .lt. lo_freqtol ) then
                    ! Don't care for acoustic modes
                    m0(:,b1)=0.0_r8
                    modefixed(b1)=.true.
                elseif ( ompoint%degeneracy(b1) .eq. 1 ) then
                    ! Not degenerate, nothing to worry about
                    m0(:,b1)=ompoint%egv(:,b1)
                    modefixed(b1)=.true.
                else
                    allocate(m1(ompoint%degeneracy(b1),ompoint%degeneracy(b1)))
                    allocate(m2(nb,ompoint%degeneracy(b1)))
                    allocate(m3(nb,ompoint%degeneracy(b1)))
                    allocate(v0(ompoint%degeneracy(b1)))
                    m1=0.0_r8
                    m2=0.0_r8
                    m3=0.0_r8
                    v0=0.0_r8

                    do i=1,ompoint%degeneracy(b1)
                    do j=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        jj=ompoint%degenmode(j,b1)
                        m1(i,j)=dot_product(ompoint%egv(:,ii),matmul(Dq(:,:,ialpha),ompoint%egv(:,jj)))
                    enddo
                    enddo
                    ! Figure out the mean eigenvalue
                    call lo_zheev(m1,v0,jobz='V')
                    f0=sum(v0)/real(ompoint%degeneracy(b1),r8)
                    ! Rotate eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        m2(:,i)=ompoint%egv(:,ii)
                    enddo
                    m3=matmul(m2,m1)
                    ! Store correctly rotated eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        b2=ompoint%degenmode(i,b1)
                        m0(:,b2)=m3(:,i)
                        v1(b2)=f0
                        modefixed(b2)=.true.
                    enddo

                    deallocate(m1)
                    deallocate(m2)
                    deallocate(m3)
                    deallocate(v0)
                endif
            enddo
            ! At this point m0 should hold the appropriately rotated eigenvectors.
            ! We could possibly do even better if we also average over the small
            ! group of q.

            ! allocate(m1(nb,nb))
            ! allocate(m2(nb,nb))
            ! m1=0.0_r8
            ! m2=0.0_r8
            ! do i=1,qpoint%n_invariant_operation
            !     iop=qpoint%invariant_operation(i)
            !     call lo_eigenvector_transformation_matrix(rotmat,p%rcart,qpoint%r,p%sym%op(iop))
            !     m1 = matmul(rotmat,m0)
            !     do b1=1,nb
            !     do b2=1,nb
            !         m2(b1,b2)=m2(b1,b2)+dot_product(m1(:,b1),matmul(Dq(:,:,ialpha),m1(:,b2)))
            !     enddo
            !     enddo
            ! enddo
            ! m2=m2/real(qpoint%n_invariant_operation,r8)

            ! Sandwich with rotated eigenvectors
            allocate(m2(nb,nb))
            m2=0.0_r8
            do b1=1,nb
            do b2=1,nb
                m2(b1,b2)=dot_product(m0(:,b1),matmul(Dq(:,:,ialpha),m0(:,b2)))
            enddo
            enddo
            ! Make sure the degeneracies hold
            do b1=1,nb
                if ( ompoint%degeneracy(b1) .eq. 1 ) cycle
                do i=1,ompoint%degeneracy(b1)
                    b2=ompoint%degenmode(i,b1)
                    if ( b1 .eq. b2 ) then
                        m2(b1,b1)=v1(b1)
                    else
                        m2(b1,b2)=0.0_r8
                    endif
                enddo
            enddo
            ! Store as generalized group velocity
            do b1=1,nb
            do b2=1,nb
                if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
                if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
                genvel(ialpha,b1,b2)=real( m2(b1,b2)/(ompoint%omega(b1)+ompoint%omega(b2)), r8)
            enddo
            enddo

            deallocate(m2)
        enddo

        ! Maybe a final average over the small group? Yes no maybe. Seems reasonable.
        do b1=1,nb
        do b2=1,nb
            if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
            if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
            w0=genvel(:,b1,b2)
            w1=0.0_r8
            do i=1,qpoint%n_invariant_operation
                iop=qpoint%invariant_operation(i)
                if ( iop .gt. 0 ) then
                    w1=w1+matmul(p%sym%op(iop)%m,w0)
                else
                    w1=w1-matmul(p%sym%op(iop)%m,w0)
                endif
            enddo
            genvel(:,b1,b2)=w1/real(qpoint%n_invariant_operation)
        enddo
        enddo
end subroutine

!> returns the \mathbb{L}_{\alpha} thingy I defined
subroutine lo_return_angmom_Lalpha_thingy(p,qpoint,ompoint,Lalpha)
        !> structure
        type(lo_crystalstructure), intent(in) :: p
        !> qpoint
        class(lo_qpoint), intent(in) :: qpoint
        !> harmonic properties at this q
        type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
        !> generalized group velocities (xyz,mode,mode)
        complex(r8), dimension(:,:,:), intent(out) :: Lalpha


        integer, dimension(3,3,3) :: lcv
        complex(r8), dimension(:,:,:), allocatable :: Dq
        complex(r8), dimension(:,:), allocatable :: m0,m1,m2,m3,rotmat
        real(r8), dimension(:), allocatable :: v0,v1
        complex(r8), dimension(3) :: cw0,cw1
        real(r8), dimension(3) :: w0,w1
        real(r8) :: f0
        logical, dimension(:), allocatable :: modefixed
        integer :: i,j,ii,jj,ialpha,ibeta,igamma,b1,b2,nb,iop

        ! Will use the levi-cevita tensor since I have that
        ! in my notes. Not the most efficient but it's fine.
        lcv=0
        lcv(1,1,1)=0
        lcv(1,1,2)=0
        lcv(1,1,3)=0
        lcv(1,2,1)=0
        lcv(1,2,2)=0
        lcv(1,2,3)=1
        lcv(1,3,1)=0
        lcv(1,3,2)=-1
        lcv(1,3,3)=0
        lcv(2,1,1)=0
        lcv(2,1,2)=0
        lcv(2,1,3)=-1
        lcv(2,2,1)=0
        lcv(2,2,2)=0
        lcv(2,2,3)=0
        lcv(2,3,1)=1
        lcv(2,3,2)=0
        lcv(2,3,3)=0
        lcv(3,1,1)=0
        lcv(3,1,2)=1
        lcv(3,1,3)=0
        lcv(3,2,1)=-1
        lcv(3,2,2)=0
        lcv(3,2,3)=0
        lcv(3,3,1)=0
        lcv(3,3,2)=0
        lcv(3,3,3)=0

        ! Hmm. For the group velocities I could sandwich the gradient of the dynamical
        ! matrix between eigenvectors to find the correct rotation. I wonder what the
        ! equivalent matrix should be for the angular momentum stuff? I think it's just
        ! the generator thingy:
        !
        ! Dq_{i\alpha,j\beta,\gamma) = \delta_{ij} \eps_{\alpha\beta\gamma}
        !
        ! At least 60% sure. Better than nothing.
        nb=p%na*3
        allocate(Dq(nb,nb,3))
        Dq=0.0_r8
        do i=1,p%na
            do ialpha=1,3
            do ibeta=1,3
            do igamma=1,3
                ii=(i-1)*3+ialpha
                jj=(i-1)*3+ibeta
                Dq(ii,jj,igamma)=lcv(ialpha,ibeta,igamma)
            enddo
            enddo
            enddo
        enddo

        allocate(m0(nb,nb))
        allocate(rotmat(nb,nb))
        allocate(modefixed(nb))
        allocate(v1(nb))
        m0=0.0_r8
        v1=0.0_r8
        modefixed=.false.
        do ialpha=1,3
            ! Find the properly rotated eigenvectors
            modefixed=.false.
            v1=0.0_r8
            m0=0.0_r8
            do b1=1,nb
                if ( modefixed(b1) ) cycle

                if ( ompoint%omega(b1) .lt. lo_freqtol ) then
                    ! Don't care for acoustic modes
                    m0(:,b1)=0.0_r8
                    modefixed(b1)=.true.
                elseif ( ompoint%degeneracy(b1) .eq. 1 ) then
                    ! Not degenerate, nothing to worry about
                    m0(:,b1)=ompoint%egv(:,b1)
                    modefixed(b1)=.true.
                else
                    allocate(m1(ompoint%degeneracy(b1),ompoint%degeneracy(b1)))
                    allocate(m2(nb,ompoint%degeneracy(b1)))
                    allocate(m3(nb,ompoint%degeneracy(b1)))
                    allocate(v0(ompoint%degeneracy(b1)))
                    m1=0.0_r8
                    m2=0.0_r8
                    m3=0.0_r8
                    v0=0.0_r8

                    do i=1,ompoint%degeneracy(b1)
                    do j=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        jj=ompoint%degenmode(j,b1)
                        m1(i,j)=dot_product(ompoint%egv(:,ii),matmul(Dq(:,:,ialpha),ompoint%egv(:,jj)))
                    enddo
                    enddo
                    ! Figure out the mean eigenvalue
                    call lo_zheev(m1,v0,jobz='V')
                    f0=sum(v0)/real(ompoint%degeneracy(b1),r8)
                    ! Rotate eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        ii=ompoint%degenmode(i,b1)
                        m2(:,i)=ompoint%egv(:,ii)
                    enddo
                    m3=matmul(m2,m1)
                    ! Store correctly rotated eigenvectors
                    do i=1,ompoint%degeneracy(b1)
                        b2=ompoint%degenmode(i,b1)
                        m0(:,b2)=m3(:,i)
                        v1(b2)=f0
                        modefixed(b2)=.true.
                    enddo

                    deallocate(m1)
                    deallocate(m2)
                    deallocate(m3)
                    deallocate(v0)
                endif
            enddo

            ! At this point m0 should hold the appropriately rotated eigenvectors.
            ! I'm like 50% sure at least. Better than nothing at least.
            ! In addition to that we can average over the small group,
            ! because I'm a grown up and I do whatever I want.


            allocate(m1(nb,nb))
            allocate(m2(nb,nb))
            m1=0.0_r8
            m2=0.0_r8
            do i=1,qpoint%n_invariant_operation
                iop=qpoint%invariant_operation(i)
                call lo_eigenvector_transformation_matrix(rotmat,p%rcart,qpoint%r,p%sym%op(iop))
                m1 = matmul(rotmat,m0)

                do b1=1,nb
                do b2=1,nb
                    if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
                    if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle

                    do j=1,p%na
                        cw0=conjg(m1( (j-1)*3+1:j*3,b1 ))
                        cw1=m1( (j-1)*3+1:j*3,b1 )
                        do ibeta=1,3
                        do igamma=1,3
                            m2(b1,b2)=m2(b1,b2) + lcv(ialpha,ibeta,igamma)*cw0(ibeta)*cw1(igamma)
                        enddo
                        enddo
                    enddo
                enddo
                enddo
            enddo
            m2=m2/real(qpoint%n_invariant_operation,r8)
            ! and add the frequency factors
            do b1=1,nb
            do b2=1,nb
                if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
                if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
                m2(b1,b2) = m2(b1,b2) / sqrt( ompoint%omega(b1)/ompoint%omega(b2) )
            enddo
            enddo
            ! Not sure how to store .. m2 should be purely imaginary at this point. Possibly.
            Lalpha(ialpha,:,:)=m2

            deallocate(m1)
            deallocate(m2)
        enddo

        ! Maybe a final average over the small group? Yes no maybe. Seems reasonable.
        ! Commented out for now because symmetry is complicated.
        do b1=1,nb
        do b2=1,nb
            if ( ompoint%omega(b1) .lt. lo_freqtol ) cycle
            if ( ompoint%omega(b2) .lt. lo_freqtol ) cycle
            cw0=Lalpha(:,b1,b2)
            cw1=0.0_r8
            do i=1,qpoint%n_invariant_operation
                iop=qpoint%invariant_operation(i)
                if ( iop .gt. 0 ) then
                    ! pseudovector so det(op) in front.
                    f0=lo_determ(p%sym%op(iop)%m)
                    w1=w1+f0*matmul(p%sym%op(iop)%m,w0)
                else
                    ! this means time reversal, not sure what to do. Skip for now.
                    w1=w1+w0
                endif
            enddo
            Lalpha(:,b1,b2)=w1/real(qpoint%n_invariant_operation)
        enddo
        enddo
end subroutine
