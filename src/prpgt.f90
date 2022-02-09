module prpgt
    use constants, only: dp
    use para, only: Nstep, &
                    Natom, q, Nwalker, maxNwalker, &
                    Eref, w, &
                    dt, alpha, sigma
    implicit none
    integer :: Npoint ! current number of walkers

    contains

    subroutine propagate
        use units, only: au2wn
        implicit none
        integer :: step
        real(dp) :: Erefnew
        integer :: i, j

        ! initialize the Eref
        Npoint = Nwalker(1)
        call get_ave_energy(q(:,:,1:Npoint), Erefnew)
        ! write(*, *) 0, Npoint, Erefnew*au2wn

        do step = 1, Nstep
            Nwalker(step) = Npoint
            Eref(step) = Erefnew

            call random_move
            call born_kill(Eref(step), Erefnew)

            write(*, '(I5,I8,2F12.3)') step, Npoint, Erefnew*au2wn, sum(w(1:Npoint))
            write(90, '(I5,I8,2F12.3)') step, Npoint, Erefnew*au2wn, sum(w(1:Npoint))

            if (Npoint < 3 ) then
                write(*, '("Too few walkers. STOP! ", I8)') Npoint
                stop
            end if
        end do

        do j = 1, Npoint
            write(999, *) Natom
            write(999, *)
            do i = 1, Natom
                write(999, '(I3, 3F12.6)') i*2, q(:, i, j)
            end do
        end do

    end subroutine propagate


    subroutine born_kill(Erefin, Erefout)
        implicit none
        real(dp), intent(in) :: Erefin
        real(dp) :: v(maxNwalker)
        real(dp) :: qNew(3, Natom, maxNwalker)
        real(dp) :: wNew(maxNwalker)
        real(dp), intent(out) :: Erefout
        integer :: i, count
        real(dp) :: nume, denom ! numerator and denominator for Eref
        real(dp) :: u(Npoint)

        do i = 1, Npoint
            call pot(q(:,:,i), v(i))
        end do
        ! w(1:Npoint) = dexp( - (v(1:Npoint) - Erefin) * dt ) * w(1:Npoint)
        w(1:Npoint) = -1d0 + dexp( - (v(1:Npoint) - Erefin) * dt )

        count = 0
        nume = 0d0
        denom = 0d0
        call random_number( u )
        do i = 1, Npoint
            ! if ( w(i) >= 0.1d0 ) then
            !     count = count + 1
            !     qNew(:, :, count) = q(:, :, i)
            !     wNew(count) = w(i)
            !     nume = nume + v(i) * w(i)
            !     denom = denom + w(i)
            !     if ( w(i) > 20d0 ) then
            !         count = count + 1
            !         qNew(:, :, count) = q(:, :, i)
            !         wNew(count) = w(i) * 0.5d0
            !         wNew(count-1) = wNew(count)
            !         ! nume = nume + v(i) * w(i)
            !         ! denom = denom + w(i)
            !     end if
            ! end if
            !! below commentted codes are another simple method
            !! which does not depend on the weights
            if ( v(i) < Erefin ) then
                if ( u(i) < w(i) ) then
                    count = count + 2
                    qNew(:,:,count-1) = q(:,:,i)
                    qNew(:,:,count) = q(:,:,i)
                    nume = nume + v(i) * 2
                else
                    count = count+1
                    qNew(:,:,count) = q(:,:,i)
                    nume = nume + v(i)
                end if
            else
                if ( u(i) < -w(i) ) then
                    continue
                else
                    count = count + 1
                    qNew(:,:,count) = q(:,:,i)
                    nume = nume + v(i)
                end if
            end if
            if ( count >= maxNwalker - 2 ) exit
        end do

        q(:,:,1:count) = qNew(:,:,1:count)
        w(1:count) = wNew(1:count)

        Npoint = count
        ! Erefout = nume / denom - alpha * dlog( denom / dble(Nwalker(1)) )
        Erefout = nume / dble(count) - alpha * dble(count-Nwalker(1)) / dble(Nwalker(1))

    end subroutine


    subroutine random_move
        use constants, only: pi
        implicit none
        real(dp) :: u1(3, Natom, Npoint), u2(3, Natom, Npoint)
        integer :: i

        call random_number(u1)
        call random_number(u2)
        do i = 1, Natom
            q(:,i,1:Npoint) = q(:,i,1:Npoint) + &
                sigma(i) &
                * dsqrt( -2d0 * dlog(u1(:,i,1:Npoint)) ) & 
                * dcos( 2d0 * pi * u2(:,i,1:Npoint) )
        end do
    end subroutine


    subroutine get_ave_energy(q, ave_enrg)
        implicit none
        real(dp), intent(in) :: q(3, Natom, Npoint)
        real(dp), intent(out) :: ave_enrg
        integer :: i
        real(dp) :: v(Npoint)

        do i = 1, Npoint
            call pot(q(:,:,i), v(i))
        end do
        ave_enrg = sum(v) / Npoint

    end subroutine

end module