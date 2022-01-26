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
        call get_ave_energy(q(:,:,1:Nwalker(1)), Eref(1))
        Erefnew = Eref(1)

        do step = 1, Nstep
            Nwalker(step) = Npoint
            Eref(step) = Erefnew
            call random_move

            call born_kill(Eref(step), Erefnew)
            write(*, *) step, Npoint, Erefnew*au2wn

            if (Npoint < int(0.05 * Nwalker(1)) ) then
                write(*, '("Too small walkers. STOP! ", I8)') Npoint
                stop
            end if
        end do

        do j = 1, Npoint
            write(999, *) Natom
            write(999, *)
            do i = 1, Natom
                write(999, '(I3, 3F12.6)') i, q(:, i, j)
            end do
        end do

    end subroutine propagate


    subroutine born_kill(Erefin, Erefout)
        implicit none
        real(dp), intent(in) :: Erefin
        real(dp) :: v(Npoint)
        real(dp) :: qNew(3, Natom, maxNwalker)
        real(dp), intent(out) :: Erefout
        integer :: i, count

        do i = 1, Npoint
            call pot(q(:,:,i), v(i))
        end do
        w(1:Npoint) = dexp( - (v(:) - Erefin) * dt ) * w(1:Npoint)
        write(88, '(F12.6)') w(1:Npoint)
        stop

        count = 0
        Erefout = 0d0
        do i = 1, Npoint
            if ( w(i) >= 0.1d0 ) then
                count = count + 1
                qNew(:, :, count) = q(:, :, i)
                Erefout = Erefout + v(i)
                if ( w(i) > 20d0 ) then
                    count = count + 1
                    qNew(:, :, count) = q(:, :, i)
                    Erefout = Erefout + v(i)
                end if
            end if
        end do
        q = -1d16
        q(:,:,1:count) = qNew(:,:,1:count)

        Npoint = count
        Erefout = dot_product( v, w ) / sum(w) - alpha * &
                    dlog( sum(w) / dble(Nwalker(1)) )

    end subroutine


    subroutine random_move
        use constants, only: pi
        implicit none
        real(dp) :: u1(3, Natom, Npoint), u2(3, Natom, Npoint)

        call random_number(u1)
        call random_number(u2)
        q(:,:,1:Npoint) = q(:,:,1:Npoint) + &
            sigma * dsqrt( -2d0 * dlog(u1) ) * dcos( 2d0 * pi * u2 )

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