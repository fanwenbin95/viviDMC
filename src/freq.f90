subroutine vib_ana(hess, freq, V)
    use constants, only: dp, pi
    use para, only: Natom, massMat
    use units, only: au2J, au2ang, ang2m, ang2m, amu2kg, clight, au2wn
    implicit none
    real(dp), intent(in) :: hess(3*Natom, 3*Natom)
    real(dp), intent(out) :: freq(3*Natom), V(3*Natom, 3*Natom)
    real(dp) :: E(3*Natom)
    integer :: i, j, bubble
    real(dp) :: freq_temp, V_temp(3*Natom)

    call DEEV(3*Natom, hess/massMat, E, V)

    freq = E * au2J / ang2m**2 / au2ang**2 / amu2kg ! to SI unit
    do i = 1, 3*Natom
        if ( freq(i) >= 0d0 ) then
            freq(i) = dsqrt(freq(i))
        else
            freq(i) = -dsqrt(abs(freq(i)))
        end if
    end do
    freq = freq / 2d0 / pi / clight / 100d0 ! to per cm
    freq = freq / au2wn ! to atomic unit ! Hartree

    V = V * massMat ! unscale the mass weight

    ! sort frequencies by absolute value using bubble sort
    i = 3*Natom
    do while ( i > 1 )
        bubble = 0
        do j = 1, i-1
            if ( abs(freq(j)) < abs(freq(j+1)) ) then
                freq_temp = freq(j)
                V_temp = V(:, j)

                freq(j) = freq(j+1)
                V(:,j) = V(:,j+1)

                freq(j+1) = freq_temp
                V(:,j+1) = V_temp

                bubble = j
            end if
        end do
        i = bubble
    end do

end subroutine vib_ana


! a more convenient routine for evaluating frequency from coordinates
subroutine q_vib_ana(q, freq, V)
    use constants, only: dp
    use para, only: Natom
    implicit none
    real(dp), intent(in) :: q(3, Natom)
    real(dp), intent(out) :: freq(3*Natom), V(3*Natom, 3*Natom)
    real(dp) :: hess(3*Natom, 3*Natom)

    call calc_hess(q, hess)
    call vib_ana(hess, freq, V)

end subroutine


! ! remove extra modes of translation and rotation
! subroutine remove_modes(f_in, v_in, f_out, v_out)
!     use constants, only: dp
!     use para, only: Natom, Nmode
!     implicit none
!     real(dp), intent(in) :: f_in(3*Natom, 3*Natom), v_in(3*Natom, 3*Natom)
!     real(dp), intent(out) :: f_out(3*Natom, Nmode), v_out(3*Natom, 3*Nmode)
!     integer :: mode

!     ! rank the frequencies 
!     f_out = 0d0
!     v_out = 0d0

    

! end subroutine