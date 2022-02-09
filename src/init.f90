module init
    use constants, only: dp
    use para, only: Natom, Nframe, batch, q0, q, Nwalker
    ! use para, only: q => q0
    implicit none
    
contains

    subroutine pgm_init
        implicit none

        call random_seed()
        call pot_initial()

    end subroutine

    subroutine walker_init
        implicit none
        integer :: frame, i

        do frame = 1, Nframe
            ! change the initialization way as you want
            ! * leave a switch on desk file, 
            !   like 'v' for vib, 'r' for completely random, ...
            call vib_init(frame)
        end do

        ! check that walker coordinates have been sampled in code level
        do i = 1, Nwalker(1)
            if ( sum(q(:,:,i)) < -100d0 ) then
                write(*, &
                    '("The coordinate of ", I8, &
                    "th walker has not been sampled! ")') i
                write(*, '(3F12.6)') q(:,:,i)
            end if
        end do

    end subroutine walker_init


    ! initialize the walker in a vibrational way
    subroutine vib_init(ith)
        use units, only: au2wn
        use para, only: Nmode, Ngrid, vibrc, vibpot, &
                        kmat, pmat, vibwf, vibcdf, &
                        vibrand
        implicit none
        integer, intent(in) :: ith
        real(dp) :: c(3, Natom)
        real(dp) :: freq(3*Natom), evec(3*Natom, 3*Natom) ! for results of vibrational analysis
        integer :: i
        integer :: mode, grid
        ! temporary variable
        real(dp) :: v, c0(3, Natom)
        real(dp) :: dvr_energy(0:Ngrid), dvr_wf(0:Ngrid,0:Ngrid)
        integer :: gs ! index of DVR ground state

        c(1:3, 1:Natom) = q0(1:3, 1:Natom, ith)
        ! do vibrational analysis
        ! * all frequencies and vectors has been ranked by the absolute value of frequencies
        call q_vib_ana(c, freq, evec)
        ! normalize the un-mass-weighted eigenvectors
        ! * clearly the vibrational modes are in dimension 2, say A(:, i), which was verified by Gaussian 16's output
        call matrix_norm(3*Natom, evec, 2, evec)
        write(*, *) ith
        do i = 1, 3*Natom
            write(*, '(F12.3, 9F7.2)') freq(i)*au2wn, evec(:,i)
        end do

        ! first Nmode modes do not include 3 trans and 3 rot modes usually
        do mode = 1, Nmode
            ! scan the whole vib energy
            c0 = q0(1:3, 1:Natom, ith)
            do grid = 0, Ngrid
                c = c0 + reshape( evec(:,mode), (/3,Natom/) ) * vibrc(grid)
                call pot(c, v)
                ! v = ( vibrc(grid) ) **2 * 2d0
                vibpot(grid, mode) = v
                pmat(grid, grid) = v
                write(86, *) mode, vibrc(grid), v*au2wn
            end do
            ! DVR
            call DEEV(Ngrid+1, kmat+pmat, dvr_energy, dvr_wf)

            ! write dvr
            gs = minloc(dvr_energy, dim=1) - 1
            ! "- 1" deals with the first 0 index
            ! * which trapped me for one day. QAQ
            do grid = 0, Ngrid
                write(87, '(I3,5F13.8)') mode, vibrc(grid), vibpot(grid, mode), dvr_wf(grid, gs)
            end do
            ! write(*, '(I3,F13.6,I5,4F13.6)') mode, freq(mode)*au2wn, gs, dvr_energy(gs)*au2wn
            vibwf(0:Ngrid, mode) = dvr_wf(0:Ngrid, gs) ** 2
        end do

        ! get cummulative distribution function CDF
        vibcdf(0, 1:Nmode) = vibwf(0, 1:Nmode)
        do i = 1, Ngrid
            vibcdf(i, 1:Nmode) = vibcdf(i-1, 1:Nmode) + vibwf(i, 1:Nmode)
        end do
        ! normalize cdf to [0, 1]
        do mode = 1, Nmode
            vibcdf(:, mode) = ( vibcdf(:, mode) - vibcdf(0, mode) ) / &
                              ( vibcdf(Ngrid, mode) - vibcdf(0, mode) )
        end do

        ! assign random vibrational initial coordinates
        call random_number( vibrand )
        do i = batch(ith-1)+1, batch(ith)
            q(:, :, i) = q0(:, :, ith)
            do mode = 1, Nmode
                vibrand(i, mode) = &
                    vibrc( &
                        minloc( &
                            abs(vibcdf(:, mode) - vibrand(i, mode)), dim=1 &
                            ) &
                        )
                q(:, :, i) = q(:, :, i) + & 
                    reshape(evec(:, mode), (/3, Natom/) ) * vibrand(i, mode)
            end do

            ! write(999, *) Natom
            ! write(999, *)
            ! do mode = 1, Natom
            !     write(999, '(I3, 3F12.6)') mode, q(:, mode, i)
            ! end do
            !! The vibrational coordinates have been checked. ! 2022-01-26 16:57:28 Wenbin FAN @FDU
        end do

    end subroutine

end module init