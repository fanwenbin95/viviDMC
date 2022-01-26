module input
    implicit none

    contains

    subroutine get_input
        implicit none
        
        call log_copyright
        call get_geom
        call get_para
        call allocate_para

    end subroutine


    ! read parameters from 'desk.inp'
    subroutine get_para
        use constants, only: dp
        use para, only: N0, maxNwalker, dt, Nstep, Neq, &
                        sample_type
        implicit none
        character(128) :: inputFileName
        logical inputFileStatus

        integer :: Nmax
        real(dp) :: Neq_ratio
        character(1) :: sample_type_char

        inputFileName = 'desk.inp'
        inquire(file=trim(inputFileName), exist=inputFileStatus)
        if (inputFileStatus) then
            open(375, file=trim(inputFileName), status='old', action='read')
            ! desk -> dsk -> 375
        else
            write(*, *) 'The parameter file "desk.inp" did not exist! '
            stop
        end if

        ! initial number of walkers
        read(375, *)
        read(375, *) N0

        ! maximum number of walkers
        read(375, *)
        read(375, *) Nmax
        if (Nmax >= N0) then
            maxNwalker = Nmax
        else
            maxNwalker = Nmax * N0
        end if

        ! time step
        read(375, *)
        read(375, *) dt

        ! number of steps
        read(375, *)
        read(375, *) Nstep

        ! number of equilibrium steps
        read(375, *)
        read(375, *) Neq_ratio
        if (Neq_ratio < 0d0 ) Neq_ratio = 0d0
        if (Neq_ratio < 1d0 ) then
            Neq = int(Nstep * Neq_ratio)
        else
            Neq = int(Neq_ratio)
        end if

        ! sampling method
        read(375, *)
        read(375, *) sample_type_char
        sample_type = 1
        if ( trim(sample_type_char) == 'v' .or. trim(sample_type_char) == 'V') then
            sample_type = 1
        else
            write(*, *) 'The sampling method not supported! '
            write(*, *) 'Defaultly use the vibrational sampling! '
        end if

        ! sample options
        if (sample_type == 1) then
            call read_vib_option(375)
        end if

        close(375)

    end subroutine


    subroutine read_vib_option(unit)
        use constants, only: dp, pi
        use para, only: Nmode, eleMass, &
                        vib_scale, Ngrid, vibrc, vibwf, vibcdf, &
                        vibpot, kmat, pmat, &
                        vibrand, N0
        use units, only: amu2me
        implicit none
        integer, intent(in) :: unit ! file id
        real(dp) :: space
        integer :: i, j

        read(unit, *)
        read(unit, *) vib_scale
        read(unit, *)
        read(unit, *) Ngrid

        ! vibrational reaction coordinate s
        ! * the s, including both ends, would be with length Ngrid+1
        allocate( vibrc(0:Ngrid) )
        space = 2 * vib_scale / dble(Ngrid)
        do i = 0, Ngrid
            vibrc(i) = - vib_scale + i * space
        end do

        ! potential energy in each mode
        allocate( vibpot(0:Ngrid, Nmode) )

        ! initialize kinetic matrix
        ! REF: D. T. Colbert and W. H. Miller, The Journal of Chemical Physics 96, 1982 (1992).
        allocate( kmat(0:Ngrid, 0:Ngrid) )
        do i = 0, Ngrid
            do j = 0, Ngrid
                if ( i == j ) then
                    kmat(i, j) = pi * pi / 3d0
                else
                    kmat(i, j) = 2d0 / ( i - j ) / ( i - j )
                end if
                kmat(i, j) = kmat(i, j) * (-1d0) ** (i-j)
            end do
        end do
        kmat = kmat / 2d0 / space / space / sum(eleMass) / amu2me! note hbar === 1

        allocate( pmat(0:Ngrid, 0:Ngrid) )
        pmat = 0d0

        allocate( vibwf(0:Ngrid, Nmode) )
        allocate( vibcdf(0:Ngrid, Nmode) )
        allocate( vibrand(N0, Nmode) )

    end subroutine


    subroutine allocate_para
        use constants, only: dp
        use para, only: N0, maxNwalker, Nstep, Neq, & 
                        dt, alpha, sigma, &
                        batch, Nframe, &
                        Nwalker, q, Eref, w, &
                        Natom, eleMass
        implicit none
        integer :: i
        integer :: batch_size

        allocate( Nwalker(Nstep) )
        Nwalker = -1
        Nwalker(1) = N0
        alpha = 0.5d0 / dt
        sigma = dsqrt(dt / sum(eleMass))

        allocate( q(3, Natom, maxNwalker) )
        q = -1d16

        allocate( Eref(Nstep) )
        Eref = 0d0

        allocate( w(maxNwalker) )
        w = 0d0
        w(1:N0) = 1d0

        ! split all walkers into Nframe groups
        ! * starting location is batch(i-1)+1
        ! * ending location is batch(i)
        allocate( batch(0:Nframe) )
        batch_size = floor( float( N0 / Nframe ) )
        batch(0) = 0
        do i = 1, Nframe - 1
            batch(i) = i * batch_size
        end do
        batch(Nframe) = N0

        write(*, '(" Initial number of walker :    ", I8)') Nwalker(1)
        write(*, '(" Max number of walker :        ", I8)') maxNwalker
        write(*, '(" time-step (a.u.) :            ", F12.6)') dt
        write(*, '(" Number of total steps :       ", I8)') Nstep
        write(*, '(" Number of equilibrium steps : ", I8)') Neq
        write(*, '(" alpha === 1 / (2 dt) :        ", F12.6)') alpha
        write(*, '(" sigma === sqrt(dt/M) :        ", F12.6)') sigma

    end subroutine allocate_para


    ! read geometries from 'EQ.xyz'
    subroutine get_geom
        use constants, only: dp
        use para, only: Natom, Nframe, eleName, q0
        use units, only : au2ang
        implicit none
        character(128) :: inputFileName
        logical inputFileStatus
        ! character(1024) commandLine
        integer i, j
        integer NatomCheck
        integer lineStatus

        inputFileName = 'EQ.xyz'
        inquire(file=trim(inputFileName), exist=inputFileStatus)
        if (inputFileStatus) then
            open(478, file=trim(inputFileName), status='old', action='read')
        else
            write(*, *) 'Input file "EQ.xyz" not exists! '
            stop
        end if

        ! read number of geometries `Nframe`
        read(478, *) Natom ! the first number of atom
        rewind(478)
        Nframe = 0
        lineStatus = 0
        do while (lineStatus == 0)
            read(478, *, iostat=lineStatus) NatomCheck
            if (lineStatus /= 0) exit
            ! read following geometry if Natom is equal to the first one
            if (NatomCheck == Natom) then
                Nframe = Nframe + 1
                read(478, *)
                do i = 1, Natom
                    read(478, *)
                end do
            ! finish reading file if not equals
            else
                lineStatus = -1
            end if
        end do
        rewind(478)
        
        call allocate_geometry

        ! read all geometries
        do j = 1, Nframe
            read(478, *)
            read(478, *)! commandLine
            do i = 1, Natom
                read(478, *) eleName(i), q0(1:3, i, j)
            end do
        end do

        close(478)
        write(*, '(" Number of geometries : ", I3)') Nframe
        
        q0 = q0 / au2ang ! to Bohr

        ! allocate all values in program
        call get_mass
        call get_mass_mat

        ! call initial_energy(q(1:3, 1:Natom, 1))

    end subroutine


    subroutine allocate_geometry
        use para
        implicit none

        allocate( eleName(Natom) )
        allocate( eleMass(Natom), massMat(3*Natom, 3*Natom) )
        allocate( q0(3,Natom,Nframe), grad(3,Natom), hess(3*Natom, 3*Natom) )

        ! number of vibraional modes
        ! * only suitable for non-linear molecule
        Nmode = 3*Natom - 6
        Nbond = Natom * (Natom - 1) / 2

    end subroutine allocate_geometry


    subroutine initial_energy(q)
        use constants, only: dp
        use para, only: Natom
        use units, only: au2wn, au2ang
        implicit none
        real(dp), intent(in) :: q(3, Natom)
        real(dp) :: v, g(3, Natom), h(3*Natom, 3*Natom), f(3*Natom)

        write(*, '(A)') ' Initial geometry (Bohr) :'
        write(*, '(3F12.6)') q
        call pot(q, v)
        write(*, '(A, F12.8)') ' Energy (Hartree) :', v

        call calc_grad(q, g)
        write(*, '(A)') ' Gradient (Hartree/Bohr) :'
        write(*, '(3F12.6)') g

        call calc_hess(q, h)
        write(*, '(A)') ' Hessian (per cm/Ang^2) :'
        write(*, '(6F12.2)') h * au2wn * au2ang
        ! write(*, '(A)') ' Hessian (Hartree/Bohr^2) :'
        ! write(*, '(6F12.2)') h

        call vib_ana(h, f)
        write(*, '(A)') ' Frequency (per cm) :'
        write(*, '(6F12.4)') f * au2wn

    end subroutine


    subroutine get_mass
        use constants, only: dp
        use para, only: Natom, eleName, eleMass
        use elements, only: supNum, allName, allMass
        implicit none
        integer :: atom, ele

        do atom = 1, Natom
            do ele = 1, supNum
                if ( trim(eleName(atom)) == trim(allName(ele))) then
                    eleMass(atom) = allMass(ele)
                    exit
                end if
            end do
        end do

    end subroutine get_mass


    ! prepare mass matrix for vibrational analysis
    ! - squared-rooted for the mass matrix
    subroutine get_mass_mat
        use constants, only: dp
        use para, only: Natom, eleMass, massMat
        implicit none
        integer :: i, j
        
        do i = 1, Natom
            do j = 1, Natom
                massMat( (i-1)*3+1:i*3, (j-1)*3+1:j*3 ) = eleMass(i) * eleMass(j)
            end do
        end do
        massMat = dsqrt(massMat)

    end subroutine get_mass_mat


    subroutine log_copyright
    implicit none

        write(*, '(A)') ' viviDMC'
        write(*, '(A)') ''
        write(*, '(A)') ' - Author : Wenbin FAN (wbfan21@m.fudan.edu.cn)'
        write(*, '(A)') ' - Date : Jan. 20, 2020'
        write(*, '(A)') ''
        write(*, '(A)') ' * Atomic unit in program. '
        write(*, '(A)') ''
        

    end subroutine log_copyright

end module