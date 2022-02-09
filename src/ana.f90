module ana
    use constants, only:dp
    use para, only: Eref, Neq, Nstep
    use units, only: au2wn
    use para, only: eleMass
    implicit none
    
contains

    subroutine analyze
    implicit none

    ! average energy
    write(*, *) 'Average reference energy (per cm) : ', & 
        sum(Eref(Neq:Nstep)) / (Nstep - Neq) * au2wn
    
    end subroutine analyze
    
end module ana