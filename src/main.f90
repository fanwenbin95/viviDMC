! GEAR : global reaction maps

! Program notes:
! - All values are in atomic units

program main
    use input
    use init
    use prpgt
    use ana
    implicit none
    
    call get_input
    call pgm_init
    call walker_init
    call propagate
    call analyze

end program main


