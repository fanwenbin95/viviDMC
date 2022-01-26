! GEAR : global reaction maps

! Program notes:
! - All values are in atomic units

program main
    use input
    use init
    use prpgt
    implicit none
    
    call get_input
    call walker_init
    call propagate

end program main


