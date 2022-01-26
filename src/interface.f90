subroutine pot_initial
    implicit none
    ! Initialize your PES here
end subroutine


subroutine pot(q, v)
    use constants, only: dp
    use para, only: Natom, Nbond
    use units, only: au2ang
    implicit none
    real(dp), intent(in) :: q(3, Natom)
    real(dp), intent(out) :: v
    real(dp) :: r(Nbond)

    ! hook your PES here, please. 
    ! The units are in atomic always. 

    call cart2bond(q, r)
    r(3) = dacos( (r(1)*r(1) + r(2)*r(2) - r(3)*r(3)) / (2d0*r(1)*r(2)) )
    call vibpot(r, v, 1)
    ! Note that the author commented all detailed outputs of PES. 

end subroutine pot


! Replace the Hessian or gradient routine if they are analytical for your PES. 
! * programmingly the q0 avoids changing the value of coordinates
subroutine calc_hess(q0, h)
    use constants, only: dp, fd, quadru_fd_square
    use para, only: Natom
    implicit none
    real(dp), intent(in) :: q0(3, Natom)
    real(dp), intent(out) :: h(3*Natom, 3*Natom)
    real(dp) :: q(3, Natom), vpp, vpm, vmm, vmp ! plus, minus
    integer i, j, k, l ! xyz, Natom, xyz, Natom for loops
    integer nx, ny ! one index for xyz and Natom

    q = q0
    do i = 1, Natom ! outer is y
        do j = 1, 3
            do k = 1, Natom ! inner loop is x
                do l = 1, 3
                    ! naming two index for Hessian matrix
                    nx = (k-1)*3+l
                    ny = (i-1)*3+j
                    ! evaluate if element is in upper triangle
                    if ( nx >= ny ) then ! 
                        ! + +
                        q(j, i) = q(j, i) + fd
                        q(l, k) = q(l, k) + fd
                        call pot(q, vpp)
                        ! + -
                        q(l, k) = q(l, k) - fd - fd
                        call pot(q, vpm)
                        ! - -
                        q(j, i) = q(j, i) - fd - fd
                        call pot(q, vmm)
                        ! - +
                        q(l, k) = q(l, k) + fd + fd
                        call pot(q, vmp)
                        ! evaluate second derivative
                        h( nx, ny ) = &
                            (vpp - vpm + vmm - vmp) / quadru_fd_square
                        ! restore two coordinate components
                        q(j, i) = q0(j, i)
                        q(l, k) = q0(l, k)
                    else ! else using the upper elements
                        h( nx, ny ) = h( ny, nx )
                    end if
                end do
            end do
        end do
    end do
    

end subroutine calc_hess


subroutine calc_grad(q0, g)
    use constants, only: dp, fd, double_fd
    use para, only: Natom
    implicit none
    real(dp), intent(in) :: q0(3, Natom)
    real(dp), intent(out) :: g(3, Natom)
    real(dp) :: q(3, Natom), vf, vb ! forward, backward
    integer i, j

    q = q0
    do i = 1, Natom
        do j = 1, 3
            q(j,i) = q(j,i) + fd
            call pot( q, vf )
            q(j,i) = q(j,i) - fd - fd
            call pot( q, vb )
            q(j,i) = q(j,i) + fd
            g(j,i) = (vf - vb) / double_fd
        end do
    end do
    
end subroutine


subroutine cart2bond(q, r)
    use constants, only: dp
    use para, only: Natom, Nbond
    implicit none
    real(dp), intent(in) :: q(3, Natom)
    real(dp), intent(out) :: r(Nbond)
    real(dp) :: tmp(3)
    integer :: i, j, k

    k = 0
    do i = 1, Natom-1
        do j = i+1, Natom
            k = k + 1
            tmp = q(1:3,j) - q(1:3,i)
            r(k) = dsqrt( dot_product( tmp, tmp ) )
        end do
    end do
    if ( k /= Nbond ) then
        write(*, *) 'The bond amount did not equal to computed value! '
        write(*, *) k, Nbond
        stop
    end if


end subroutine