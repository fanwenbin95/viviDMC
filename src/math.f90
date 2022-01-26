! simple interface for eigenvalue decomposition
! * N  size of square matrix N
! * A  N*N double precision read matrix
! * E  N eigenvalues
! * V  N*N eigenvectors
subroutine DEEV(N, A, E, V)
    use constants, only: dp
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: A(N, N)
    real(dp), intent(out) :: E(N), V(N, N)
    real(dp) :: W(N*N*2)
    integer :: info
    real(dp) :: Vi(N, N) ! imaginary part of eigenvector
    real(dp) :: Er(N) ! right part of eigenvalue

    ! dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    call dgeev('V', 'V', & ! jobvl, jobvr, 
               N, A, N, & ! n, a, lda, 
               E, Er, & ! wr, wi, 
               V, N, & ! vl, ldvl, 
               Vi, N, & ! vr, ldvr, 
               W, N*N*2, & ! work, lwork, 
               info &
               )

    if (info /= 0) then
        write(*, *) 'ERROR! The eigen- decomposition failed. '
        write(*, *) ' * Dimension : ', N
        write(*, *) ' * info : ', info
        write(*, *) 'Check the following matrix. '
        write(*, '(6F12.6)') A
    end if

end subroutine

subroutine SVD(N, A, U, S, VT)
    use constants, only: dp
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: A(N,N)
    real(dp), intent(out) :: U(N,N), S(N), VT(N,N)
    real(dp) :: W(N*N)
    integer :: info !, lw

    call dgesvd('A', 'A', & ! jobu, jobvt
                N, N, & ! m, n
                A, N, & ! a, lda
                S, U, N, & ! s, u, ldu
                VT, N, & ! vt, ldvt
                W, N*N, & ! work, lwork
                info & ! info
                )
    ! lwork=45, for Natom=3
    
    if ( info /= 0 ) then
        write(*, '(A)') 'WARNING! The SVD results are not reliable! '
        write(*, '(A, I5)') ' * Dimension :', N
        write(*, '(A, I5)') ' * INFO flag :', info
        write(*, '(A)') 'Check the matrix below. '
        write(*, '(6F12.6)') A
    end if

end subroutine SVD


! normalize a square matrix
! * n is the size
! * axis is the dimension which will be normalized
! * A for input, Z for output
subroutine matrix_norm(n, A, axis, Z)
    use constants, only: dp
    implicit none
    integer, intent(in) :: n, axis
    real(dp), intent(in) :: A(n, n)
    real(dp), intent(out) :: Z(n, n)
    integer :: i

    if ( axis == 1 ) then
        do i = 1, n
            Z(i,:) = A(i,:) / norm2(A(i,:))
        end do
    else if ( axis == 2 ) then
        do i = 1, n
            Z(:,i) = A(:,i) / norm2(A(:,i))
        end do
    else
        write(*, *) 'ERROR! Not supported axis when normalized a matrix! '
        write(*, *) ' * axis = ', axis
    end if

end subroutine
