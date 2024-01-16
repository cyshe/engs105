    program test_solver
    implicit none


    integer :: i, j, n, io
    real, allocatable :: a_mat(:, :), x(:), b(:)
    
    
    n = 20
    
    allocate(a_mat(n, 3), b(n))

    do i = 1, n
        a_mat(i, 1) = 2.0
        a_mat(i, 2) = 2.0
        a_mat(i, 3) = 0.0
        b(i) = i
    end do

    call solve(3, a_mat, b, n, 1, n, 3)

    do i = 1, n
        print *, b(i)
    end do


    end program