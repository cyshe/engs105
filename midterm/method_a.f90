    program method_a
    implicit none

    integer :: i
    real :: L, D, h, dt, r
    real, allocatable :: a_mat(:, :), b_mat1(:, :), b_mat2(:, :), b(:), b_prev(:)

    L = 10.0
    D = 0.5
    h = 0.1
    dt = 0.05
    r = D*dt/h**2

    n = L/h

    ! form matrix systems
    allocate(a_mat(n, 3), b_mat1(n, 3), b_mat2(n, 3))
    allocate(b(n), b_prev(n))

    a_mat = 0.0
    b_mat1 = 0.0
    b_mat2 = 0.0
    b = 0.0
    b_prev = 0.0

    do i = 1, n
        a_mat(i, 2) = 1.0
        
        if (i == 1) then
            b_mat1(i, 1) = 0.0
            b_mat1(i, 2) = 1 - r
            b_mat1(i, 3) = 0.0

            b_mat2(i, 1) = 0.0
            b_mat2(i, 2) = r
            b_mat2(i, 3) = 0.0

        else if (i == n) then
            b_mat1(i, 1) = 0.0
            b_mat1(i, 2) = 1 - r
            b_mat1(i, 3) = 0.0

            b_mat2(i, 1) = 0.0
            b_mat2(i, 2) = r
            b_mat2(i, 3) = 0.0
        else
            b_mat1(i, 1) = 3 * r / 2
            b_mat1(i, 2) = 1 - r
            b_mat1(i, 3) = 3 * r / 2

            b_mat2(i, 1) = -r / 2
            b_mat2(i, 2) = r
            b_mat2(i, 3) = -r / 2
        end if
    end do

    call solve(1, a_mat, b, n*n, n, n*n, 2*n+1)
    ! form right hand side vector b1 * u^l + b2 * u^l-1

    do i = 1, n
        if (i == 1) then
            b(i) = b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) &
                   + b_mat1(i, 3) * b(i+1) + b_mat2(i, 3) * b_prev(i+1)          
        
        else if (i == n) then
            b(i) = b_mat1(i, 1) * b(i-1) + b_mat2(i, 1) * b_prev(i-1) &
                   + b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) &
        else
            b(i) = b_mat1(i, 1) * b(i-1) + b_mat2(i, 1) * b_prev(i-1) &
                   + b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) &
                   + b_mat1(i, 3) * b(i+1) + b_mat2(i, 3) * b_prev(i+1)
        end if
    end do


    print *, "Solve matrix equation"

    call solve(2, a_mat, b, n*n, n, n*n, 2*n+1)
    

    end program method_a