    program method_b
    implicit none

    integer :: i, j, n, io
    real :: L, D, h, dt, r, sigma
    real, allocatable :: a_mat(:, :), b_mat1(:, :), b_mat2(:, :), b(:), b_prev(:), rhs(:)

    L = 10.0
    D = 0.5
    h = 0.1
    dt = 0.01
    r = D*dt/h**2
    sigma = 0.1

    n = int(L/h) + 1

    ! form matrix systems
    allocate(a_mat(n, 3), b_mat1(n, 3), b_mat2(n, 3))
    allocate(b(n), b_prev(n), rhs(n))

    a_mat = 0.0
    b_mat1 = 0.0
    b_mat2 = 0.0
    b = 0.0
    b_prev = 0.0

    do i = 1, n
        a_mat(i, 2) = 1.0
        
        if (i == 1) then
            b_mat1(i, 1) = 0.0
            b_mat1(i, 2) = 0.0
            b_mat1(i, 3) = 0.0

            b_mat2(i, 1) = 0.0
            b_mat2(i, 2) = 0.0
            b_mat2(i, 3) = 0.0

        else if (i == n) then
            b_mat1(i, 1) = 0.0
            b_mat1(i, 2) = 0.0
            b_mat1(i, 3) = 0.0

            b_mat2(i, 1) = 0.0
            b_mat2(i, 2) = 0.0
            b_mat2(i, 3) = 0.0
        else
            a_mat(i, 1) = -5 * r / 12
            a_mat(i, 2) = 1 + 5 * r / 6
            a_mat(i, 3) = -5 * r / 12

            b_mat1(i, 1) = 2 * r / 3
            b_mat1(i, 2) = 1 - 4 * r/3
            b_mat1(i, 3) = 2 * r / 3

            b_mat2(i, 1) = -r / 12
            b_mat2(i, 2) = r / 6
            b_mat2(i, 3) = -r / 12
        end if
    end do

    call solve(1, a_mat, b, n, 1, n, 3)
    print *, "Solve matrix equation"

    ! initialize b and b_prev
    do i = 1, n
        b(i) = exp(-((i-1) * h - 5)**2/ (2.0 * sigma **2))
        b_prev(i) = exp(-((i-1) * h - 5)**2/ (2.0 * sigma **2))
    end do




    do j = 1, 150
        ! calculate rhs vector from previous time steps
        do i = 1, n
            if (i == 1) then
                rhs(i) = b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) &
                           + b_mat1(i, 3) * b(i+1) + b_mat2(i, 3) * b_prev(i+1)          

            else if (i == n) then
                rhs(i) = b_mat1(i, 1) * b(i-1) + b_mat2(i, 1) * b_prev(i-1) &
                           + b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) 
            else
                 rhs(i) = b_mat1(i, 1) * b(i-1) + b_mat2(i, 1) * b_prev(i-1) &
                           + b_mat1(i, 2) * b(i) + b_mat2(i, 2) * b_prev(i) &
                           + b_mat1(i, 3) * b(i+1) + b_mat2(i, 3) * b_prev(i+1)
            end if
        end do
        
        call solve(2, a_mat, rhs, n, 1, n, 3)

        ! update b_prev and b
        b_prev = b
        b = rhs

        if (j == 1) then
            open(unit=io, file="b_t1.dat", status="unknown")
            do i = 1, n
                write(io, *) b(i)
            end do
            close(io)
        end if
        
        if (j == 4) then
            open(unit=io, file="b_t4.dat", status="unknown")
            do i = 1, n
                write(io, *) b(i)
            end do
            close(io)
        end if 

        if (j == 15) then
            open(unit=io, file="b_t15.dat", status="unknown")
            do i = 1, n
                write(io, *) b(i)
            end do
            close(io)
        end if


        if (j == 150) then
            open(unit=io, file="b_t30.dat", status="unknown")
            do i = 1, n
                write(io, *) b(i)
            end do
            close(io)
        end if

    end do

    end program method_b