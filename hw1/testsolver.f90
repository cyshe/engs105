    program test_solver
    implicit none


    integer :: i, n
    real :: sum
    real, allocatable :: a_mat(:, :), x(:), b(:)
    
    
    n = 5
    
    allocate(a_mat(n, 3), b(n), x(n))
    a_mat = 0.0
    b = 0.0
    x = 0.0

    do i = 1, n
        
        a_mat(i, 1) = 1.0
        a_mat(i, 2) = 1.0
        a_mat(i, 3) = 1.0
        b(i) = 3
        x(i) = 3

        if (i == 1) then
            a_mat(i, 1) = 0.0
            b(i) = 2
            x(i) = 2
        end if
        if (i == n) then
            a_mat(i, 3) = 0.0
            b(i) = 2
            x(i) = 2
        end if
        

        
    end do
    do i = 1, n
        print *, b(i), a_mat(i, 1), a_mat(i, 2), a_mat(i, 3) 
    end do
    print *, "-------------------"
    call solve(3, a_mat, b, n, 1, n, 3)

    do i = 1, n
        print *, b(i), a_mat(i, 1), a_mat(i, 2), a_mat(i, 3) 
    end do

    do i = 1, n
        sum = 0
        if (i == 1) then 
            sum = sum + a_mat(i, 2) * b(i) + a_mat(i, 3) * b(i + 1)
        else if (i == n) then
            sum = sum + a_mat(i, 1) * b(i - 1) + a_mat(i, 2) * b(i)
        else
            sum = sum + a_mat(i, 1) * b(i - 1) + a_mat(i, 2) * b(i) + a_mat(i, 3) * b(i + 1)
        end if
        print *, sum - x(i)
        print *, b(i), a_mat(i, 1), a_mat(i, 2), a_mat(i, 3) 
    end do

    deallocate(a_mat, b, x)

    end program