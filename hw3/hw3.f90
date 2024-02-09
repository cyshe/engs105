    program hw3
    implicit none

    integer :: i, j, n, idx, io = 6, it, timesteps = 5000
    real::h, k, a, R, beta, rij, pi, a_const, theta, r_mult, temp
    real, allocatable :: a_mat(:, :), b(:), rhs_mat(:, :), b_prev(:), p1(:), p2(:), p3(:), p4(:), p5(:)
    
    pi=4.D0*DATAN(1.D0)

    n = 80
    a = 0.1
    R = 1.0
    a_const = 30.0
    theta = 0.0
    r_mult = 0.25        !r = Dk/h^2 in FD equation

    h = (R-a)/(n-1)     !dr
    k = 0.5*pi/(n-1)    !dtheta
    beta = (h/k) * (h/k)

    print *, "pi = ", pi
    
    allocate(a_mat(n*n, 2*n+1), b(n*n), rhs_mat(n*n, 2*n+1), b_prev(n*n))
    allocate(p1(timesteps), p2(timesteps), p3(timesteps), p4(timesteps), p5(timesteps))
    a_mat = 0.0
    b = 0.0
    rhs_mat = 0.0
    b_prev = 0.0

    p1 = 0.0
    p2 = 0.0
    p3 = 0.0
    p4 = 0.0
    p5 = 0.0

    ! assemble matrix .and. rhs
    ! i in direction r, j in direction theta
    print *, "Assemble matrix equation"
    do i = 1, n
        ! print *, "i = ", i
        do j = 1, n
            idx = n*(j-1) + i
            rij = a + (i-1)*h
            ! print *, "idx = ", idx
            if (i == 1 .and. j == 1) then
                a_mat(idx, n+2) = 2.0 * rij * rij

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = 2.0 *  beta

            else if (i == 1 .and. j == n) then
                a_mat(idx, n+1) = 1 
                
            else if (i == n .and. j == 1) then
                a_mat(idx, n+1) = 1
                b(idx) = -R

            else if (i == n .and. j == n) then
                a_mat(idx, n+1) = 1

    
            else if (i == 1) then
                a_mat(idx, n+2) = 2 * rij * rij 

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = beta
                a_mat(idx, 1) = beta

            else if (i == n) then
                a_mat(idx, n+1) = 1

                b(idx) = -R * cos(3 * k * (j-1))

            else if (j == 1) then
                a_mat(idx, n) = -h * rij/2.0 + rij * rij
                a_mat(idx, n+2) = h * rij/2.0 + rij * rij 

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = 2.0 * beta

            else if (j == n) then
                a_mat(idx, n+1) = 1

            else
                a_mat(idx, n) = -h * rij/2.0 + rij * rij
                a_mat(idx, n+2) = h * rij/2.0 + rij * rij

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = beta
                a_mat(idx, 1) = beta

            end if
        end do
    end do
    

    do i = 1, n*n
        do j = 1, 2*n+1
            temp = a_mat(i, j)
            rhs_mat(i, j) = temp * r_mult * (1.0 - theta)
            if (j == n+1) then
                rhs_mat(i, j) = 1 + temp * r_mult * (1.0 - theta)
            end if

            a_mat(i, j) = temp * (-1.0) * r_mult * theta
            if (j == n+1) then
                a_mat(i, j) = 1 + (temp * (-1.0) * r_mult * theta)
            end if
        end do
    end do
    
    do i = 1, n
        do j = 1, n
            idx = n*(j-1) + i
            if (i == n .or. j == n) then
                a_mat(idx, n+1) = 1 
                rhs_mat(idx, n+1) = 1
            end if
        end do
    end do

    

    call solve(1, a_mat, b, n*n, n, n*n, 2*n+1)
    print *, "Solve matrix equation"
    
    do i = 1, timesteps
        do j = 1, n*n
            b_prev(j) = b(j)
        end do
        p1(i) = b(1245)
        p2(i) = b(2245)
        p3(i) = b(3245)
        p4(i) = b(4245)
        p5(i) = b(5245)

        b = 0.0

        do j = 1, n*n
            do it = 1, 2*n+1
                if (j + it - n - 1 > 0 .and. j + it - n - 1 <= n*n) then
                    b(j) = b(j) + rhs_mat(j, it) * b_prev(j + it - n - 1)
                end if
            end do
        end do
        call solve(2, a_mat, b, n*n, n, n*n, 2*n+1)
    end do

    print *, "Write solution to file"    

    open(unit=io, file="output.dat", status="unknown")
    do i = 1, n*n
        write(io, *) b(i)
    end do
    close(io)

    open(unit=io, file="p1.dat", status="unknown")
    do i = 1, timesteps
        write(io, *) p1(i)
    end do
    close(io)

    open(unit=io, file="p2.dat", status="unknown")
    do i = 1, timesteps
        write(io, *) p2(i)
    end do
    close(io)

    open(unit=io, file="p3.dat", status="unknown")
    do i = 1, timesteps
        write(io, *) p3(i)
    end do
    close(io)

    open(unit=io, file="p4.dat", status="unknown")
    do i = 1, timesteps
        write(io, *) p4(i)
    end do
    close(io)

    open(unit=io, file="p5.dat", status="unknown")
    do i = 1, timesteps
        write(io, *) p5(i)
    end do
    close(io)
    
    deallocate(a_mat, b, rhs_mat, b_prev, p1, p2, p3, p4, p5)
    end program hw3