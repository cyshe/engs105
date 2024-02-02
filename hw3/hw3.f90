    program hw3
    implicit none

    integer :: i, j, n, idx, io
    real::h, k, a, R, beta, rij, pi, a_const, theta, r_mult
    real, allocatable :: a_mat(:, :), b(:), rhs_mat(:, :)
    
    pi=4.D0*DATAN(1.D0)

    n = 80
    a = 0.1
    R = 1.0
    a_const = 30.0
    theta = 0.7
    r_mult = 1.0        !r = Dk/h^2 in FD equation

    h = (R-a)/(n-1)     !dr
    k = 0.5*pi/(n-1)    !dtheta
    beta = (h/k) * (h/k)

    print *, "pi = ", pi
    
    allocate(a_mat(n*n, 2*n+1), b(n*n), rhs_mat(n*n, 2*n+1))
    a_mat = 0.0
    b = 0.0
    rhs_mat = 0.0
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
                a_mat(idx, n+2) = 
                
                a_mat(idx, n+1) = 
                
                a_mat(idx, 2*n+1) = 
                

            else if (i == 1 .and. j == n) then
                a_mat(idx, n+2) =   

                a_mat(idx, n+1) = 
                
                a_mat(idx, 1) = 
                
                
            else if (i == n .and. j == 1) then
                a_mat(idx, n+1) = 1


            else if (i == n .and. j == n) then
                a_mat(idx, n+1) = 1


    
            else if (i == 1) then
                a_mat(idx, n+2) = 2.0 * rij * rij  

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = beta
                a_mat(idx, 1) = beta


            else if (i == n) then
                a_mat(idx, n+1) = 1
                

            else if (j == 1) then
                a_mat(idx, n) = -h * rij/2 + rij * rij
                a_mat(idx, n+2) = h * rij/2 + rij * rij  

                a_mat(idx, n+1) = -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = 2.0 * beta

            else if (j == n) then
                a_mat(idx, n+1) = 1

                rhs_mat(idx, n+1) = 1

            else
                a_mat(idx, n) = -r_mult * theta
                a_mat(idx, n+2) = -r_mult * theta  

                a_mat(idx, n+1) = 1 + 4 * r_mult * theta
                
                a_mat(idx, 2*n+1) = r_mult * theta
                a_mat(idx, 1) = r_mult * theta


                rhs_mat(idx, n) = (1 -theta) * r_mult
                rhs_mat(idx, n+2) = (1 - theta) * r_mult

                rhs_mat(idx, n+1) = 1 + 4 * (1-theta) * r_mult

                rhs_mat(idx, 2*n+1) = (1 - theta) * r_mult
                rhs_mat(idx, 1) = (1 - theta) * r_mult
            end if
        end do
    end do
    
    

    print *, "Solve matrix equation"
    call solve(3, a_mat, b, n*n, n, n*n, 2*n+1)

    print *, "Write solution to file"    

    open(unit=io, file="200_a30.dat", status="unknown")
    do i = 1, n*n
        write(io, *) b(i)
    end do
    
    deallocate(a_mat, b)
    end program hw3