    program hw3
    implicit none

    integer :: i, j, n, idx, io, it
    real::h, k, a, R, beta, rij, pi, a_const, theta, r_mult
    real, allocatable :: a_mat(:, :), b(:), rhs_mat(:, :), b_prev(:)
    
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
    
    allocate(a_mat(n*n, 2*n+1), b(n*n), rhs_mat(n*n, 2*n+1), b_prev(n*n))
    a_mat = 0.0
    b = 0.0
    rhs_mat = 0.0
    b_prev = 0.0
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
                a_mat(idx, n+2) = -r_mult * theta * 2.0 * rij * rij

                a_mat(idx, n+1) = 1 - r_mult * theta * -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = -r_mult * theta * 2.0 *  beta


                rhs_mat(idx, n+2) = (1 - theta)  * r_mult * 2.0 * rij * rij

                rhs_mat(idx, n+1) = 1 + (1-theta) * r_mult * -2.0 * (rij * rij + beta)

                rhs_mat(idx, 2*n+1) = (1 - theta) * r_mult * 2 * beta
                

            else if (i == 1 .and. j == n) then
                a_mat(idx, n+1) = 1 

                rhs_mat(idx, n+1) = 1
                
                
            else if (i == n .and. j == 1) then
                a_mat(idx, n+1) = 1

                rhs_mat(idx, n+1) = 1

                b(idx) = -R

            else if (i == n .and. j == n) then
                a_mat(idx, n+1) = 1

                rhs_mat(idx, n+1) = 1



    
            else if (i == 1) then
                a_mat(idx, n+2) = -r_mult * theta * 2 * rij * rij 

                a_mat(idx, n+1) = 1 - r_mult * theta * -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = -r_mult * theta * beta
                a_mat(idx, 1) = -r_mult * theta * beta


                rhs_mat(idx, n+2) = (1 - theta) * r_mult * 2.0 * rij * rij

                rhs_mat(idx, n+1) = 1 + (1-theta) * r_mult * -2.0 * (rij * rij + beta)

                rhs_mat(idx, 2*n+1) = (1 - theta) * r_mult * beta
                rhs_mat(idx, 1) = (1 - theta) * r_mult * beta


            else if (i == n) then
                a_mat(idx, n+1) = 1

                rhs_mat(idx, n+1) = 1
                b(idx) = -R * cos(3 * k * (j-1))

            else if (j == 1) then
                a_mat(idx, n) = -r_mult * theta * (-h * rij/2.0 + rij * rij)
                a_mat(idx, n+2) = -r_mult * theta * (h * rij/2.0 + rij * rij) 

                a_mat(idx, n+1) = 1 - r_mult * theta * -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = -r_mult * theta * 2 * beta


                rhs_mat(idx, n) = (1 -theta) * r_mult * (-h * rij/2.0 + rij * rij)
                rhs_mat(idx, n+2) = (1 - theta) * r_mult * (h * rij/2.0 + rij * rij)

                rhs_mat(idx, n+1) = 1 + (1-theta) * r_mult * -2.0 * (rij * rij + beta)

                rhs_mat(idx, 2*n+1) = (1 - theta) * r_mult * 2 * beta

            else if (j == n) then
                a_mat(idx, n+1) = 1

                rhs_mat(idx, n+1) = 1

            else
                a_mat(idx, n) = -r_mult * theta * (-h * rij/2.0 + rij * rij)
                a_mat(idx, n+2) = -r_mult * theta * (h * rij/2.0 + rij * rij)

                a_mat(idx, n+1) = 1 - r_mult * theta * -2.0 * (rij * rij + beta)
                
                a_mat(idx, 2*n+1) = -r_mult * theta * beta
                a_mat(idx, 1) = -r_mult * theta * beta


                rhs_mat(idx, n) = (1 -theta) * r_mult * (-h * rij/2.0 + rij * rij)
                rhs_mat(idx, n+2) = (1 - theta) * r_mult * (h * rij/2.0 + rij * rij)

                rhs_mat(idx, n+1) = 1 + (1-theta) * r_mult * -2.0 * (rij * rij + beta)

                rhs_mat(idx, 2*n+1) = (1 - theta) * r_mult * beta
                rhs_mat(idx, 1) = (1 - theta) * r_mult * beta
            end if
        end do
    end do 
    

    print *, "Solve matrix equation"
    
    do i = 1, 100
        do j = 1, n*n
            b_prev(j) = b(j)
        end do
        b = 0.0
        do j = 1, n*n
            do it = 1, 2*n+1
                if (j + it - n - 1 > 0 .and. j + it - n - 1 <= n*n) then
                    b(j) = b(j) + rhs_mat(j, it) * b_prev(j + it - n - 1)
                end if
            end do
        end do
        call solve(3, a_mat, b, n*n, n, n*n, 2*n+1)
    end do

    print *, "Write solution to file"    

    open(unit=io, file="output.dat", status="unknown")
    do i = 1, n*n
        write(io, *) b(i)
    end do
    
    deallocate(a_mat, b)
    end program hw3