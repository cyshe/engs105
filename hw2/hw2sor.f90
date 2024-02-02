program hw2_sor
    implicit none
    integer :: i, j, n, idx, io, iters
    real*8 :: R = 1.0, pi, h, k, a, a_const, rij, beta, delta, delta_prev, omega, rhoj
    real*8, allocatable :: b(:), b_prev(:)
    
    n = 40
    allocate(b(n*n), b_prev(n*n))
    pi = 4.0 * atan(1.0)

    a_const = 10.0
    a = 0.1
    delta_prev = 0.0
    delta = 1.0
    b_prev = 0.0
    b = 10.0
    h = 0.9/(n-1)
    k = 0.5 * pi / (n-1)
    beta = (h/k) * (h/k)
    rhoj =  0.999
    omega = 2.0/ (1 + (1-rhoj*rhoj)**0.5)
    iters = 0
    
    do while (maxval(abs(b-b_prev)) > 1e-5)
        delta_prev = delta
        delta = NORM2(b_prev - b)
        print *, delta/delta_prev
        b_prev = b
        iters = iters + 1

        do i = 1, n
            do j = 1, n
                idx = (j-1)*n + i
                rij = a + (i-1)*h

                if (i == n) then
                    b(idx) = - R * cos(3 * k * (j-1))
                
                else if (i == 1 .and. j == 1) then
                    b(idx) = (2.0 * rij * rij) * b_prev(idx+1) & 
                            + 2.0 * beta * b_prev(idx+n)

                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta))+ (1-omega) * b_prev(idx)

                else if (i == 1 .and. j == n) then
                    b(idx) = (2.0 * rij * rij) * b_prev(idx+1) & 
                            + 2.0 * beta * b(idx-n) &
                            - 2.0 * 9.0 * beta

                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta + 3*beta*a_const))+ (1-omega) * b_prev(idx)

                else if (i == 1) then
                    b(idx) = (2.0 * rij * rij) * b_prev(idx+1) &
                            + beta * b(idx-n) &
                            + beta * b_prev(idx+n)

                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta))+ (1-omega) * b_prev(idx)

                else if (j == 1) then
                    b(idx) = (-h*rij/2.0 + rij * rij)* b(idx-1) &
                            +(h* rij/2.0 + rij * rij) * b_prev(idx+1) &
                            + 2.0 * beta * b_prev(idx+n)

                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta))+ (1-omega) * b_prev(idx)
                
                else if (j == n) then
                    b(idx) = (-h * rij/2.0 + rij * rij)* b(idx-1) &
                            +(h * rij/2.0 + rij * rij) * b_prev(idx+1) &
                            + 2.0 * beta * b(idx-n) &
                            - 2 * 9 * beta

                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta + 3 * beta * a_const)) + (1-omega) * b_prev(idx)

                else
                    b(idx) = (-h * rij/2.0 + rij * rij)* b(idx-1) &
                            +(h * rij/2.0 + rij * rij) * b_prev(idx+1) &
                            + beta * b(idx-n) &
                            + beta * b_prev(idx+n)
                    
                    b(idx) = b(idx) * omega / (2.0 * (rij * rij + beta))+ (1-omega) * b_prev(idx)
                end if
            end do
        end do
    end do

    print *, "iteration: ", iters
    print *, "omega: ", omega
    open(unit=io, file="40sor_a10.dat", status="unknown")
    do i = 1, n*n
        write(io, *) b(i)
    end do

    end program hw2_sor