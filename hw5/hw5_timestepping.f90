    program hw5
    implicit none

    integer :: io=6, n_ele=1089, n_node=969, n_bounds=64, l, i, j, ii, jj, jb, hbw=50, m,it
    integer, allocatable :: ele(:, :)
    real, allocatable :: node(:, :), bc(:, :), heating_rate(:, :), a(:, :), b(:), a_ele(:, :), b_ele(:), dx(:), dy(:), a_rhs(:, :)
    real, allocatable :: a_ele_r(:, :), b_prev(:), rhs(:)
    real :: x1, x2, x3, x4, y1, y2, y3, y4, area, dx1, dx2, dx3, dy1, dy2, dy3, ke = 0, km, mm, heating_rate_m, l1, l2
    real :: material(6,3), gauss_points(4, 2), z, e, phi(4), dpx(4), dpy(4), dj, xs(4), ys(4), theta, dt, pc

    allocate(ele(n_ele, 6), node(n_node, 3), bc(n_bounds, 7), heating_rate(n_ele, 3))
    allocate(a(n_node, 2*hbw+1), b(n_node), a_ele(4, 4), b_ele(4), dx(n_ele), dy(n_ele))
    allocate(a_rhs(n_node, 2*hbw+1), a_ele_r(4, 4), b_prev(n_node), rhs(n_node))

    ele = 0
    node = 0
    bc = 0
    heating_rate = 0
    a = 0
    b = 0
    a_ele = 0
    b_ele = 0
    a_ele_r = 0
    a_rhs = 0
    b_prev = 0
    rhs = 0

    dx = 0
    dy = 0
    theta = 0.5
    dt = 0.05

    ! kappa, c* rho, m for each material
    material(1, 1) = 0.210
    material(2, 1) = 0.642
    material(3, 1) = 0.436
    material(4, 1) = 0.561
    material(5, 1) = 0.515
    material(6 ,1) = 0.642

    material(1, 2) = 2.12e6
    material(2, 2) = 3.72e6
    material(3, 2) = 2.25e6
    material(4, 2) = 3.98e6
    material(5, 2) = 4.35e6
    material(6, 2) = 3.72e6

    material(1, 3) = 200
    material(2, 3) = 2001.4
    material(3, 3) = 0
    material(4, 3) = 0
    material(5, 3) = 1482.5
    material(6, 3) = 0
    
    gauss_points(1, 1) = -0.5773502
    gauss_points(1, 2) = -0.5773502
    gauss_points(2, 1) = 0.5773502
    gauss_points(2, 2) = -0.5773502
    gauss_points(3, 1) = 0.5773502
    gauss_points(3, 2) = 0.5773502
    gauss_points(4, 1) = -0.5773502
    gauss_points(4, 2) = 0.5773502

    ! Read in mesh data
    open(newunit=io, file="epeltr4.dat", status="old", action="read")
    do i = 1, n_ele
        read(io, *) ele(i, :)
    end do
    close(io)

    open(newunit=io, file="npeltr4.dat", status="old", action="read")
    do i = 1, n_node
        read(io, *) node(i, :)
    end do
    close(io)
    
    open(newunit=io, file="bpeltr4.dat", status="old", action="read")
    do i = 1, n_bounds
        read(io, *) bc(i, :)
    end do
    close(io)

    open(newunit=io, file="ppelt4.dat", status="old", action="read")
    do i = 1, n_ele
        read(io, *) heating_rate(i, :)
    end do
    close(io)

    ! element assembly
    do l = 1, n_ele
        if (ele(l,4) == ele(l,5)) then
            ! triangular element
            mm = material(ele(l, 6)-2, 3)
            ke = material(ele(l, 6)-2, 1)
            pc = material(ele(l, 6)-2, 2)


            x1 = node(ele(l, 2), 2)
            x2 = node(ele(l, 3), 2)
            x3 = node(ele(l, 4), 2)

            y1 = node(ele(l, 2), 3)
            y2 = node(ele(l, 3), 3)
            y3 = node(ele(l, 4), 3)

            dx1 = x2 - x3
            dx2 = x3 - x1
            dx3 = x1 - x2

            dy1 = y2 - y3
            dy2 = y3 - y1
            dy3 = y1 - y2

            area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

            a_ele(1, 1) = (-ke * (dy1 * dy1)/(4 * area) -ke * (dx1 * dx1)/(4 * area) - (mm * area)/6)*dt*theta - (pc * area)/6
            a_ele(1, 2) = (-ke * (dy1 * dy2)/(4 * area) -ke * (dx1 * dx2)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(1, 3) = (-ke * (dy1 * dy3)/(4 * area) -ke * (dx1 * dx3)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(2, 1) = (-ke * (dy2 * dy1)/(4 * area) -ke * (dx2 * dx1)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(2, 2) = (-ke * (dy2 * dy2)/(4 * area) -ke * (dx2 * dx2)/(4 * area) - (mm * area)/6)*dt*theta - (pc * area)/6
            a_ele(2, 3) = (-ke * (dy2 * dy3)/(4 * area) -ke * (dx2 * dx3)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(3, 1) = (-ke * (dy3 * dy1)/(4 * area) -ke * (dx3 * dx1)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(3, 2) = (-ke * (dy3 * dy2)/(4 * area) -ke * (dx3 * dx2)/(4 * area) - (mm * area)/12)*dt*theta- (pc * area)/12
            a_ele(3, 3) = (-ke * (dy3 * dy3)/(4 * area) -ke * (dx3 * dx3)/(4 * area) - (mm * area)/6)*dt*theta- (pc * area)/6

            a_ele_r(1, 1)=(-ke * (dy1 * dy1)/(4*area)-ke * (dx1 * dx1)/(4*area) - (mm*area)/6)*dt*(1-theta) - (pc * area)/6
            a_ele_r(1, 2)=(-ke * (dy1 * dy2)/(4*area)-ke * (dx1 * dx2)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(1, 3)=(-ke * (dy1 * dy3)/(4*area)-ke * (dx1 * dx3)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(2, 1)=(-ke * (dy2 * dy1)/(4*area)-ke * (dx2 * dx1)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(2, 2)=(-ke * (dy2 * dy2)/(4*area)-ke * (dx2 * dx2)/(4*area) - (mm*area)/6)*dt*(1-theta) - (pc * area)/6
            a_ele_r(2, 3)=(-ke * (dy2 * dy3)/(4*area)-ke * (dx2 * dx3)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(3, 1)=(-ke * (dy3 * dy1)/(4*area)-ke * (dx3 * dx1)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(3, 2)=(-ke * (dy3 * dy2)/(4*area)-ke * (dx3 * dx2)/(4*area) - (mm*area)/12)*dt*(1-theta)- (pc * area)/12
            a_ele_r(3, 3)=(-ke * (dy3 * dy3)/(4*area)-ke * (dx3 * dx3)/(4*area) - (mm*area)/6)*dt*(1-theta) - (pc * area)/6

            b_ele(1) = -heating_rate(l, 3) * area * dt / 3
            b_ele(2) = -heating_rate(l, 3) * area * dt / 3
            b_ele(3) = -heating_rate(l, 3) * area * dt / 3
            
            do i = 1, 3
                ii = ele(l, i+1)
                b(ii) = b(ii) + b_ele(i)

                do j = 1, 3
                    jj = ele(l, j+1)

                    jb = (hbw+ 1) + jj - ii

                    a(ii, jb) = a(ii, jb) + a_ele(i, j)
                    a_rhs(ii, jb) = a_rhs(ii, jb) + a_ele_r(i, j)
                end do
            end do
        else 
            ! quadrilateral element
            x1 = node(ele(l, 2), 2)
            x2 = node(ele(l, 3), 2)
            x3 = node(ele(l, 4), 2)
            x4 = node(ele(l, 5), 2)
            xs(1) = x1
            xs(2) = x2
            xs(3) = x3
            xs(4) = x4

            y1 = node(ele(l, 2), 3)
            y2 = node(ele(l, 3), 3)
            y3 = node(ele(l, 4), 3)
            y4 = node(ele(l, 5), 3)

            ys(1) = y1
            ys(2) = y2
            ys(3) = y3
            ys(4) = y4

            a_ele = 0
            b_ele = 0

            do m = 1, 4
                z = gauss_points(m, 1)
                e = gauss_points(m, 2)
                phi = 0
                dpx = 0
                dpy = 0
                dj = 0
                heating_rate_m = 0
                km = 0
                mm = 0
                pc = 0
                call basis(phi, dpx, dpy, dj, z, e, xs, ys)

                ! assemble coefficients 
                do i = 1, 4
                    km = km + material(ele(l, 6)-2, 1) * phi(i)
                    mm = mm + material(ele(l, 6)-2, 3) * phi(i)
                    pc = pc + material(ele(l, 6)-2, 2) * phi(i)
                    heating_rate_m = heating_rate_m + heating_rate(l, 3) * phi(i)
                end do

                do i = 1, 4
                    do j = 1, 4
                    a_ele(i,j) = a_ele(i,j) + dj*(-km *dt*theta*(dpx(i)*dpx(j)+dpy(i)*dpy(j)) &
                        -((dt*theta*mm+pc)*phi(i)*phi(j)))
                    a_ele_r(i,j)=a_ele_r(i,j) +dj*(km*dt*(1-theta)*(dpx(i)*dpx(j)+dpy(i)*dpy(j)) &
                        +((dt*(1-theta)*mm-pc)*phi(i)*phi(j)))
                    end do
                    b_ele(i) = b_ele(i) - dj * heating_rate_m * phi(i) * dt
                end do
            end do 
            do i = 1, 4
                ii = ele(l, i+1)
                b(ii) = b(ii) + b_ele(i)
                do j = 1, 4
                    jj = ele(l, j+1)

                    jb = (hbw+ 1) + jj - ii

                    a(ii, jb) = a(ii, jb) + a_ele(i, j)
                    a_rhs(ii, jb) = a_rhs(ii, jb) + a_ele_r(i, j)
                end do
            end do
        end if 
       
    end do

     do i = 1, n_bounds
         ii = int(bc(i, 2))
         do j = 1, 2 * hbw + 1
             a(ii, j) = 0
             a_rhs(ii, j) = 0
         end do
         a(ii, hbw+1) = 1
         a_rhs(ii, hbw+1) = 1
         b(ii) = 0
     end do 
    ! apply boundary conditions 

    !do i = 1, n_bounds
    !    ii = int(bc(i, 2))
    !    l1 = sqrt((node(ii,2)-node(int(bc(i,4)),2))**2+(node(ii,3)-node(int(bc(i,4)),3))**2)
    !    l2 = sqrt((node(ii,2)-node(int(bc(i,5)),2))**2+(node(ii,3)-node(int(bc(i,5)),3))**2)
    !    b(ii) = b(ii)-0.5*l1*bc(i,6)*bc(i,7)/material(1,1) * dt
    !    b(ii) = b(ii)-0.5*l2*bc(i,6)*bc(i,7)/material(1,1) * dt
    !    a(ii, hbw+1) = a(ii, hbw+1) - (l1+l2)*bc(i,6)/material(1,1)/3 *dt * theta

    !    jb = (hbw+1) + int(bc(i,4)) - ii
    !    a(ii, jb) = a(ii, jb) - l1*bc(i,6)/material(1,1)/6 *dt * theta

    !    jb = (hbw+1) + int(bc(i,5)) - ii
    !    a(ii, jb) = a(ii, jb) - l2*bc(i,6)/material(1,1)/6 * dt * theta
    !end do

    ! solve for the unknowns
    call solve(1, a, b, n_node, hbw, n_node, 2*hbw+1)
    

    do i = 1, n_node
        rhs(i) = b(i)
    end do

    b = 0.0

    do i = 1, 3
        do j = 1, n_node
            b_prev(j) = b(j)
        end do

        b = 0.0
        do j = 1, n_node
            do it = 1, 2*hbw+1
                if (j + it - hbw - 1 > 0 .and. j + it - hbw - 1 <= n_node) then
                    b(j) = b(j) + a_rhs(j, it) * b_prev(j + it - hbw - 1) 
                end if
            end do
        end do

        do j = 1, n_node
            b(j) = b(j) + rhs(j)
        end do


        call solve(2, a, b, n_node, hbw, n_node, 2*hbw+1)
    end do        
    open(unit=7, file="output_3", status="unknown")
    do j = 1, n_node
        write(7, *) b(j)
    end do 
    close(7)

    deallocate(ele, node, a, b, a_ele, b_ele, dx, dy, a_rhs, a_ele_r, b_prev, rhs)
    end program hw5


    subroutine basis(phi, dpx, dpy, dj, z, e, xl, yl)
    
    real, intent(in) :: z, e, xl(4), yl(4)
    real, intent(out) :: phi(4), dpx(4), dpy(4), dj
    real :: dxz, dyz, dxe, dye, dpz(4), dpe(4)
    integer :: i
    
    dxz = 0
    dyz = 0
    dxe = 0
    dye = 0
    dpz = 0
    dpe = 0
    phi = 0
    dpx = 0
    dpy = 0
    dj = 0

    phi(1) = (1-z)*(1-e)/4
    phi(2) = (1+z)*(1-e)/4
    phi(3) = (1+z)*(1+e)/4
    phi(4) = (1-z)*(1+e)/4

    dpz(1) = -(1-e)/4
    dpz(2) = (1-e)/4
    dpz(3) = (1+e)/4
    dpz(4) = -(1+e)/4

    dpe(1) = -(1-z)/4
    dpe(2) = -(1+z)/4
    dpe(3) = (1+z)/4
    dpe(4) = (1-z)/4


    do i = 1, 4
        dxz = dxz + xl(i) * dpz(i)
        dxe = dxe + xl(i) * dpe(i)
        dyz = dyz + yl(i) * dpz(i)
        dye = dye + yl(i) * dpe(i)
    end do

    dj = dxz * dye - dyz * dxe

    do i = 1, 4
        dpx(i) = (dye * dpz(i) - dyz * dpe(i))/dj
        dpy(i) = (-dxe * dpz(i) + dxz * dpe(i))/dj
    end do

    end subroutine basis