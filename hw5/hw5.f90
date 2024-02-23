    program hw5
    implicit none

    integer :: io=6, n_ele=1089, n_node=969, n_bounds=64, l, i, j, ii, jj, jb, hbw=6
    integer, allocatable :: ele(:, :)
    real, allocatable :: node(:, :), bc(:, :), heating_rate(:, :), a(:, :), b(:), a_ele(:, :), b_ele(:), dx(:), dy(:)
    real :: x1, x2, x3, x4, y1, y2, y3, y4, area, dx1, dx2, dx3, dy1, dy2, dy3, ke = 0, material(6), gauss_points(4, 2)

    allocate(ele(n_ele, 6), node(n_node, 3), bc(n_bounds, 7), heating_rate(n_ele, 3))
    allocate(a(n_node, n_node), b(n_node), a_ele(4, 4), b_ele(4), dx(n_ele), dy(n_ele))

    ele = 0
    node = 0
    a = 0
    b = 0
    a_ele = 0
    b_ele = 0
    dx = 0
    dy = 0

    material(1) = 200
    material(2) = 2001.4
    material(3) = 0
    material(4) = 0
    material(5) = 1482.5
    material(6) = 0
    
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
    do i = 1, n_node
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
        if (ele(l, 4) == ele(l,5)) then
            ! triangular element
            ke = material(ele(l, 6))

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

            a_ele(1, 1) = - (dy1 * dy1)/(4 * area) - (dx1 * dx1)/(4 * area) + (ke * ke * area)/6
            a_ele(1, 2) = - (dy1 * dy2)/(4 * area) - (dx1 * dx2)/(4 * area) + (ke * ke * area)/12
            a_ele(1, 3) = - (dy1 * dy3)/(4 * area) - (dx1 * dx3)/(4 * area) + (ke * ke * area)/12


            a_ele(2, 1) = - (dy2 * dy1)/(4 * area) - (dx2 * dx1)/(4 * area) + (ke * ke * area)/12
            a_ele(2, 2) = - (dy2 * dy2)/(4 * area) - (dx2 * dx2)/(4 * area) + (ke * ke * area)/6
            a_ele(2, 3) = - (dy2 * dy3)/(4 * area) - (dx2 * dx3)/(4 * area) + (ke * ke * area)/12

            a_ele(3, 1) = - (dy3 * dy1)/(4 * area) - (dx3 * dx1)/(4 * area) + (ke * ke * area)/12
            a_ele(3, 2) = - (dy3 * dy2)/(4 * area) - (dx3 * dx2)/(4 * area) + (ke * ke * area)/12
            a_ele(3, 3) = - (dy3 * dy3)/(4 * area) - (dx3 * dx3)/(4 * area) + (ke * ke * area)/6

            b_ele(1) = 0
            b_ele(2) = 0
            b_ele(3) = 0

            do i = 1, 3
                ii = ele(l, i+1)
                b(ii) = b(ii) + b_ele(i)

                do j = 1, 3
                    jj = ele(l, j+1)

                    jb = (hbw+ 1) + jj - ii

                    a(ii, jb) = a(ii, jb) + a_ele(i, j)
                end do
            end do
        else 
            ! quadrilateral element
            ke = material(ele(l, 6))

            x1 = node(ele(l, 2), 2)
            x2 = node(ele(l, 3), 2)
            x3 = node(ele(l, 4), 2)
            x4 = node(ele(l, 5), 2)

            y1 = node(ele(l, 2), 3)
            y2 = node(ele(l, 3), 3)
            y3 = node(ele(l, 4), 3)
            y4 = node(ele(l, 5), 3)

            do m = 1, 4
                z = gauss_points(m, 1)
                e = gauss_points(m, 2)

                call basis(phi, dpx, dpy, dj, z, e, [x1, x2, x3, x4], [y1, y2, y3, y4])

                do i = 1, 4
                    do j = 1, 4
                        
                    end do      
                end do

                


            end do 
        end if 
       
    end do

    ! apply boundary conditions
    do i = 1, n_node
        if then
        end if 
    end do

    ! solve for the unknowns
    call solve(1, a, b, n_node, hbw, n_node, 2*hbw+1)
    call solve(2, a, b, n_node, hbw, n_node, 2*hbw+1)

    ! output solution
    open(unit=7, file="output_2", status="unknown")
    do i = 1, n_node
        write(7, *) b(i)
    end do 
    close(7)





    deallocate(ele, node, a, b, a_ele, b_ele, dx, dy)
    end program hw5


    subroutine basis(phi, dpx, dpy, dj, z, e, xl, yl)
    
    real, intent(in) :: z, e, xl, yl
    real, intent(out) :: phi(4), dpx(4), dpy(4), dj(4, 4)
    real :: x, y, dxz, dxe, dye, dx, dy

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
