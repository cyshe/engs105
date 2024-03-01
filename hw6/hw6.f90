    program hw6
    implicit none

    integer :: io=6, n_ele = 38, n_node=28, l, i, j, ii, jj, jb, hbw=6
    integer, allocatable :: ele(:, :)
    real, allocatable :: node(:, :), a(:, :), b(:), a_ele(:, :), b_ele(:), dx(:), dy(:)
    real :: x1, x2, x3, y1, y2, y3, area, dx1, dx2, dx3, dy1, dy2, dy3, ke = 0

    allocate(ele(n_ele, 5), node(n_node, 3), a(n_node, n_node), b(n_node), a_ele(3, 3), b_ele(3), dx(n_ele), dy(n_ele))
    
    ele = 0
    node = 0
    a = 0
    b = 0
    a_ele = 0
    b_ele = 0
    dx = 0
    dy = 0

    ! Read in mesh data
    open(newunit=io, file="hw44.ele", status="old", action="read")
    do i = 1, n_ele
        read(io, *) ele(i, :)
    end do
    close(io)

    open(newunit=io, file="hw44.nod", status="old", action="read")
    do i = 1, n_node
        read(io, *) node(i, :)
    end do
    close(io)

   
    ! element assembly
    do l = 1, n_ele
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
    end do

    ! apply boundary conditions
    do i = 1, n_node
        if (i == 1 .or. i == 6 .or. i == 11 .or. i == 16) then
            do j = 1, hbw*2 +1
                a(i, j) = 0
            end do
            a(i, hbw + 1) = 1
            b(i) = 0

        else if (i == 5 .or. i == 10 .or. i == 15 .or. i == 20 .or. i == 24 .or. i == 27.or. i == 28 .or. i == 25) then
            do j = 1, hbw*2 +1  
                a(i, j) = 0
            end do
            a(i, hbw + 1) = 1
            b(i) = 1
        
        else if (i < 5) then
            do j = 1, hbw*2 +1  
                a(i, j) = 0
            end do
            a(i, hbw + 1) = 1
            b(i) = (i-1)/4.0
        !else if ( ) then
        !    b(i) = 0
        end if 
    end do

    ! solve for the unknowns
    call solve(1, a, b, n_node, hbw, n_node, 2*hbw+1)
    call solve(2, a, b, n_node, hbw, n_node, 2*hbw+1)

    ! output solution
    open(unit=7, file="output_1", status="unknown")
    do i = 1, n_node
        write(7, *) b(i)
    end do 
    close(7)


    do l = 1, n_ele
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
        
        dx(l) = b(ele(l, 2)) * dy1 + b(ele(l, 3)) * dy2 + b(ele(l, 4)) * dy3
        dy(l) = -b(ele(l, 2)) * dx1 - b(ele(l, 3)) * dx2 - b(ele(l, 4)) * dx3
    end do


    ! output derivatives
    open(unit=8, file="flow_x_1.dat", status="unknown")
    do i = 1, n_ele
        write(8, *) dy(i)
    end do
    close(8)

    open(unit=9, file="flow_y_1.dat", status="unknown")
    do i = 1, n_ele
        write(9, *) -dx(i)
    end do
    close(9)


    deallocate(ele, node, a, b, a_ele, b_ele, dx, dy)
    end program hw6