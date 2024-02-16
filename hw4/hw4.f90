    program hw4
    implicit none

    integer :: io, n_ele = 38, n_node=28, l, i, j, ii, jj, jb
    real :: 
    real, allocatable :: ele(:, :), node(:, :), a(:, :), b(:), a_ele(:, :), b_ele(:)

    allocate(ele(n_ele, 5), node(n_node, 3), a_ele(3, 3), b_ele(3))
    
    ! Read in mesh data
    open(newunit=io, file="hw4.ele", status="old", action="read")
    read(io, *) ele
    close(io)

    open(newunit=io, file="hw4.node", status="old", action="read")
    read(io, *) node
    close(io)

    


    ! element assembly
    do l = 1, n_ele
        do i = 1, 3
            ii = ele(l, i+1)
            b(ii) = b(ii) + b_ele(i)

            do j = 1, 3
                jj = ele(l, j+1)

                jb = () + jj - ii

                a(ii, jb) = a(ii, jb) + a_ele(i, j)
            end do
        end do
    end do

    ! apply boundary conditions

    ! solve for the unknowns

    ! output solution


    end program hw4

    subroutine element_matrix(l, ele, node, a_ele, b_ele)
        implicit none
        integer, intent(in) :: l
        real, intent(in) :: ele(5, 38), node(3, 28)
        real, intent(out) :: a_ele(3, 3), b_ele(3)

        ! local variables
        real :: x1, x2, x3, y1, y2, y3, area, k
        
        x1 = node(ele(l, 2), 2)
        x2 = node(ele(l, 3), 2)
        x3 = node(ele(l, 4), 2)

        y1 = node(ele(l, 2), 3)
        y2 = node(ele(l, 3), 3)
        y3 = node(ele(l, 4), 3)

        area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        
        

        
    end subroutine element_matrix