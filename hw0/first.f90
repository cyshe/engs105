    program hw0
    implicit none

    integer :: i, j, io, k
    real :: sum
    real, dimension(1200, 6) :: A
    real, dimension(6, 6) :: B, symB, symB_norm
    real, dimension(1200, 1200) :: C, symC, symC_norm

    open(unit=io, file='hw0.dat', status='old', action='read')
    read(io, *)
    read(io, *)
    do i = 1, 1200
        read(io, *) A(i, :)
    end do
    close(io)
    
    print*, A(1, :) ! print the first row of A
    print*, A(1200, :) ! print the last row of A
    

    B(:,:) = 0.0
    
    do i = 1, 6
        do j = 1, 6
            sum = 0.0
            do k = 1, 1200
                sum = sum + A(k, i) * A(k, j)
            end do
            B(i, j) = sum
        end do
    end do
    print*, "B: "
    call print_matrix(B, 6, 6) ! print the matrix

    open(newunit=io, file="B.dat", status="new", action="write")
    do i = 1, 6
        write(io, *) B(i, :)
    end do
    close(io)

    ! extract the diagonal elements of B
    open (newunit=io, file="B_diag.dat", status="new", action="write")
    do i = 1, 6
        write(io, *) B(i, i)
    end do
    close(io)


    ! Symmetry check (not normalized)
    do i = 1, 6
        do j = 1, 6
            symB(i, j) = abs(B(i, j) - B(j,i))
        end do
    end do 
    
    print*, "Symmetry Check (not normalized): "
    call report_symmetry(symB, 6)

    ! symmetry check (normalized)
    do i = 1, 6
        do j = 1, 6
            symB_norm(i, j) = abs((B(i, j) - B(j,i)) / ((B(i, j) + B(j, i))/2))
        end do
    end do
    print *, "Symmetry Check (normalized): "
    call report_symmetry(symB_norm, 6)

    print *, "Second Subdiagonal of B:"
    do i = 3, 6
        print *, B(i, i - 2)
    end do

    print *, "Third Row of B:"
    do i = 1, 6
        print *, B(3, i)
    end do

    print *, "======================================"
    print *, "C: "
    
    C(:,:) = 0.0
    
    do i = 1, 1200
        do j = 1, 1200
            sum = 0.0
            do k = 1, 6
                sum = sum + A(i, k) * A(j, k)
            end do
            C(i, j) = sum
        end do
    end do

    ! extract the diagonal elements of C
    open (newunit=io, file="C_diag.dat", status="new", action="write")
    do i = 1, 1200
        write(io, *) C(i, i)
    end do
    close(io)

    ! Symmetry check (not normalized)
    do i = 1, 6
        do j = 1, 6
            symC(i, j) = abs(C(i, j) - C(j,i))
        end do
    end do 
    
    print*, "Symmetry Check (not normalized): "
    call report_symmetry(symC, 1200)

    ! symmetry check (normalized)
    do i = 1, 1200
        do j = 1, 1200
            symC_norm(i, j) = abs((C(i, j) - C(j,i)) / ((C(i, j) + C(j, i))/2))
        end do
    end do
    print *, "Symmetry Check (normalized): "
    call report_symmetry(symC_norm, 1200)


    open (newunit=io, file="C.dat", status="new", action="write")

    print *, "Second Subdiagonal of C:"
    print *, "Write to file C.dat"
    do i = 3, 1200
        write(io, *) C(i, i - 2)
    end do

    print *, "Third Row of C: "
    print *, "Write to file C.dat"
    do i = 1, 1200
        write(io, *) C(3, i)
    end do
    close(io)
    end program hw0



    subroutine print_matrix(array, n, m)
    implicit none 
    real, intent(in) :: array(n,m)
    integer, intent(in) :: n,m
    integer :: i
    do i = 1,n
    print*, array(:,i)
    end do
    end subroutine print_matrix
    
    subroutine report_symmetry(symB, n)
    implicit none
    real, intent(in) :: symB(n, n)
    integer, intent(in) :: n
    real :: mean, rms, max
    integer :: i, j

    mean = 0.0
    rms = 0.0
    max = symB(1, 1)
    do i = 1, n
        do j = 1, n
            mean = mean + symB(i, j)
            rms = rms + symB(i, j) ** 2
            if (symB(i, j) > max) then
                max = symB(i, j)
            end if
        end do
    end do
    mean = mean / (n * n)
    rms = sqrt(rms / (n * n))

    print*, "mean: ", mean
    print*, "rms: ", rms
    print*, "max: ", max
    end subroutine report_symmetry
