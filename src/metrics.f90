module metrics
    use types
    implicit none

contains

subroutine calculate_rmsd(struct1, struct2, rmsd)
    type(ProteinStructure), intent(in) :: struct1, struct2
    real(8), intent(out) :: rmsd
    integer :: i
    real(8) :: sum_sq_diff

    sum_sq_diff = 0.0d0
    do i = 1, struct1%num_atoms
        sum_sq_diff = sum_sq_diff + &
            (struct1%residues(i)%atoms(1)%x - struct2%residues(i)%atoms(1)%x)**2 + &
            (struct1%residues(i)%atoms(1)%y - struct2%residues(i)%atoms(1)%y)**2 + &
            (struct1%residues(i)%atoms(1)%z - struct2%residues(i)%atoms(1)%z)**2
    end do
    rmsd = sqrt(sum_sq_diff / real(struct1%num_atoms, 8))
end subroutine calculate_rmsd

subroutine calculate_rmsd_ca(struct1, struct2, rmsd)
    type(ProteinStructure), intent(in) :: struct1, struct2
    real(8), intent(out) :: rmsd
    integer :: i
    real(8) :: sum_sq_diff

    sum_sq_diff = 0.0d0
    do i = 1, struct1%num_residues
        sum_sq_diff = sum_sq_diff + &
            (struct1%residues(i)%atoms(1)%x - struct2%residues(i)%atoms(1)%x)**2 + &
            (struct1%residues(i)%atoms(1)%y - struct2%residues(i)%atoms(1)%y)**2 + &
            (struct1%residues(i)%atoms(1)%z - struct2%residues(i)%atoms(1)%z)**2
    end do
    rmsd = sqrt(sum_sq_diff / real(struct1%num_residues, 8))
end subroutine calculate_rmsd_ca

subroutine calculate_overlap(struct1, struct2, displacement, derivs, overlap)
    type(ProteinStructure), intent(in) :: struct1, struct2
    real(8), intent(in) :: displacement(:), derivs(:, :, :)
    real(8), intent(out) :: overlap
    integer :: i, j, k
    real(8) :: diff(3), proj

    overlap = 0.0d0
    do i = 1, struct1%num_residues
        diff(1) = struct2%residues(i)%atoms(1)%x - struct1%residues(i)%atoms(1)%x
        diff(2) = struct2%residues(i)%atoms(1)%y - struct1%residues(i)%atoms(1)%y
        diff(3) = struct2%residues(i)%atoms(1)%z - struct1%residues(i)%atoms(1)%z
        proj = 0.0d0
        do j = 1, size(displacement)
            do k = 1, 3
                proj = proj + derivs(i, j, k) * displacement(j) * diff(k)
            end do
        end do
        overlap = overlap + proj
    end do
    overlap = overlap / real(struct1%num_residues, 8)
end subroutine calculate_overlap

subroutine check_clashes(structure, threshold, num_clashes)
    type(ProteinStructure), intent(in) :: structure
    real(8), intent(in) :: threshold
    integer, intent(out) :: num_clashes
    integer :: i, j
    real(8) :: dist

    num_clashes = 0
    do i = 1, structure%num_residues - 1
        do j = i + 1, structure%num_residues
            dist = sqrt( &
                (structure%residues(i)%atoms(1)%x - structure%residues(j)%atoms(1)%x)**2 + &
                (structure%residues(i)%atoms(1)%y - structure%residues(j)%atoms(1)%y)**2 + &
                (structure%residues(i)%atoms(1)%z - structure%residues(j)%atoms(1)%z)**2 &
            )
            if (dist < threshold) then
                num_clashes = num_clashes + 1
            end if
        end do
    end do
end subroutine check_clashes

subroutine align_structures(struct1, struct2)
    type(ProteinStructure), intent(inout) :: struct1
    type(ProteinStructure), intent(in) :: struct2
    real(8) :: rotation(3, 3), translation(3)
    integer :: i, j
    real(8) :: new_coords(3)

    call kabsch(struct1, struct2, rotation, translation)

    ! Применяем поворот и трансляцию к struct1
    do i = 1, struct1%num_residues
        do j = 1, size(struct1%residues(i)%atoms)
            new_coords(1) = struct1%residues(i)%atoms(j)%x
            new_coords(2) = struct1%residues(i)%atoms(j)%y
            new_coords(3) = struct1%residues(i)%atoms(j)%z
            new_coords = matmul(rotation, new_coords) + translation
            struct1%residues(i)%atoms(j)%x = new_coords(1)
            struct1%residues(i)%atoms(j)%y = new_coords(2)
            struct1%residues(i)%atoms(j)%z = new_coords(3)
        end do
    end do
end subroutine align_structures

subroutine kabsch(struct1, struct2, rotation, translation)
    type(ProteinStructure), intent(in) :: struct1, struct2
    real(8), intent(out) :: rotation(3, 3), translation(3)
    real(8) :: centroid1(3), centroid2(3)
    real(8) :: coords1(3, struct1%num_atoms), coords2(3, struct2%num_atoms)
    real(8) :: A(3, 3), U(3, 3), VT(3, 3), S(3)
    real(8) :: d
    integer :: i, j

    ! Вычисляем центроиды
    centroid1 = 0.0d0
    centroid2 = 0.0d0
    do i = 1, struct1%num_atoms
        centroid1(1) = centroid1(1) + struct1%residues(i)%atoms(1)%x
        centroid1(2) = centroid1(2) + struct1%residues(i)%atoms(1)%y
        centroid1(3) = centroid1(3) + struct1%residues(i)%atoms(1)%z
        centroid2(1) = centroid2(1) + struct2%residues(i)%atoms(1)%x
        centroid2(2) = centroid2(2) + struct2%residues(i)%atoms(1)%y
        centroid2(3) = centroid2(3) + struct2%residues(i)%atoms(1)%z
    end do
    centroid1 = centroid1 / real(struct1%num_atoms, 8)
    centroid2 = centroid2 / real(struct2%num_atoms, 8)

    ! Центрируем координаты
    do i = 1, struct1%num_atoms
        coords1(1, i) = struct1%residues(i)%atoms(1)%x - centroid1(1)
        coords1(2, i) = struct1%residues(i)%atoms(1)%y - centroid1(2)
        coords1(3, i) = struct1%residues(i)%atoms(1)%z - centroid1(3)
        coords2(1, i) = struct2%residues(i)%atoms(1)%x - centroid2(1)
        coords2(2, i) = struct2%residues(i)%atoms(1)%y - centroid2(2)
        coords2(3, i) = struct2%residues(i)%atoms(1)%z - centroid2(3)
    end do

    ! Вычисляем ковариационную матрицу A
    A = matmul(coords1, transpose(coords2))

    ! SVD разложение A = U * S * VT
    call svd(A, U, S, VT)

    ! Вычисляем матрицу поворота R = U * VT
    rotation = matmul(U, VT)

    ! Проверяем определитель R
    d = determinant(rotation)
    if (d < 0.0d0) then
        ! Корректируем U, умножая последнюю колонку на -1
        U(:, 3) = -U(:, 3)
        rotation = matmul(U, VT)
    end if

    ! Вычисляем вектор трансляции
    translation = centroid2 - matmul(rotation, centroid1)
end subroutine kabsch

function determinant(matrix) result(det)
    real(8), intent(in) :: matrix(3, 3)
    real(8) :: det

    det = matrix(1,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) - &
          matrix(1,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) + &
          matrix(1,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))
end function determinant

subroutine svd(A, U, S, VT)
    real(8), intent(inout) :: A(3, 3)
    real(8), intent(out) :: U(3, 3), S(3), VT(3, 3)
    integer :: info, lwork
    real(8) :: work_query(1)
    real(8), allocatable :: work(:)

    U = 0.0d0
    S = 0.0d0
    VT = 0.0d0
    info = 0

    lwork = -1
    call dgesvd('A', 'A', 3, 3, A, 3, S, U, 3, VT, 3, work_query, lwork, info)
    if (info /= 0) then
        print *, "Error in DGESVD workspace query, INFO = ", info
        stop
    end if

    lwork = int(work_query(1))
    allocate(work(lwork))

    call dgesvd('A', 'A', 3, 3, A, 3, S, U, 3, VT, 3, work, lwork, info)
    if (info /= 0) then
        print *, "Error in DGESVD, INFO = ", info
        stop
    end if

    deallocate(work)
end subroutine svd

end module metrics