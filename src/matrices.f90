module matrices
    use types
    use potential
    use internal_coordinates
    implicit none

contains

subroutine calculate_hessian_ic(structure, hessian)
    type(ProteinStructure), intent(in) :: structure
    real(8), intent(out) :: hessian(:, :)
    real(8) :: delta_q = 1.0d-6
    real(8) :: e_plus, e_minus, e_0
    real(8), allocatable :: temp_coords(:)
    integer :: i, j, n_dof
    type(ProteinStructure) :: temp_struct

    n_dof = 2 * structure%num_residues
    hessian = 0.0d0

    ! Создаём копию структуры
    temp_struct = structure
    allocate(temp_coords(n_dof))
    temp_coords = structure%internal_coords

    ! Вычисляем энергию в начальной точке
    call calculate_energy(temp_struct, e_0)

    ! Численное дифференцирование
    do i = 1, n_dof
        ! Положительный шаг
        temp_coords(i) = temp_coords(i) + delta_q
        temp_struct%internal_coords = temp_coords
        call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
        call calculate_energy(temp_struct, e_plus)

        ! Отрицательный шаг
        temp_coords(i) = temp_coords(i) - 2.0d0 * delta_q
        temp_struct%internal_coords = temp_coords
        call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
        call calculate_energy(temp_struct, e_minus)

        ! Восстанавливаем координаты
        temp_coords(i) = temp_coords(i) + delta_q
        temp_struct%internal_coords = temp_coords

        ! Вторая производная
        hessian(i, i) = (e_plus + e_minus - 2.0d0 * e_0) / (delta_q**2)

        do j = i + 1, n_dof
            ! Код для смешанных производных (если нужно)
        end do
    end do

    deallocate(temp_coords)
end subroutine calculate_hessian_ic

subroutine calculate_kinetic_ic(structure, kinetic)
    type(ProteinStructure), intent(in) :: structure
    real(8), intent(out) :: kinetic(:, :)
    integer :: n_dof

    n_dof = 2 * structure%num_residues
    kinetic = 0.0d0

    ! Для простоты предполагаем единичную кинетическую матрицу
    call identity_matrix(kinetic)
end subroutine calculate_kinetic_ic

subroutine identity_matrix(matrix)
    real(8), intent(out) :: matrix(:, :)
    integer :: i, n

    n = size(matrix, 1)
    matrix = 0.0d0
    do i = 1, n
        matrix(i, i) = 1.0d0
    end do
end subroutine identity_matrix

end module matrices