module internal_coordinates
    use types
    implicit none
    real(8), parameter :: bond_length = 3.8d0  ! Длина связи Cα-Cα (Å)
    real(8), parameter :: bond_angle = 1.911d0  ! Угол Cα-Cα-Cα (≈109.5° в радианах)
contains
    subroutine define_internal_coordinates(structure, internal_coords)
        type(ProteinStructure), intent(in) :: structure
        real(8), allocatable, intent(out) :: internal_coords(:)
        ! Уже определено в read_pdb (заглушка)
        allocate(internal_coords(2 * structure%num_residues))
        internal_coords = structure%internal_coords
    end subroutine define_internal_coordinates

    subroutine internal_to_cartesian(structure, internal_coords)
        type(ProteinStructure), intent(inout) :: structure
        real(8), intent(in) :: internal_coords(:)
        integer :: i, n_residues
        real(8) :: phi, psi
        real(8) :: x(3), y(3), z(3)  ! Локальная система координат
        real(8) :: r(3)              ! Позиция текущего атома
        real(8) :: r_prev(3), r_prev2(3)  ! Позиции предыдущих атомов
        real(8) :: theta, norm

        n_residues = structure%num_residues

        ! Инициализация первых трёх атомов (нужны для задания начальной системы координат)
        if (n_residues < 3) stop "Слишком мало остатков для реконструкции"

        ! Атом 1: в начале координат
        structure%residues(1)%atoms(1)%x = 0.0d0
        structure%residues(1)%atoms(1)%y = 0.0d0
        structure%residues(1)%atoms(1)%z = 0.0d0
        r_prev2 = [0.0d0, 0.0d0, 0.0d0]

        ! Атом 2: вдоль оси X
        structure%residues(2)%atoms(1)%x = bond_length
        structure%residues(2)%atoms(1)%y = 0.0d0
        structure%residues(2)%atoms(1)%z = 0.0d0
        r_prev = [bond_length, 0.0d0, 0.0d0]

        ! Атом 3: в плоскости XY с учетом bond_angle
        structure%residues(3)%atoms(1)%x = bond_length * (1.0d0 + cos(bond_angle))
        structure%residues(3)%atoms(1)%y = bond_length * sin(bond_angle)
        structure%residues(3)%atoms(1)%z = 0.0d0
        r = [bond_length * (1.0d0 + cos(bond_angle)), bond_length * sin(bond_angle), 0.0d0]

        ! Цикл по остальным атомам
        do i = 4, n_residues
            ! Определяем локальную систему координат на основе трёх предыдущих атомов
            r_prev2 = [structure%residues(i-3)%atoms(1)%x, &
                       structure%residues(i-3)%atoms(1)%y, &
                       structure%residues(i-3)%atoms(1)%z]
            r_prev = [structure%residues(i-2)%atoms(1)%x, &
                      structure%residues(i-2)%atoms(1)%y, &
                      structure%residues(i-2)%atoms(1)%z]
            r = [structure%residues(i-1)%atoms(1)%x, &
                 structure%residues(i-1)%atoms(1)%y, &
                 structure%residues(i-1)%atoms(1)%z]

            ! Ось x: направление от (i-2) к (i-1)
            x = r - r_prev
            norm = sqrt(sum(x**2))
            x = x / norm

            ! Ось z: перпендикуляр к плоскости (i-3)-(i-2)-(i-1)
            y = r_prev - r_prev2
            z(1) = x(2) * y(3) - x(3) * y(2)
            z(2) = x(3) * y(1) - x(1) * y(3)
            z(3) = x(1) * y(2) - x(2) * y(1)
            norm = sqrt(sum(z**2))
            if (norm < 1.0d-6) then
                z = [0.0d0, 0.0d0, 1.0d0]  ! Если плоскость вырождена, выбираем z произвольно
            else
                z = z / norm
            end if

            ! Ось y: перпендикуляр к x и z
            y(1) = z(2) * x(3) - z(3) * x(2)
            y(2) = z(3) * x(1) - z(1) * x(3)
            y(3) = z(1) * x(2) - z(2) * x(1)
            norm = sqrt(sum(y**2))
            y = y / norm

            ! Используем phi и psi для поворота
            phi = internal_coords(2*(i-1)-1) * 3.14159265359d0 / 180.0d0  ! phi в радианах
            psi = internal_coords(2*(i-1)) * 3.14159265359d0 / 180.0d0      ! psi в радианах

            ! Позиция нового атома с учётом phi и psi
            theta = bond_angle
            ! Сначала поворот вокруг оси z на phi
            r = r + bond_length * (cos(theta) * x + sin(theta) * (cos(phi) * y + sin(phi) * z))
            ! Затем поворот вокруг новой оси y на psi
            ! Обновляем локальную систему координат для psi
            x = r - [structure%residues(i-1)%atoms(1)%x, &
                     structure%residues(i-1)%atoms(1)%y, &
                     structure%residues(i-1)%atoms(1)%z]
            norm = sqrt(sum(x**2))
            x = x / norm
            y(1) = z(2) * x(3) - z(3) * x(2)
            y(2) = z(3) * x(1) - z(1) * x(3)
            y(3) = z(1) * x(2) - z(2) * x(1)
            norm = sqrt(sum(y**2))
            y = y / norm
            z(1) = x(2) * y(3) - x(3) * y(2)
            z(2) = x(3) * y(1) - x(1) * y(3)
            z(3) = x(1) * y(2) - x(2) * y(1)
            norm = sqrt(sum(z**2))
            z = z / norm
            ! Применяем поворот на psi
            r = r + bond_length * (cos(psi) * x + sin(psi) * z)

            structure%residues(i)%atoms(1)%x = r(1)
            structure%residues(i)%atoms(1)%y = r(2)
            structure%residues(i)%atoms(1)%z = r(3)
        end do
    end subroutine internal_to_cartesian

    subroutine compute_derivatives_cartesian(structure, derivs)
        type(ProteinStructure), intent(in) :: structure
        real(8), allocatable, intent(out) :: derivs(:, :, :)
        real(8) :: delta_q = 1.0d-3  ! Шаг для численного дифференцирования
        type(ProteinStructure) :: temp_struct
        real(8), allocatable :: coords_plus(:, :), coords_minus(:, :)
        integer :: i, j, k, n_dof, n_atoms

        n_dof = 2 * structure%num_residues
        n_atoms = structure%num_atoms
        allocate(derivs(n_atoms, 3, n_dof))
        allocate(coords_plus(n_atoms, 3), coords_minus(n_atoms, 3))

        ! Копируем структуру
        temp_struct = structure
        do i = 1, n_dof
            ! q_i + delta_q
            temp_struct%internal_coords = structure%internal_coords
            temp_struct%internal_coords(i) = temp_struct%internal_coords(i) + delta_q
            call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
            do j = 1, structure%num_residues
                coords_plus(j, 1) = temp_struct%residues(j)%atoms(1)%x
                coords_plus(j, 2) = temp_struct%residues(j)%atoms(1)%y
                coords_plus(j, 3) = temp_struct%residues(j)%atoms(1)%z
            end do

            ! q_i - delta_q
            temp_struct%internal_coords = structure%internal_coords
            temp_struct%internal_coords(i) = temp_struct%internal_coords(i) - delta_q
            call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
            do j = 1, structure%num_residues
                coords_minus(j, 1) = temp_struct%residues(j)%atoms(1)%x
                coords_minus(j, 2) = temp_struct%residues(j)%atoms(1)%y
                coords_minus(j, 3) = temp_struct%residues(j)%atoms(1)%z
            end do

            ! Численное дифференцирование
            do j = 1, n_atoms
                do k = 1, 3
                    derivs(j, k, i) = (coords_plus(j, k) - coords_minus(j, k)) / (2.0d0 * delta_q)
                end do
            end do
        end do
    end subroutine compute_derivatives_cartesian
end module internal_coordinates