program ReFolder
    use types
    use io
    use internal_coordinates
    use matrices
    use eigenproblem
    use modify_structure
    use metrics
    implicit none

    type(ProteinStructure) :: unbound_struct, bound_struct, modified_struct, best_struct
    real(8), allocatable :: internal_coords(:), hessian(:, :), kinetic(:, :)
    real(8), allocatable :: modes(:, :), frequencies(:), derivs(:, :, :)
    real(8), allocatable :: trial_coords(:)  ! Для хранения временных внутренних координат
    integer :: i, j, k, n_dof, mode1, mode2, num_clashes, min_clashes
    real(8) :: beta1, beta2, overlap, best_overlap, best_beta1, best_beta2
    real(8) :: rmsd, best_rmsd
    real(8), parameter :: beta_step = 0.1d0
    real(8), parameter :: beta_max = 2.0d0
    real(8), parameter :: clash_threshold = 2.5d0

    ! Загружаем структуры
    call read_pdb('unbound.pdb', unbound_struct)
    call read_pdb('bound.pdb', bound_struct)

    ! Определяем внутренние координаты
    call define_internal_coordinates(unbound_struct, internal_coords)

    ! Вычисляем Hessian и кинетическую матрицу
    n_dof = 2 * unbound_struct%num_residues
    call calculate_hessian_ic(unbound_struct, hessian)
    call calculate_kinetic_ic(unbound_struct, kinetic)

    ! Решаем задачу на собственные значения
    call solve_eigenproblem(hessian, kinetic, modes, frequencies)

    ! Вычисляем производные декартовых координат по внутренним
    call compute_derivatives_cartesian(unbound_struct, derivs)

    ! Копируем начальную структуру для модификации
    modified_struct = unbound_struct
    best_struct = unbound_struct
    allocate(trial_coords(n_dof))

    ! Перебираем комбинации двух мод
    best_overlap = 0.0d0
    best_rmsd = 1.0d+10
    min_clashes = 1000000
    do mode1 = 1, min(5, n_dof - 1)
        do mode2 = mode1 + 1, min(6, n_dof)
            print *, "Testing modes:", mode1, mode2
            do i = -int(beta_max / beta_step), int(beta_max / beta_step)
                beta1 = i * beta_step
                do j = -int(beta_max / beta_step), int(beta_max / beta_step)
                    beta2 = j * beta_step

                    ! Модифицируем внутренние координаты
                    do k = 1, n_dof
                        trial_coords(k) = unbound_struct%internal_coords(k) + &
                                          beta1 * modes(k, mode1) + &
                                          beta2 * modes(k, mode2)
                    end do

                    ! Переводим в декартовы координаты только для текущей попытки
                    modified_struct%internal_coords = trial_coords
                    call internal_to_cartesian(modified_struct, modified_struct%internal_coords)

                    ! Вычисляем перекрытие
                    call calculate_overlap(modified_struct, bound_struct, &
                        beta1 * modes(:, mode1) + beta2 * modes(:, mode2), derivs, overlap)

                    ! Вычисляем RMSD_Ca-Ca
                    call calculate_rmsd_ca(modified_struct, bound_struct, rmsd)

                    ! Проверяем стерические столкновения
                    call check_clashes(modified_struct, clash_threshold, num_clashes)

                    ! Сохраняем лучшие значения (но не структуру)
                    if (overlap > best_overlap) then
                        best_overlap = overlap
                        best_beta1 = beta1
                        best_beta2 = beta2
                        ! Сохраняем только внутренние координаты
                        best_struct%internal_coords = trial_coords
                    end if
                    if (rmsd < best_rmsd) then
                        best_rmsd = rmsd
                    end if
                    if (num_clashes < min_clashes) then
                        min_clashes = num_clashes
                    end if
                end do
            end do

            print *, "Modes:", mode1, mode2
            print *, "Best overlap:", best_overlap
            print *, "Best beta1, beta2:", best_beta1, best_beta2
            print *, "Best RMSD_Ca-Ca:", best_rmsd
            print *, "Min clashes:", min_clashes
        end do
    end do

    ! Переводим в декартовы координаты только лучшую структуру
    call internal_to_cartesian(best_struct, best_struct%internal_coords)

    ! Сохраняем лучшую структуру
    call write_pdb('best_structure.pdb', best_struct)
    print *, "Лучшая структура сохранена в best_structure.pdb"
end program ReFolder