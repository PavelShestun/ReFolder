program ReFolder
    use types
    use io
    use internal_coordinates
    use matrices
    use eigenproblem
    use modify_structure
    use metrics
    implicit none

    type(ProteinStructure) :: unbound_struct, bound_struct, modified_struct
    real(8), allocatable :: internal_coords(:), hessian(:, :), kinetic(:, :)
    real(8), allocatable :: eigenvalues(:), modes(:, :), derivs(:, :, :)
    real(8) :: beta_factor, min_rmsd_m_b, optimal_beta, overlap
    integer :: num_modes_to_compute, mode_index
    real(8), parameter :: beta_min = 0.0d0, beta_max = 2.0d0, beta_step = 0.1d0

    ! 1. Загрузить структуры
    call read_pdb('unbound.pdb', unbound_struct)
    call read_pdb('bound.pdb', bound_struct)

    ! 2. Определить внутренние координаты
    call define_internal_coordinates(unbound_struct, internal_coords)

    ! 3. Рассчитать H и T
    call calculate_hessian_ic(unbound_struct, hessian)
    call calculate_kinetic_ic(unbound_struct, kinetic)

    ! 4. Решить задачу на собственные значения
    num_modes_to_compute = min(20, 2 * unbound_struct%num_residues)
    call solve_eigenproblem_ic(hessian, kinetic, num_modes_to_compute, eigenvalues, modes)

    ! 5. Вычислить производные для перекрытия
    call compute_derivatives_cartesian(unbound_struct, derivs)

    ! 6. Для каждой моды найти оптимальное RMSD_m-b
    do mode_index = 1, num_modes_to_compute
        min_rmsd_m_b = huge(1.0d0)
        optimal_beta = 0.0d0

        ! Цикл по амплитудам beta
        do beta_factor = beta_min, beta_max, beta_step
            call modify_structure_ic(unbound_struct, modes(:, mode_index), &
                                     sqrt(eigenvalues(mode_index)), beta_factor, 0.0d0, &
                                     modified_struct)
            call calculate_rmsd(modified_struct, bound_struct, min_rmsd_m_b)
            if (min_rmsd_m_b < min_rmsd_m_b) then
                min_rmsd_m_b = min_rmsd_m_b
                optimal_beta = beta_factor
            end if
        end do

        ! Рассчитать перекрытие
        call calculate_overlap(unbound_struct, bound_struct, modes(:, mode_index), derivs, overlap)

        ! Вывод результатов
        print *, "Mode:", mode_index, "Min RMSD_m-b:", min_rmsd_m_b, &
                 "Optimal Beta:", optimal_beta, "Overlap:", overlap
    end do

    ! TODO: Добавить цикл для комбинации двух мод с варьированием фаз
end program ReFolder