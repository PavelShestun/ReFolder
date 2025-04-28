module matrices
    use types
    use potential
    use internal_coordinates
    implicit none
contains
    subroutine calculate_hessian_ic(structure, hessian)
        type(ProteinStructure), intent(in) :: structure
        real(8), allocatable, intent(out) :: hessian(:, :)
        real(8) :: delta_q = 1.0d-3
        real(8), allocatable :: derivs_plus(:), derivs_minus(:)
        integer :: i, j, n_dof

        n_dof = 2 * structure%num_residues
        allocate(hessian(n_dof, n_dof))
        allocate(derivs_plus(n_dof), derivs_minus(n_dof))

        do i = 1, n_dof
            ! q_i + delta_q
            structure%internal_coords(i) = structure%internal_coords(i) + delta_q
            call calculate_first_derivatives(structure, structure%internal_coords, derivs_plus)

            ! q_i - delta_q
            structure%internal_coords(i) = structure%internal_coords(i) - 2.0d0 * delta_q
            call calculate_first_derivatives(structure, structure%internal_coords, derivs_minus)

            ! Восстановить q_i
            structure%internal_coords(i) = structure%internal_coords(i) + delta_q

            ! Численное дифференцирование
            do j = 1, n_dof
                hessian(j, i) = (derivs_plus(j) - derivs_minus(j)) / (2.0d0 * delta_q)
            end do
        end do
    end subroutine calculate_hessian_ic

    subroutine calculate_kinetic_ic(structure, kinetic)
        type(ProteinStructure), intent(in) :: structure
        real(8), allocatable, intent(out) :: kinetic(:, :)
        real(8), allocatable :: derivs(:, :, :)
        integer :: i, j, k, a, n_dof, n_atoms

        n_dof = 2 * structure%num_residues
        n_atoms = structure%num_atoms
        allocate(kinetic(n_dof, n_dof))
        call compute_derivatives_cartesian(structure, derivs)

        kinetic = 0.0d0
        do i = 1, n_dof
            do j = 1, n_dof
                do a = 1, n_atoms
                    do k = 1, 3
                        kinetic(i, j) = kinetic(i, j) + &
                            structure%residues(a)%atoms(1)%mass * derivs(a, k, i) * derivs(a, k, j)
                    end do
                end do
            end do
        end do
    end subroutine calculate_kinetic_ic
end module matrices