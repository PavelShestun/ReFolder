module potential
    use types
    implicit none
    real(8), parameter :: characteristic_distance = 3.8d0
    real(8), parameter :: s_torsion = 1.0d0
contains
    subroutine calculate_force_constant(r0, Fij)
        real(8), intent(in) :: r0
        real(8), intent(out) :: Fij
        Fij = 1.0d0 / (1.0d0 + (r0 / characteristic_distance)**6)
    end subroutine calculate_force_constant

    subroutine calculate_potential_energy(structure, internal_coords, energy)
        type(ProteinStructure), intent(in) :: structure
        real(8), intent(in) :: internal_coords(:)
        real(8), intent(out) :: energy
        integer :: i, j, k
        real(8) :: r0, rij, Fij, dx, dy, dz
        real(8) :: term1, term2

        term1 = 0.0d0
        term2 = 0.0d0

        ! Первый член: упругая сетка
        do i = 1, structure%num_residues - 1
            do j = i + 1, structure%num_residues
                dx = structure%residues(i)%atoms(1)%x - structure%residues(j)%atoms(1)%x
                dy = structure%residues(i)%atoms(1)%y - structure%residues(j)%atoms(1)%y
                dz = structure%residues(i)%atoms(1)%z - structure%residues(j)%atoms(1)%z
                r0 = sqrt(dx**2 + dy**2 + dz**2)

                ! Текущее расстояние (упрощение: используем r0)
                rij = r0
                call calculate_force_constant(r0, Fij)
                term1 = term1 + Fij * (rij - r0)**2
            end do
        end do
        term1 = 0.5d0 * term1

        ! Второй член: торсионная жесткость
        do i = 1, structure%num_residues
            do k = 1, 2  ! phi и psi
                term2 = term2 + s_torsion * (internal_coords(2*i-2+k) - structure%internal_coords(2*i-2+k))**2
            end do
        end do
        term2 = 0.5d0 * term2

        energy = term1 + term2
    end subroutine calculate_potential_energy

    subroutine calculate_first_derivatives(structure, internal_coords, derivs)
        type(ProteinStructure), intent(in) :: structure
        real(8), intent(in) :: internal_coords(:)
        real(8), allocatable, intent(out) :: derivs(:)
        real(8) :: delta_q = 1.0d-3
        real(8) :: energy_plus, energy_minus
        integer :: i, n_dof
        type(ProteinStructure) :: temp_struct

        n_dof = 2 * structure%num_residues
        allocate(derivs(n_dof))
        temp_struct = structure

        do i = 1, n_dof
            ! q_i + delta_q
            temp_struct%internal_coords = internal_coords
            temp_struct%internal_coords(i) = temp_struct%internal_coords(i) + delta_q
            call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
            call calculate_potential_energy(temp_struct, temp_struct%internal_coords, energy_plus)

            ! q_i - delta_q
            temp_struct%internal_coords = internal_coords
            temp_struct%internal_coords(i) = temp_struct%internal_coords(i) - delta_q
            call internal_to_cartesian(temp_struct, temp_struct%internal_coords)
            call calculate_potential_energy(temp_struct, temp_struct%internal_coords, energy_minus)

            derivs(i) = (energy_plus - energy_minus) / (2.0d0 * delta_q)
        end do
    end subroutine calculate_first_derivatives
end module potential