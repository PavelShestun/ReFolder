module modify_structure
    use types
    use internal_coordinates
    implicit none
contains
    subroutine modify_structure_ic(unbound_struct, mode, omega, beta, phase, modified_struct)
        type(ProteinStructure), intent(in) :: unbound_struct
        real(8), intent(in) :: mode(:), omega, beta, phase
        type(ProteinStructure), intent(out) :: modified_struct
        real(8) :: alpha, kB_Teff
        real(8), allocatable :: new_coords(:)
        integer :: i, n_dof

        n_dof = 2 * unbound_struct%num_residues
        modified_struct = unbound_struct

        ! Вычисляем амплитуду
        kB_Teff = 1.0d0  ! Упрощение
        alpha = beta * sqrt(2.0d0 * kB_Teff / (omega**2 + 1.0d-10))

        ! Новые внутренние координаты
        allocate(new_coords(n_dof))
        new_coords = unbound_struct%internal_coords + alpha * mode * cos(phase)

        ! Обновляем структуру
        modified_struct%internal_coords = new_coords
        call internal_to_cartesian(modified_struct, new_coords)
    end subroutine modify_structure_ic
end module modify_structure