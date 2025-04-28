module metrics
    use types
    implicit none
contains
    subroutine calculate_overlap(unbound_struct, bound_struct, mode_ic, derivs, overlap)
        type(ProteinStructure), intent(in) :: unbound_struct, bound_struct
        real(8), intent(in) :: mode_ic(:)
        real(8), intent(in) :: derivs(:, :, :)
        real(8), intent(out) :: overlap
        real(8), allocatable :: delta_r_exp(:, :), delta_r_mode(:, :)
        real(8) :: norm_exp, norm_mode, dot_product
        integer :: i, j, k, n_atoms

        n_atoms = unbound_struct%num_atoms
        allocate(delta_r_exp(n_atoms, 3), delta_r_mode(n_atoms, 3))

        ! Вектор реального перехода
        do i = 1, n_atoms
            delta_r_exp(i, 1) = bound_struct%residues(i)%atoms(1)%x - unbound_struct%residues(i)%atoms(1)%x
            delta_r_exp(i, 2) = bound_struct%residues(i)%atoms(1)%y - unbound_struct%residues(i)%atoms(1)%y
            delta_r_exp(i, 3) = bound_struct%residues(i)%atoms(1)%z - unbound_struct%residues(i)%atoms(1)%z
        end do

        ! Вектор смещения моды (линейное приближение)
        delta_r_mode = 0.0d0
        do i = 1, n_atoms
            do k = 1, size(mode_ic)
                do j = 1, 3
                    delta_r_mode(i, j) = delta_r_mode(i, j) + derivs(i, j, k) * mode_ic(k)
                end do
            end do
        end do

        ! Вычисляем перекрытие
        norm_exp = 0.0d0
        norm_mode = 0.0d0
        dot_product = 0.0d0
        do i = 1, n_atoms
            do j = 1, 3
                norm_exp = norm_exp + delta_r_exp(i, j)**2
                norm_mode = norm_mode + delta_r_mode(i, j)**2
                dot_product = dot_product + delta_r_exp(i, j) * delta_r_mode(i, j)
            end do
        end do
        norm_exp = sqrt(norm_exp)
        norm_mode = sqrt(norm_mode + 1.0d-10)
        overlap = dot_product / (norm_exp * norm_mode)
    end subroutine calculate_overlap

    subroutine calculate_rmsd(struct1, struct2, rmsd)
        type(ProteinStructure), intent(in) :: struct1, struct2
        real(8), intent(out) :: rmsd
        integer :: i
        real(8) :: dx, dy, dz, sum_sq

        sum_sq = 0.0d0
        do i = 1, struct1%num_residues
            dx = struct1%residues(i)%atoms(1)%x - struct2%residues(i)%atoms(1)%x
            dy = struct1%residues(i)%atoms(1)%y - struct2%residues(i)%atoms(1)%y
            dz = struct1%residues(i)%atoms(1)%z - struct2%residues(i)%atoms(1)%z
            sum_sq = sum_sq + dx**2 + dy**2 + dz**2
        end do
        rmsd = sqrt(sum_sq / struct1%num_residues)
    end subroutine calculate_rmsd
end module metrics