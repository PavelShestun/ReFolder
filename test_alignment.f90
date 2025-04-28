program TestAlignment
    use types
    use io
    use metrics
    implicit none

    type(ProteinStructure) :: unbound_struct, bound_struct, unbound_orig, bound_orig
    real(8) :: rmsd_before, rmsd_after
    integer :: i

    ! Загружаем структуры
    call read_pdb('unbound.pdb', unbound_struct)
    call read_pdb('bound.pdb', bound_struct)

    ! Сохраняем копии исходных структур
    unbound_orig = unbound_struct
    bound_orig = bound_struct

    ! Выводим информацию о структурах
    print *, "Количество остатков в unbound:", unbound_struct%num_residues
    print *, "Количество остатков в bound:", bound_struct%num_residues
    print *, "Количество атомов в unbound:", unbound_struct%num_atoms
    print *, "Количество атомов в bound:", bound_struct%num_atoms

    ! Проверяем, совпадает ли количество остатков
    if (unbound_struct%num_residues /= bound_struct%num_residues) then
        print *, "Ошибка: структуры имеют разное количество остатков!"
        stop
    end if

    ! Вычисляем RMSD до выравнивания, используя исходные структуры
    call calculate_rmsd(unbound_orig, bound_orig, rmsd_before)
    print *, "RMSD до выравнивания:", rmsd_before

    ! Выравниваем структуры
    call align_structures(unbound_struct, bound_struct)

    ! Вычисляем RMSD после выравнивания
    call calculate_rmsd(unbound_struct, bound_struct, rmsd_after)
    print *, "RMSD после выравнивания:", rmsd_after

    ! Сохраняем выровненные структуры
    call write_pdb('unbound_aligned.pdb', unbound_struct)
    call write_pdb('bound_aligned.pdb', bound_struct)

    print *, "Выровненные структуры сохранены в unbound_aligned.pdb и bound_aligned.pdb"
end program TestAlignment