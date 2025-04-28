program TestAlignment
    use types
    use io
    use metrics
    implicit none

    type(ProteinStructure) :: unbound_struct, bound_struct
    real(8) :: rmsd_before, rmsd_after
    integer :: i

    ! Загружаем структуры
    call read_pdb('unbound.pdb', unbound_struct)
    call read_pdb('bound.pdb', bound_struct)

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

    ! Выводим координаты сразу после загрузки
    print *, "Координаты unbound_struct после загрузки:"
    do i = 1, unbound_struct%num_residues
        print *, i, unbound_struct%residues(i)%atoms(1)%x, unbound_struct%residues(i)%atoms(1)%y, unbound_struct%residues(i)%atoms(1)%z
    end do
    print *, "Координаты bound_struct после загрузки:"
    do i = 1, bound_struct%num_residues
        print *, i, bound_struct%residues(i)%atoms(1)%x, bound_struct%residues(i)%atoms(1)%y, bound_struct%residues(i)%atoms(1)%z
    end do

    ! Вычисляем RMSD до выравнивания
    call calculate_rmsd(unbound_struct, bound_struct, rmsd_before)
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