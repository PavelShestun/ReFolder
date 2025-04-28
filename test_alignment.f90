program TestAlignment
    use types
    use io
    use metrics
    implicit none

    type(ProteinStructure) :: unbound_struct, bound_struct
    real(8) :: rmsd_before, rmsd_after

    ! Загружаем структуры
    call read_pdb('unbound.pdb', unbound_struct)
    call read_pdb('bound.pdb', bound_struct)

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