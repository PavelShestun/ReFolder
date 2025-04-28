# Загружаем структуры
load unbound.pdb, unbound
load bound.pdb, bound
load unbound_aligned.pdb, unbound_aligned
load bound_aligned.pdb, bound_aligned

# Устанавливаем цвета
color blue, unbound
color red, bound
color cyan, unbound_aligned
color magenta, bound_aligned

# Отображаем структуры как cartoon
as cartoon, unbound
as cartoon, bound
as cartoon, unbound_aligned
as cartoon, bound_aligned

# Центрируем вид
center

# Увеличиваем масштаб для лучшего обзора
zoom
