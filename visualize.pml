load unbound_aligned.pdb
load bound_aligned.pdb
color blue, unbound_aligned
color red, bound_aligned
show cartoon
align unbound_aligned, bound_aligned
center
zoom
png alignment_visualization.png, width=800, height=600, dpi=300
quit
