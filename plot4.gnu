set terminal png
set output "P8-20-21-ﬁg4.png"

set title  "Gràfica d'autovectors normalitzats per diferents β" 
set title font "Times New Roman,20"
set xlabel font "Times New Roman,12"
set ylabel font "Times New Roman,12"
set xlabel "X[Å]" 
set ylabel "ɸ_E(1)"
!set key outside
set grid xtics
set grid ytics


set key top left
set key font "Times New Roman,12"
set key samplen 1


plot "auxiliar.dat" i 7 u 1:2 w l lw 2 t"β_1= 0  ","auxiliar.dat" i 8 u 1:2 w l lw 2  t"β_2 = 0.1","auxiliar.dat" i 9 u 1:2 w l lw 2  t"β_3=0.25"