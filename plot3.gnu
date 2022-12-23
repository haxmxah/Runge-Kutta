set terminal png
set output "P8-20-21-ﬁg3.png"
set term png size 1100, 600


set title  "Gràfica d'autovectors normalitzats" 
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


plot "auxiliar.dat" i 4 u 1:2 w l lw 2 t"E_1[eV] = 1.38  ","auxiliar.dat" i 5 u 1:2 w l lw 2  t"E_2[eV] = 4.14","auxiliar.dat" i 6 u 1:2 w l lw 2  t"E_3[eV]=6.9"