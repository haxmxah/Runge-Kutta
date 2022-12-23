set terminal png
set output "P8-20-21-ﬁg1.png"
set term png size 1100, 600


set title  "Gràfica de les solucions ɸ_{1,2,3,4} sense normalitzar" 
set title font "Times New Roman,20"
set xlabel font "Times New Roman,12"
set ylabel font "Times New Roman,12"
set xlabel "X[Å]" 
set ylabel "ɸ_E(1)"
!set key outside
set grid xtics
set grid ytics
L = 16.
set xrange [-L/2.:L/4.]


set key top left
set key font "Times New Roman,12"
set key samplen 1


plot "auxiliar.dat" i 0 u 2:3 w l lw 2 t"E_1[eV]","auxiliar.dat" i 1 u 2:3 w l lw 2  t"E_2[eV]","auxiliar.dat" i 2 u 2:3 w l lw 2  t"E_3[eV]","auxiliar.dat" i 3 u 2:3 w l lw 2 t"E_4[eV]"