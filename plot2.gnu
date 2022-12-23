set terminal png
set output "P8-20-21-ﬁg2.png"
set term png size 1100, 600


set multiplot title  "Convergència del mètode" font "Times New Roman,20" layout 1,3

set title font "Times New Roman,30"
set xlabel font "Times New Roman,18"
set ylabel font "Times New Roman,18"
set xlabel "Iteracions" 
set ylabel "E[eV]"
!set key outside
set grid xtics
set grid ytics
set bmargin 5
set tmargin 2



set key top left
set key font "Times New Roman,12"
set key samplen 1


plot "auxiliar2.dat" i 0 u 1:2 w l lw 2 t"E_1[eV]"
plot "auxiliar2.dat" i 1 u 1:2 w l lw 2 t"E_2[eV]"
plot "auxiliar2.dat" i 2 u 1:2 w l lw 2 t"E_3[eV]"