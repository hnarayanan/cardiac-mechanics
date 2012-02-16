set log xy
set xrange [20:900]
set xtics (20, 50, 100, 200, 500, 900)
set yrange [1:10]
set ytics (1, 2, 5, 10)
set xlabel "Time (s)"                           
set ylabel "Engineering Stress (MPa)"
#set terminal postscript enhanced color fontfile 'sfss2074.pfb' "SFSS2074" 20.74 linewidth 8.0
set terminal postscript enhanced color fontfile 'sfrm2074.pfb' "SFRM2074" 20.74 linewidth 8.0
#set terminal postscript enhanced color fontfile 'sfbx2074.pfb' "SFBX2074" 20.74 linewidth 8.0
set output "data-comparison-3.eps"
plot 'provenzano-expt-3.txt' title 'Experimental Data', 'nonlin-visco-3.txt' title 'Nonlinear Viscoelasticity' with lines

