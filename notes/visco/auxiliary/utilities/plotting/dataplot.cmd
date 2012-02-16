set yrange [-1.2:0.0001]
set xlabel "Time (days)"                           
set ylabel "Horizontal Stress (MPa)"
set terminal postscript enhanced color fontfile 'sfrm2074.pfb' "SFRM2074" 20.74 linewidth 8.0
set output "constrained-stress-evolution.eps"
plot 'constrained-stress-evolution.data' notitle with lines