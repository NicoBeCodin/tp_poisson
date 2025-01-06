set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'resvec_plot.png'
set title "Residual Convergence"
set xlabel "Iteration"
set ylabel "Residual"
set logscale y
set grid
plot "RESVEC.dat" using ($0+1):1 with linespoints title "Residual"
