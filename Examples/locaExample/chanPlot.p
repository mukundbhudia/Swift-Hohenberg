set term post eps enhanced color
set output "locaChanProblem.eps"
set style line 1 lt 1 lw 6 linecolor rgb "#3D21A3"
set ylabel "T_{max}" font ",19"
set xlabel "{/Symbol a}" font ",19"
set xtics font ",16"
set ytics font ",16"
set xrange[0.0:5.0]
set yrange[-2.0:28.0]
set key top left font ",16"

plot "ChanContinuation.dat" using 1:2:(1.0) title 'Chan bifurcation' smooth bezier with line ls 1