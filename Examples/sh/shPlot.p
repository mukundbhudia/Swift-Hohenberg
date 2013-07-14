set term post eps enhanced color
set output "sh-example.eps"
set style line 1 lt 1 lw 6 linecolor rgb "#9f1a01"
set style line 2 lt 1 lw 6 linecolor rgb "#9f1a01"
set ylabel "T"
set xlabel "{/Symbol a}"
plot "shContinuation.dat" using 1:2 title 'sh bifurcation' with line ls 1, \
"shContinuation.dat" using 1:3 title '' with line ls 2, \
"shContinuation.dat" using 1:4 title '' with line ls 2, \
"shContinuation.dat" using 1:5 title '' with line ls 2