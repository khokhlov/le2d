in  = "d"
out = "double.svg"

set term svg
set output out

set grid
set style data histogram
set style fill solid 0.3
set boxwidth 0.8
set ylabel "Time, s"
#unset xtics
set yrange[0:100]
set auto x
plot in using 2 ti col, '' u 3 ti col, '' u 4:xtic(1) ti col

