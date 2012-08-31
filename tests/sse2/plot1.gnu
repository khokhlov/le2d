#set term png
#set output "plot1.png"
set style data lines
set grid
set xlabel "Grid size, MB"
set ylabel "Cycles rdtsc per node (one step)"

plot "le2d_f.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "float", \
	"le2d_d.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "double", \
	"le2d_d_sse.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "double, sse", \
	"le2d_f_sse.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "float, sse", \
	"le2d_d_avx.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "double, avx", \
	"le2d_f_avx.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "float, avx"
