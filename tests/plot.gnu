set term png
set output "plot.png"
set style data lines
set grid
set xlabel "Grid size, MB"
set ylabel "Cycles rdtsc per node (one step)"

plot "le2d.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "No cache friendly", \
	"le2d_cf.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "Cache friendly", \
	"le2d_soa.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "SOA", \
	"le2d_sse.log" using ($1*$2/1024/1024*8*5):($5/$1/$2/$3) ti "SSE"
