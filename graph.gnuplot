set term pdf enhanced

set output "graph.pdf"

set title "Random Brownian path with 1 million steps"
set xlabel "p"
set ylabel "p-variation"
set logscale y
set ytics nomirror
set y2label "time"
unset logscale y2
set y2tics
plot "graph.log" using 1:2 with lines title "p-var",\
	"graph.log" using 1:3 with lines title "p-var^{1/p}",\
	"graph.log" using 1:4 with lines title "secs" axes x1y2,\
	"graph.log" using 1:5 with lines title "R mthd secs" axes x1y2
