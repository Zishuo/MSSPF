#set xlabel ""
set xlabel "Link indexes"
#set width 0.2
set xrange [0:26]
set yrange [0:22]
plot "cost239load.txt" u 1:2 w lp t "Normal", "" u 1:3 w lp t "Recalculation", "" u 1:4 w lp t "MSSPF", "" u 1:5 w lp t "Notvia"