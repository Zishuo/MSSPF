set style fill pattern
#set style histograms clustered
set style data histograms
plot "nsfnetdstr.txt" u 2 t "Normal", "" u 3 t "Recalculation", "" u 4 t "MSSPF", "" u 5 t "Notvia"