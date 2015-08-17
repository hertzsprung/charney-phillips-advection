set multiplot layout 3,1
set view map
set cbrange [0:1]
set xrange [-150e3:150e3]
set yrange [0:25e3]
splot 'results' index 0 using 2:3:4 with image
splot 'results' index 1 using 2:3:4 with image
splot 'results' index 2 using 2:3:4 with image
