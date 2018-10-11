#!/usr/bin/gnuplot
#
# Generate an animated spiral
#
# AUTHOR: Hagen Wierstorf

reset

# png
set terminal gif animate
set output 'test.gif'
set xrange[0:1.5]
set yrange[0:1.5]
n=5
unset key
do for [i=1:300] {
    plot  'output.dat' every ::n*(i-1)+1::n*i 
}
