f(x) = 8*16*2500000000
g(x) = 8*68000000000*x
h(x) = 16*2500000000
i(x) = 68000000000*x
set xrange [0:1000000]
set yrange [0:100000000000000000]
set logscale x 2
set logscale y 2
unset key
plot h(x) with lines linestyle 1, \
     i(x) with lines linestyle 2, \
     f(x) with lines linestyle 3, \
     g(x) with lines linestyle 4, \
     "naive-parallel-one-sided.plotData"
