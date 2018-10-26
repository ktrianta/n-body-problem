# N Body Problem

## Building The Project

```bash
cd build
cmake ..
make && make install
cd ..

# Run the tests:
!TODO

# Run some binary (e.g. naive sequential):
cd bin/naive-sequential
./prog
```

## Execution Options

* -n: number of particles
* -t: duration of the simulation
* -s: duration of the simulation's timestep
* -i: input file containing the initial configuration of the particles

## Plotting Options

Run plot.text as 'gnuplot -e "n=#P" plot.text'
where #P should be the same as the -n execution option
