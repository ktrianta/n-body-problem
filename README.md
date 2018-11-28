# N Body Problem

## Building the Project

```bash
# Build out of source in build dir
mkdir build && cd build
cmake ..

# Run makefiles and install binaries in bin dir
make && make install
cd ..

# Run the tests:
!TODO

# Run some binary (e.g. naive sequential):
cd bin/naive-sequential
./prog

# If you encounter any problem detele dirs and try again
rm -rf build
rm -rf bin
```

## Execution Options

* -n: number of particles
* -t: duration of the simulation
* -s: duration of the simulation's timestep
* -i: input file containing the initial configuration of the particles

## Plotting Options

Run plot.text as 'gnuplot -e "n=#P" plot.text'
where #P should be the same as the -n execution option

## Setting up the Environment on the Cluster
```bash
# Install OpenMPI in $HOME/opt
cd $HOME
mkdir opt

wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.0.tar.gz
tar xf openmpi-4.0.0.tar.gz

cd openmpi-4.0.0
./configure --prefix=$HOME/opt
make -j 4 all
make install

# Set the new path in .bashrc
PATH=$HOME/opt/bin:$PATH

# Source .bashrc
source .bashrc
```

## Loading required Modules on the Cluster
```bash
./cluster.sh
```
## Installing LibLSB
```bash
#Download Library from https://spcl.inf.ethz.ch/Research/Performance/LibLSB/
tar -xf liblsb-0.2.2.tar.gz
cd liblsb-0.2.2
./configure && make && make install

#Now you want to go to the source code of our project 
#and add the following to the local CMakeLists.txt e.g.:
vi src/naive/parallel/CMakeLists.txt

#add target_link_libraries(naive-parallel lsb) after target_link_libraries(naive-parallel ${MPI_LIBRARIES})

```
