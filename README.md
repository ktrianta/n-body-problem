# N Body Problem

## Building the Project


# codes
src/naive baseline implementation
src/barnes barnes-hut implementation


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

* -n [#particles] : number of particles
* -s [#steps] : number of simulation timesteps
* -t [time] : duration of a timestep (default 0.0001)
* -i [file] : input file containing the initial configuration of the particles
* -w : write particle position in output file every 200 timesteps
* -e : print the error between initial and final energy 


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
./configure --with-mpi --enable-sync && make && make install

#Now you want to go to the source code of our project 
#and add the following to the local CMakeLists.txt e.g.:
vi src/naive/parallel/CMakeLists.txt

#add target_link_libraries(naive-parallel lsb) 
#and set(CMAKE_CXX_COMPILER mpicxx)
#after target_link_libraries(naive-parallel ${MPI_LIBRARIES})

```
