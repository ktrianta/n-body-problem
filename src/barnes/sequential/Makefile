objects = seqBarnesHut.o initialization.o quadtree.o
 
CXX = g++
CXXFLAGS = -std=c++11
CXXFLAGS += -Wall -Wextra -Wpedantic
CXXFLAGS += -g3 -O3
 
seqBarnesHut: ${objects}
	${CXX} ${CXXFLAGS} -o $@ $^
quadtree.o: quadtree.cpp
initialization.o: ../utils/initialization.cpp
	${CXX} ${CXXFLAGS} -c $<
seqBarnesHut.o: seqBarnesHut.cpp
	${CXX} ${CXXFLAGS} -c $<
 
clean:
	rm *.o seqBarnesHut 
clean-all:
	rm *.o seqBarnesHut *.gif *.dat
