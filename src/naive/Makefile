objects = sequential.o initialization.o
 
CXX = g++
CXXFLAGS = -std=c++11
CXXFLAGS += -Wall -Wextra -Wpedantic
CXXFLAGS += -O3
 
sequential: ${objects}
	${CXX} ${CXXFLAGS} -o $@ $^
initialization.o: ../utils/initialization.cpp
	${CXX} ${CXXFLAGS} -c $<
sequential.o: sequential.cpp
	${CXX} ${CXXFLAGS} -c $<
 
clean:
	rm *.o sequential
clean-all:
	rm *.o sequential *.gif *.dat
