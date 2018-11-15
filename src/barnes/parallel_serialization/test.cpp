// classes example
#include <iostream>
#include <mpi.h>
using namespace std;

class Rectangle {
    int width, height;
    int s, r;
  public:
    void set_values (int, int, int, int);
    int area() {return width*height;}
};

void Rectangle::set_values (int x, int y, int size, int rank) {
if (rank == 0)
{
  width = x;
  height = y;
  }
}

int main (int argc, char** argv) {

  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  Rectangle rect;
  rect.set_values (3,4, size, rank);
  cout << "area: " << rect.area() <<  endl;

  MPI_Finalize();
  return 0;
}
