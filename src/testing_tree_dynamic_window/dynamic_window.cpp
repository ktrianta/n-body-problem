#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
  MPI_Win win;
  int rank, comm_size;
  int *a;
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  /* int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) */
  MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  MPI_Alloc_mem(2*sizeof(int), MPI_INFO_NULL, &a);

  //a = (int *) malloc(2*sizeof(int));

  if (rank==0) { a[0] = 1; a[1] = 3; }
  if (rank==1) { a[0] = 5; a[1] = 7; }
  /* Here the size is known */
  int n = 2;

  MPI_Win_attach(win, a, n*sizeof(int));

  /*----------------------------------------------------------*/
  //the others do not know the address of the target memory!!
  //distribute it
  MPI_Aint my_displ;
  MPI_Aint * target_displ = (MPI_Aint *) malloc(sizeof(MPI_Aint)*comm_size);
  
  MPI_Get_address(a, &my_displ);
      
  MPI_Allgather(&my_displ, 1, MPI_AINT, target_displ, 1, MPI_AINT, MPI_COMM_WORLD);
  /*----------------------------------------------------------*/

    
  MPI_Win_fence(0,win);
  if (rank == 1)
  {
  /* MPI_Get(void *origin_addr, int origin_count,MPI_Datatype origin_dtype, int target_rank,
                    MPI_Aint target_disp, int target_count,MPI_Datatype target_dtype, MPI_Win win)*/
    MPI_Put(a, 2 ,MPI_INT, 0, target_displ[0], 2, MPI_INT, win);
  }

  printf("win: %li\n", (long int) win);
  MPI_Win_fence(0, win);

  if (rank==0) std::cout<<a[0]<<" "<<a[1]<<std::endl;

  MPI_Win_detach(win,a); free(a);
  MPI_Win_free(&win);
  free(target_displ);
  MPI_Finalize();
  return 0;
}
