#include<iostream>
#include<mpi.h>
struct Treenode {
	double x;double y;double z;
	double w;double h;double t;
	double mass;double massCenter[3];
	int index;
	bool leaf;
	size_t cum_size;size_t child;

	Treenode();
	void set(double px, double py, double pz, double pw, double ph, double pt);
	};

Treenode::Treenode()
    : x(0), y(0), z(0), w(0), h(0), t(0), index(-1), leaf(true), cum_size(0), child(0)
    {}

    
void Treenode::set(double px, double py, double pz, double pw, double ph, double pt) {
     x = px;
     y = py;
     z = pz;
     w = pw;
     h = ph;
     t = pt;
     index = -1;
     leaf = true;
     cum_size = 0;
     child = 0;
}


int main(int argc, char** argv){
	MPI_Init(&argc,&argv);

	int comm_size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	Treenode tree;
	MPI_Datatype treeNodeStruct;
	

	int blocklength[] = {10,1,1,2};
	MPI_Datatype old_types[] = {MPI_DOUBLE,MPI_INT,MPI_C_BOOL,MPI_UNSIGNED_LONG};

	MPI_Aint baseaddr,a1,a2,a3,a4;
	MPI_Get_address(&tree,&baseaddr);
	MPI_Get_address(&tree.x, &a1);
	MPI_Get_address(&tree.index,&a2);
	MPI_Get_address(&tree.leaf,&a3);
	MPI_Get_address(&tree.cum_size,&a4);

	MPI_Aint indices[] = {a1-baseaddr,a2-baseaddr,a3-baseaddr,a4-baseaddr};
		
	MPI_Type_create_struct(4,blocklength,indices,old_types,&treeNodeStruct);
	MPI_Type_commit(&treeNodeStruct);

	
	if (rank == 0){
		tree.x = 1.23; tree.y = 4.1; tree.z = 1.41; tree.w = 25.1;
		tree.h = 1.21; tree.t = 7.2; 
		tree.mass = 27.2; tree.massCenter[0] = 1.1;tree.massCenter[1] = 1.2;
		tree.massCenter[2] = 1.3;
		tree.index = 10; tree.leaf = true; tree.cum_size = 31; tree.child = 41;
	}

	MPI_Bcast(&tree,1,treeNodeStruct,0,MPI_COMM_WORLD);	
	
	if (rank == 1){
		std::cout<<tree.x<<" "<<tree.y<<" "<<tree.z<<" "<<tree.w<<" "<<std::endl;
		std::cout<<tree.h<<" "<<tree.t<<" "<<tree.mass<<" "<<tree.massCenter[0]<<" "<<std::endl;
		std::cout<<tree.massCenter[1]<<" "<<tree.massCenter[2]<<" "<<tree.index<<" "<<tree.leaf<<" "<<std::endl;
		std::cout<<tree.cum_size<<" "<<tree.child<<std::endl;
	}


	MPI_Type_free(&treeNodeStruct);
	MPI_Finalize();

	return 0;
}


