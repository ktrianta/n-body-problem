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
     x = px;y = py;z = pz;
     w = pw;h = ph;t = pt;
     index = -1;leaf = true;
     cum_size = 0;child = 0;
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


	MPI_Win win;
	MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);

	int NT;
	int NNPT;

	Treenode ** t = new Treenode*[NNPT];
	for (int i=0; i<NNPT; i++)
		Treenode*[i] = new Treenode[NT];



	MPI_Alloc_mem(sizeof(struct Treenode), MPI_INFO_NULL, &t);
	MPI_Win_attach(win,t,sizeof(struct Treenode));
 
	MPI_Alloc_mem(sizeof(struct Treenode), MPI_INFO_NULL, &t1);
	MPI_Win_attach(win,t1,sizeof(struct Treenode));
 
	MPI_Aint my_displ[2];
	MPI_Aint (*target_displ)[2] = new MPI_Aint[comm_size][2];
	MPI_Get_address(t,&my_displ[0]);
	MPI_Get_address(t1,&my_displ[1]);

	MPI_Allgather(&my_displ, 2, MPI_AINT, target_displ, 2, MPI_AINT, MPI_COMM_WORLD);
	if (rank == 0){
		t[0].x = 1.23; t[0].y = 4.1; t[0].z = 1.41; t[0].w = 25.1;
		t[0].h = 1.21; t[0].t = 7.2;
		t[0].mass = 27.2; t[0].massCenter[0] = 1.1;t[0].massCenter[1] = 1.2;
		t[0].massCenter[2] = 1.3;
		t[0].index = 10; t[0].leaf = true; t[0].cum_size = 31; t[0].child = 41;
	}

	MPI_Win_fence(0,win);
	if (rank == 1){
		MPI_Get(t,1,treeNodeStruct,0,target_displ[0][0],1,treeNodeStruct,win);
	}	
	MPI_Win_fence(0,win);
	
	if (rank == 1){
		std::cout<<t[0].x<<" "<<t[0].y<<" "<<t[0].z<<" "<<t[0].w<<" "<<std::endl;
		std::cout<<t[0].h<<" "<<t[0].t<<" "<<t[0].mass<<" "<<t[0].massCenter[0]<<" "<<std::endl;
		std::cout<<t[0].massCenter[1]<<" "<<t[0].massCenter[2]<<" "<<t[0].index<<" "<<t[0].leaf<<" "<<std::endl;
		std::cout<<t[0].cum_size<<" "<<t[0].child<<std::endl;
	}
	
	
	
	
	MPI_Win_detach(win,t);free(t);
	MPI_Win_free(&win);
	free(target_displ);
	MPI_Type_free(&treeNodeStruct);
	MPI_Finalize();

	return 0;
}
