#include <iostream>
#include <fstream>
#include <vector>
#include "../../utils/types.hpp"
#include "../../utils/initialization.hpp"
#include "quadtree.hpp"
#include <unistd.h>
#include <mpi.h>

using namespace std;

void writeDataToFile(int N, sim::data_type (*r)[3], sim::data_type (*u)[3], ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r[i][0] << "   "
             << r[i][1] << "   "
             << r[i][2] << "   "
             << u[i][0] << "   "
             << u[i][1] << "   "
             << u[i][2] << "\n";
    }
}

int main(int argc, char** argv)
{
    // ----- MPI ----- //
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
// coordinates for tab8096
    xc = -2.5;
    yc = -2.5;
    zc = -2.5;
    w2 = 22.5;
    h2 = 17.5;
    t2 = 22.5;

    int c;
    int N = 5;      // the number of particles
    double theta = 0.0;      // the number of particles
    sim::data_type T = 1;
    sim::data_type dt = 0.00001;
    string filename;

    while ((c = getopt (argc, argv, "n:t:s:i:h:")) != -1)
    {
        switch (c)
        {
            case 'n':
                N = atoi(optarg);
                break;
            case 't':
                T = atof(optarg);
                break;
            case 's':
                dt = atof(optarg);
                break;
            case 'i':
                filename = optarg;
                break;
            case 'h':
                theta = atof(optarg);
                break;
        }
    }

    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
    vector<sim::data_type> m(N, 1.0/N);

    if (rank == 0)
    {
        if (!filename.empty()) {
            ifstream ifile;
            ifile.open(filename);

            for (int i = 0; i < N; i++) {
                ifile >> m[i] >> r[i][0] >> r[i][1] >> r[i][2] >> u[i][0] >> u[i][1] >> u[i][2];
            }
        } else {
              initializePositionOnSphere(N, r);
        }
    }

    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&r[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);

    ofstream file;
    if (rank == 0)
    {
        file.open("output.dat");
        writeDataToFile(N, r, u, file);
    }

    QuadTree tree = QuadTree(r, m, N, xc, yc, zc, w2, h2, t2);

    const int Ntimesteps = T/dt + 1;
    // Local Declarations
    int *local_N = new int[size];
    int *local_Nx3 = new int[size];
    int local_N_int = N/size;
    int rem = N - local_N_int * size;
    int counter = 0;
    for (int i=0;i<size;i++)
    {
        local_N[i] = local_N_int;
        if (counter < rem)
        {
                local_N[i] += 1;
                counter ++;
        }
        local_Nx3[i] = local_N[i]*3;
    }

    sim::data_type (*u_local)[3] = new sim::data_type[local_N[rank]][3];
    sim::data_type (*r_local)[3] = new sim::data_type[local_N[rank]][3];
    sim::data_type (*a_local)[3] = new sim::data_type[local_N[rank]][3];

    int *offset = new int[size];
    offset[0] = 0;
    int *offset_x3 = new int[size];
    offset_x3[0] = 0;

    for (int i = 1; i < size; i++)
    {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x3[i] = offset_x3[i-1] + local_Nx3[i-1];
    }

    for (int i = 0; i < local_N[rank]; i++)
    {
        u_local[i][0] = u[offset[rank]+i][0];
        u_local[i][1] = u[offset[rank]+i][1];
        u_local[i][2] = u[offset[rank]+i][2];
        r_local[i][0] = r[offset[rank]+i][0];
        r_local[i][1] = r[offset[rank]+i][1];
        r_local[i][2] = r[offset[rank]+i][2];
        a_local[i][0] = a[offset[rank]+i][0];
        a_local[i][1] = a[offset[rank]+i][1];
        a_local[i][2] = a[offset[rank]+i][2];
    }

    for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
    {
        int jLocal=j-offset[rank];
        tree.computeAcceleration(j,jLocal, r, a_local, sim::g, theta); // NEEDS TO BE PARALLELIZED
    }


    double start = 0;
    double end = 0;
    double time = 0;

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
        {
            int jLocal=j-offset[rank];
            u_local[jLocal][0] += 0.5 * a_local[jLocal][0] * dt;
            u_local[jLocal][1] += 0.5 * a_local[jLocal][1] * dt;
            u_local[jLocal][2] += 0.5 * a_local[jLocal][2] * dt;
            r_local[jLocal][0] += u_local[jLocal][0] * dt;
            r_local[jLocal][1] += u_local[jLocal][1] * dt;
            r_local[jLocal][2] += u_local[jLocal][2] * dt;
        }

        MPI_Allgatherv(&(r_local[0][0]),local_N[rank]*3,MPI_DOUBLE,&(r[0][0]),local_Nx3, offset_x3, MPI_DOUBLE,MPI_COMM_WORLD);

        for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
        {
            int jLocal=j-offset[rank];
            start = MPI_Wtime();
            tree.computeAcceleration(j, jLocal, r, a_local, sim::g, theta);
            end = MPI_Wtime();
            time += end - start;

            u_local[jLocal][0] += 0.5 * a_local[jLocal][0] * dt;
            u_local[jLocal][1] += 0.5 * a_local[jLocal][1] * dt;
            u_local[jLocal][2] += 0.5 * a_local[jLocal][2] * dt;
            a_local[jLocal][0] = 0;
            a_local[jLocal][1] = 0;
            a_local[jLocal][2] = 0;
        }

        QuadTree tree = QuadTree(r, m, N, xc, yc, zc, w2, h2, t2);
        if (t % 200 == 0 && rank == 0)
        {
            writeDataToFile(N, r, u, file);
        }
    }

    printf("time= %f \n", time);
    MPI_Finalize();
    return 0;
}
