#include <unistd.h>
#include "args.hpp"


void readArgs(int argc, char** argv, sim::Parameters& params) {
    int c;

    while ((c = getopt (argc, argv, "n:t:s:i:d:")) != -1) {
        switch (c) {
            case 'n':
                params.n = atoi(optarg);
                break;
            case 't':
                params.t = atof(optarg);
                break;
            case 's':
                params.dt = atof(optarg);
                break;
            case 'i':
                params.in_filename = optarg;
                break;
            case 'd':
                params.out_dirname = optarg;
                break;
            default:
                break;
        }
    }
}
