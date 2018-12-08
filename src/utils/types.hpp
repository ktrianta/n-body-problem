#ifndef TYPES_HPP
#define TYPES_HPP

#include <string>

#ifdef MPI_IS_USED
#include <mpi.h>
#endif

namespace sim {
    // Basic data type, (probably) float or double
    typedef double data_type;

#ifdef MPI_IS_USED
    const MPI_Datatype mpi_data_type = MPI_DOUBLE;
#endif

    // Global constants
    const data_type g = 1;
    const data_type e = 0.001;
    const data_type e2 = e * e;


    // Simulation parameters
    struct Parameters {
        // Default parameters as constants. Not static as there shouldn't exist
        // more than one instance of this struct in the end.
        const size_t def_n = 10;                        // default number of particles
        const size_t def_s = 5;                         // default number of time steps
        const data_type def_dt = 0.0001;                // default time step duration
        const std::string def_out_filename = "output";
        const std::string def_out_dirname = ".";

        size_t n;                                       // number of particles
        size_t s;                                       // number of time steps
        data_type dt;                                   // time step duration
        data_type theta;                                // barnes hut theta

        std::string in_filename;
        std::string out_filename;
        std::string out_dirname;

        Parameters()
            : n(def_n), s(def_s), dt(def_dt), out_filename(def_out_filename),
              out_dirname(def_out_dirname)
        {}

    };
}

#endif  // TYPES_HPP
