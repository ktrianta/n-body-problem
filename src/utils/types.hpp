#ifndef TYPES_HPP
#define TYPES_HPP

#include <string>

namespace sim {
    // Basic data type, (probably) float or double
    typedef double data_type;

    // Global constants
    const data_type g = 1;
    const data_type e = 0.001;
    const data_type e2 = e * e;


    // Simulation parameters
    struct Parameters {
        // Default parameters as constants. Not static as there shouldn't exist
        // more than one instance of this struct in the end.
        const size_t def_n = 10;
        const data_type def_t = 10;
        const data_type def_dt = 0.00001;
        const std::string def_out_filename = "output";
        const std::string def_out_dirname = ".";

        size_t n;                                   // number of particles
        data_type t;                                // simulation duration
        data_type dt;                               // time step duration

        std::string in_filename;
        std::string out_filename;
        std::string out_dirname;

        Parameters()
            : n(def_n), t(def_t), dt(def_dt), out_filename(def_out_filename),
              out_dirname(default_out_dirname)
        {}

    };
}

#endif  // TYPES_HPP
