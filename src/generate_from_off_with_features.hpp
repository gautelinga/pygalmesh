#ifndef GENERATE_FROM_OFF_WITH_FEATURES_HPP
#define GENERATE_FROM_OFF_WITH_FEATURES_HPP

#include <string>

namespace pygalmesh {

void generate_from_off_with_features(
    const std::string &infile,
    const std::string &outfile,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb = false,
    const bool exude = false,
    const double edge_size = 0.0,
    const double facet_angle = 0.0,
    const double facet_size = 0.0,
    const double facet_distance = 0.0,
    const double cell_radius_edge_ratio = 0.0,
    const double cell_size = 0.0,
    const bool verbose = true
    );
} // namespace pygalmesh

#endif // GENERATE_FROM_OFF_WITH_FEATURES_HPP
