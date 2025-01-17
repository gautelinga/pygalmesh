#ifndef GENERATE_FROM_VOX_HPP
#define GENERATE_FROM_VOX_HPP

#include "sizing_field.hpp"

#include <memory>
#include <string>
#include <vector>

namespace pygalmesh {

void generate_mesh_from_vox(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb = true,
    const bool exude = true,
    const double edge_size = 0.0,  // std::numeric_limits<double>::max(),
    const double facet_angle = 0.0,
    const double facet_size = 0.0,
    const double facet_distance = 0.0,
    const double cell_radius_edge_ratio = 0.0,
    const double cell_size = 0.0,
    const bool verbose = true
    );

void generate_periodic_mesh_from_vox(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb = true,
    const bool exude = true,
    const double edge_size = 0.0,  // std::numeric_limits<double>::max(),
    const double facet_angle = 0.0,
    const double facet_size = 0.0,
    const double facet_distance = 0.0,
    const double cell_radius_edge_ratio = 0.0,
    const double cell_size = 0.0,
    const int number_of_copies_in_output = 1,
    const bool verbose = true
    );

void generate_periodic_mesh_from_vox_with_sizing_field(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const std::string & cell_size_file,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb = true,
    const bool exude = true,
    const double edge_size = 0.0,
    const double facet_angle = 0.0,
    const double facet_size = 0.0,
    const double facet_distance = 0.0,
    const double cell_radius_edge_ratio = 0.0,
    const int number_of_copies_in_output = 1,
    const bool verbose = true
    );

} // namespace pygalmesh

#endif // GENERATE_FROM_VOX_HPP
