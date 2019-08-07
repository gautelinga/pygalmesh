#ifndef ORIENT_SURFACE_MESH_HPP
#define ORIENT_SURFACE_MESH_HPP

#include <memory>
#include <string>

namespace pygalmesh {

void orient_surface_mesh(
    const std::string & infile,
    const std::string & outfile,
    const bool verbose = true
    );

} // namespace pygalmesh

#endif // ORIENT_SURFACE_MESH_HPP
