#define CGAL_SURFACE_MESHER_VERBOSE 1

#include "remesh_surface.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

namespace pygalmesh {

  // Domain
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

  // Polyhedron type
  typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

  // Triangulation
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // To avoid verbose function and named parameters call
  using namespace CGAL::parameters;
  
  void remesh_surface(
		      const std::string & infile,
		      const std::string & outfile,
		      const double edge_size,
		      const double facet_angle,
		      const double facet_size,
		      const double facet_distance,
		      const bool verbose
		      )
  {
    Polyhedron poly;
    std::ifstream input(infile);
    input >> poly;
    if (!input.good()) {
      std::stringstream msg;
      msg << "Cannot read input file \"" << infile << "\"" << std::endl;
      throw std::runtime_error(msg.str());
    }
    if (!CGAL::is_triangle_mesh(poly)) {
      std::stringstream msg;
      msg << "Input geometry is not triangulated." << std::endl;
      throw std::runtime_error(msg.str());
    }
    input.close();
    
    // Create vector with only one element: pointer to polyhedron
    std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);

    // Polyhedral domain with one polyhedron and no bounding polyhedron.
    Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());

    // Get sharp features??
    domain.detect_features();  // Includes detection of borders?

    // Criteria
    Mesh_criteria criteria(CGAL::parameters::edge_size = edge_size,
			   CGAL::parameters::facet_angle = facet_angle,
			   CGAL::parameters::facet_size = facet_size,
			   CGAL::parameters::facet_distance = facet_distance);

    // Mesh generation
    if (!verbose) {
      // suppress output
      std::cout.setstate(std::ios_base::failbit);
      std::cerr.setstate(std::ios_base::failbit);
    }
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(
       domain, criteria, no_perturb(), no_exude());

    // Output facets (not oriented)
    std::ofstream off_file(outfile);
    c3t3.output_boundary_to_off(off_file);
    off_file.close();

    if (!verbose) {
      std::cout.clear();
      std::cerr.clear();
    }

    return;
  }
} // namespace pygalmesh
