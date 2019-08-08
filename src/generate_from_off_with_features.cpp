#include "generate_from_off_with_features.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

namespace pygalmesh {

  // Domain
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

  // Triangulation
  typedef CGAL::Mesh_triangulation_3<
    Mesh_domain, CGAL::Default, CGAL::Sequential_tag>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // Avoid verbose function and named parameters call
  using namespace CGAL::parameters;

  void generate_from_off_with_features(const std::string &infile,
				       const std::string &outfile,
				       const bool lloyd, const bool odt,
				       const bool perturb, const bool exude,
				       const double edge_size,
				       const double facet_angle,
				       const double facet_size,
				       const double facet_distance,
				       const double cell_radius_edge_ratio,
				       const double cell_size,
				       const bool verbose) {

    // Create polyhedron
    Polyhedron poly;

    // Read file
    std::ifstream input(infile);
    input >> poly;
    if (!input.good()) {
      std::stringstream msg;
      msg << "Cannot read input file \"" << infile << "\"" << std::endl;
      throw std::runtime_error(msg.str());
    }
    input.close();

    // Domain
    Mesh_domain domain(poly);

    // Get features
    domain.detect_features();

    // Define criteria
    Mesh_criteria criteria(CGAL::parameters::edge_size = edge_size,
			   CGAL::parameters::facet_angle = facet_angle,
			   CGAL::parameters::facet_size = facet_size,
			   CGAL::parameters::facet_distance = facet_distance,
			   CGAL::parameters::cell_radius_edge_ratio = cell_radius_edge_ratio,
			   CGAL::parameters::cell_size = cell_size);

    // Mesh generation
    if (!verbose) {
      // suppress output
      std::cerr.setstate(std::ios_base::failbit);
    }
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
					lloyd ? CGAL::parameters::lloyd() : no_lloyd(),
					odt ? CGAL::parameters::odt() : no_odt(),
					perturb ? CGAL::parameters::perturb() : no_perturb(),
					exude ? CGAL::parameters::exude() : no_exude());

    if (!verbose) {
      std::cerr.clear();
    }

    // Output
    std::ofstream output(outfile);
    c3t3.output_to_medit(output);
    output.close();

    return;
  }

} // namespace pygalmesh
