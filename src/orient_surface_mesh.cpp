#define CGAL_SURFACE_MESHER_VERBOSE 1

#include "orient_surface_mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

// IO
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <vector>
#include <fstream>
//#include <iostream>

namespace pygalmesh {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;

  void orient_surface_mesh(
       const std::string & infile,
       const std::string & outfile,
       const bool verbose
       )
  {
    std::ifstream input(infile);

    std::vector<K::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;

    if (!input.good()) {
      std::stringstream msg;
      msg << "Cannot read input file \"" << infile << "\"" << std::endl;
      throw std::runtime_error(msg.str());
    }
    if (!CGAL::read_OFF(input, points, polygons) || points.empty()) {
      std::stringstream msg;
      msg << "Cannot read OFF file \"" << infile << "\"" << std::endl;
      throw std::runtime_error(msg.str());
    }
    input.close();
    
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

    Polyhedron mesh;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

    int index = 0;
    for (Polyhedron::Face_iterator fb=mesh.facets_begin(),
	   fe=mesh.facets_end(); fb != fe; ++fb)
      fb->id() = index++;

    if (CGAL::is_closed(mesh))
      CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(mesh);

    std::ofstream output(outfile);
    output << mesh;
    output.close();

    return;
  }
} // namespace pygalmesh
