#define CGAL_MESH_3_VERBOSE 1

#include "generate_from_vox.hpp"

// periodic features
#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/number_type_config.h> // CGAL_PI
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>

// non-periodic
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>


namespace pygalmesh {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Labeled_mesh_domain_3<K> Periodic_mesh_domain;

typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<K>> Mesh_domain;
  
// Triangulation
// Periodic triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type PTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<PTr> PC3t3;
// Non-periodic triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  
// Mesh Criteria
// Periodic criteria
typedef CGAL::Mesh_criteria_3<PTr> Periodic_mesh_criteria;
typedef Periodic_mesh_criteria::Facet_criteria Periodic_facet_criteria;
typedef Periodic_mesh_criteria::Cell_criteria Periodic_cell_criteria;
// Non-periodic criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Iso_cuboid_3 Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


double*** read_vox(const std::string &infile,
		   double &Lx, double &Ly, double &Lz,
		   int &Nx, int &Ny, int &Nz,
		   const bool verbose){
  std::ifstream input(infile);
  input >> Lx >> Ly >> Lz;
  input >> Nx >> Ny >> Nz;

  if (verbose) {
    std::cout << "Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz << std::endl;
    std::cout << "Nx=" << Nx << ", Ny=" << Ny << ", Nz=" << Nz << std::endl;
  }
  
  double ***C = new double**[Nx];
  for (int i=0; i < Nx; ++i){
    C[i] = new double*[Ny];
    for (int j=0; j < Ny; ++j){
      C[i][j] = new double[Nz];
      for (int k=0; k < Nz; ++k){
	input >> C[i][j][k];
      }
    }
  }
  input.close();
  
  return C;
}

int imodulo(const int a, const int b) {
  return ((a % b) + b) % b;
}

double trilinear_intp(const double x,
		      const double y,
		      const double z,
		      double*** C,
		      const int Nx,
		      const int Ny,
		      const int Nz,
		      const double Lx,
		      const double Ly,
		      const double Lz,
		      const bool verbose){
  int ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi;
  double wx, wy, wz;
  double f000, f001, f010, f011, f100, f101, f110, f111;
  double f;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dz = Lz/Nz;

  int ix_est = floor(x/dx);
  int iy_est = floor(y/dy);
  int iz_est = floor(z/dz);
 
  ix_lo = imodulo(ix_est, Nx);
  ix_hi = imodulo(ix_lo + 1, Nx);
  iy_lo = imodulo(iy_est, Ny);
  iy_hi = imodulo(iy_lo + 1, Ny);
  iz_lo = imodulo(iz_est, Nz);
  iz_hi = imodulo(iz_lo + 1, Nz);

  if (verbose)
    std::cout << int(x/dx) << " -> " << ix_lo << " " << ix_hi << std::endl;
  
  wx = (x-dx*ix_est)/dx;
  wy = (y-dy*iy_est)/dy;
  wz = (z-dz*iz_est)/dz;

  if (verbose)
    std::cout << wx << " " << wy << " " << wz << std::endl;
    
  // Nodal values
  f000 = C[ix_lo][iy_lo][iz_lo];
  f001 = C[ix_lo][iy_lo][iz_hi];
  f010 = C[ix_lo][iy_hi][iz_lo];
  f011 = C[ix_lo][iy_hi][iz_hi];
  f100 = C[ix_hi][iy_lo][iz_lo];
  f101 = C[ix_hi][iy_lo][iz_hi];
  f110 = C[ix_hi][iy_hi][iz_lo];
  f111 = C[ix_hi][iy_hi][iz_hi];

  // Trilinear interpolaton
  f = f000*(1-wx)*(1-wy)*(1-wz)
    + f001*(1-wx)*(1-wy)*wz
    + f010*(1-wx)*wy*(1-wz)
    + f011*(1-wx)*wy*wz
    + f100*wx*(1-wy)*(1-wz)
    + f101*wx*(1-wy)*wz
    + f110*wx*wy*(1-wz)
    + f111*wx*wy*wz;

  return f;
}


void
generate_mesh_from_vox(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const bool lloyd,
    const bool odt,
    const bool perturb,
    const bool exude,
    const double edge_size,
    const double facet_angle,
    const double facet_size,
    const double facet_distance,
    const double cell_radius_edge_ratio,
    const double cell_size,
    const bool verbose
    )
{
  K::Iso_cuboid_3 cuboid(
      bounding_cuboid[0],
      bounding_cuboid[1],
      bounding_cuboid[2],
      bounding_cuboid[3],
      bounding_cuboid[4],
      bounding_cuboid[5]
      );

  // Read voxel data
  int Nx, Ny, Nz;
  double Lx, Ly, Lz;
  double*** C = read_vox(infile, Lx, Ly, Lz, Nx, Ny, Nz, verbose);
  
  // wrap domain
  const auto d = [&](K::Point_3 p) {
      return trilinear_intp(p.x(), p.y(), p.z(), C, Nx, Ny, Nz, Lx, Ly, Lz, verbose);
  };

  Mesh_domain cgal_domain = Mesh_domain::create_implicit_mesh_domain(d, cuboid);

  Mesh_criteria criteria(
      CGAL::parameters::edge_size=edge_size,
      CGAL::parameters::facet_angle=facet_angle,
      CGAL::parameters::facet_size=facet_size,
      CGAL::parameters::facet_distance=facet_distance,
      CGAL::parameters::cell_radius_edge_ratio=cell_radius_edge_ratio,
      CGAL::parameters::cell_size=cell_size
      );

  // Mesh generation
  if (!verbose) {
    // suppress output
    std::cerr.setstate(std::ios_base::failbit);
  }
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(
      cgal_domain,
      criteria,
      lloyd ? CGAL::parameters::lloyd() : CGAL::parameters::no_lloyd(),
      odt ? CGAL::parameters::odt() : CGAL::parameters::no_odt(),
      perturb ? CGAL::parameters::perturb() : CGAL::parameters::no_perturb(),
      exude ? CGAL::parameters::exude() : CGAL::parameters::no_exude()
      );
  if (!verbose) {
    std::cerr.clear();
  }

  // Output
  std::ofstream medit_file(outfile);
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  return;
}

void
generate_periodic_mesh_from_vox(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const bool lloyd,
    const bool odt,
    const bool perturb,
    const bool exude,
    const double edge_size,
    const double facet_angle,
    const double facet_size,
    const double facet_distance,
    const double cell_radius_edge_ratio,
    const double cell_size,
    const int number_of_copies_in_output,
    const bool verbose
    )
{
  K::Iso_cuboid_3 cuboid(
      bounding_cuboid[0],
      bounding_cuboid[1],
      bounding_cuboid[2],
      bounding_cuboid[3],
      bounding_cuboid[4],
      bounding_cuboid[5]
      );

  // Read voxel data
  int Nx, Ny, Nz;
  double Lx, Ly, Lz;
  double*** C = read_vox(infile, Lx, Ly, Lz, Nx, Ny, Nz, verbose);
  
  // wrap domain
  const auto d = [&](K::Point_3 p) {
      return trilinear_intp(p.x(), p.y(), p.z(), C, Nx, Ny, Nz, Lx, Ly, Lz, verbose);
  };

  Periodic_mesh_domain cgal_domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(d, cuboid);

  Periodic_mesh_criteria criteria(
      CGAL::parameters::edge_size=edge_size,
      CGAL::parameters::facet_angle=facet_angle,
      CGAL::parameters::facet_size=facet_size,
      CGAL::parameters::facet_distance=facet_distance,
      CGAL::parameters::cell_radius_edge_ratio=cell_radius_edge_ratio,
      CGAL::parameters::cell_size=cell_size
      );

  // Mesh generation
  if (!verbose) {
    // suppress output
    std::cerr.setstate(std::ios_base::failbit);
  }
  PC3t3 pc3t3 = CGAL::make_periodic_3_mesh_3<PC3t3>(
      cgal_domain,
      criteria,
      lloyd ? CGAL::parameters::lloyd() : CGAL::parameters::no_lloyd(),
      odt ? CGAL::parameters::odt() : CGAL::parameters::no_odt(),
      perturb ? CGAL::parameters::perturb() : CGAL::parameters::no_perturb(),
      exude ? CGAL::parameters::exude() : CGAL::parameters::no_exude()
      );
  if (!verbose) {
    std::cerr.clear();
  }

  // Output
  std::ofstream medit_file(outfile);
  CGAL::output_periodic_mesh_to_medit(medit_file, pc3t3, number_of_copies_in_output);
  medit_file.close();

  return;
}

// same but with sizing field in cell_size
// It'd be nice if we could replace this clumsy class by a simple function wrapper (like
// domain), but CGAL expects the type FT to be present. :(
// TODO file issue for that on <https://github.com/CGAL/cgal/issues>
class Sizing_field_wrapper_vox
{
public:
  typedef K::FT FT;

  Sizing_field_wrapper_vox(const std::string & cell_size_file, bool verbose): verbose(verbose) {
    S = read_vox(cell_size_file, Lx, Ly, Lz, Nx, Ny, Nz, verbose);
  }

  virtual ~Sizing_field_wrapper_vox() = default;

  K::FT operator()(const K::Point_3& p, const int, const Mesh_domain::Index&) const
  {
    auto out = trilinear_intp(p.x(), p.y(), p.z(), S, Nx, Ny, Nz, Lx, Ly, Lz, verbose);
    return out;
  }

private:
  double ***S;
  int Nx = 0;
  int Ny = 0;
  int Nz = 0;
  double Lx = 0.0;
  double Ly = 0.0;
  double Lz = 0.0;
  bool verbose = false;
};

void
generate_periodic_mesh_from_vox_with_sizing_field(
    const std::string & infile,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const std::string & cell_size_file,
    const bool lloyd,
    const bool odt,
    const bool perturb,
    const bool exude,
    const double edge_size,
    const double facet_angle,
    const double facet_size,
    const double facet_distance,
    const double cell_radius_edge_ratio,
    const int number_of_copies_in_output,
    const bool verbose
    )
{
  K::Iso_cuboid_3 cuboid(
      bounding_cuboid[0],
      bounding_cuboid[1],
      bounding_cuboid[2],
      bounding_cuboid[3],
      bounding_cuboid[4],
      bounding_cuboid[5]
      );

  // Read voxel data
  int Nx, Ny, Nz;
  double Lx, Ly, Lz;
  double*** C = read_vox(infile, Lx, Ly, Lz, Nx, Ny, Nz, verbose);
  
  // wrap domain
  const auto d = [&](K::Point_3 p) {
      return trilinear_intp(p.x(), p.y(), p.z(), C, Nx, Ny, Nz, Lx, Ly, Lz, verbose);
  };

  Periodic_mesh_domain cgal_domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(d, cuboid);

  Sizing_field_wrapper_vox cell_size(cell_size_file, verbose);
  
  Periodic_mesh_criteria criteria(
      CGAL::parameters::edge_size=edge_size,
      CGAL::parameters::facet_angle=facet_angle,
      CGAL::parameters::facet_size=facet_size,
      CGAL::parameters::facet_distance=facet_distance,
      CGAL::parameters::cell_radius_edge_ratio=cell_radius_edge_ratio,
      CGAL::parameters::cell_size=cell_size
      );

  // Mesh generation
  if (!verbose) {
    // suppress output
    std::cerr.setstate(std::ios_base::failbit);
  }
  PC3t3 pc3t3 = CGAL::make_periodic_3_mesh_3<PC3t3>(
      cgal_domain,
      criteria,
      lloyd ? CGAL::parameters::lloyd() : CGAL::parameters::no_lloyd(),
      odt ? CGAL::parameters::odt() : CGAL::parameters::no_odt(),
      perturb ? CGAL::parameters::perturb() : CGAL::parameters::no_perturb(),
      exude ? CGAL::parameters::exude() : CGAL::parameters::no_exude()
      );
  if (!verbose) {
    std::cerr.clear();
  }

  // Output
  std::ofstream medit_file(outfile);
  CGAL::output_periodic_mesh_to_medit(medit_file, pc3t3, number_of_copies_in_output);
  medit_file.close();

  return;
}


} // namespace pygalmesh
