#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Domain
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Mesh_domain_with_polyline_features_3<
  CGAL::Implicit_mesh_domain_3<Function,K> >              Mesh_domain;

// Polyline
typedef std::vector<Point>        Polyline;
typedef std::list<Polyline>       Polylines;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function1 (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN)) - 2; }

FT sphere_function2 (const Point& p)
{ return CGAL::squared_distance(p, Point(2, 0, 0)) - 1; }

FT sphere_function (const Point& p)
{
  if(sphere_function1(p) < 0 || sphere_function2(p) < 0)
    return -1;
  else
    return 1;
}

#include <cmath>

int main()
{
  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(sphere_function, K::Sphere_3(Point(1, 0, 0), 6.));

  // Mesh criteria
  double h = 0.1;
  Mesh_criteria criteria(edge_size = h,
                         facet_angle = 25, facet_size = h,
                         cell_radius_edge_ratio = 2, cell_size = h);

  // Mesh generation without feature preservation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_features());

  std::ofstream medit_file("leftventricle.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();
  c3t3.clear();

  return 0;
}
