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

// Parameters defining the ventricular geometry
double a_endo = 2.5;     // z
double b_endo = 1.1;     // x
double c_endo = b_endo;  // y

double a_trunc = 1.4;
double t_equator_b = 0.6;
double t_equator_c = t_equator_b;
double t_apex = 0.3;

double a_epi = a_endo + t_apex;
double b_epi = b_endo + t_equator_b;
double c_epi = c_endo + t_equator_c;

double l_o_x = 0.0;
double l_o_y = 0.0;
double l_o_z = 0.0;


// Function
FT epicardium (const Point& p)
{
    return
	sqrt(b_epi*b_epi*c_epi*c_epi*p.z()*p.z() +
	     c_epi*c_epi*a_epi*a_epi*p.x()*p.x() +
	     a_epi*a_epi*b_epi*b_epi*p.y()*p.y()) -
	a_epi*b_epi*c_epi;
}

FT endocardium (const Point& p)
{
    return
	sqrt(b_endo*b_endo*c_endo*c_endo*p.z()*p.z() +
	     c_endo*c_endo*a_endo*a_endo*p.x()*p.x() +
	     a_endo*a_endo*b_endo*b_endo*p.y()*p.y()) -
	a_endo*b_endo*c_endo;
}

FT sphere_function (const Point& p)
{
    if(epicardium(p) < 0 && endocardium(p) > 0 && p.z() < a_trunc)
	return -1;
    else
	return 1;
}

#include <cmath>

int main()
{
  // Domain
  Mesh_domain domain(sphere_function, K::Sphere_3(Point(0, 0, 0), 9.));

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
