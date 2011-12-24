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

// Left ventricle
// Origin
double l_o_x = 0.0;
double l_o_y = 0.0;
double l_o_z = 0.0;

// Major endocardium axes
double l_a_endo = 2.5;       // z
double l_b_endo = 1.1;       // x
double l_c_endo = l_b_endo;  // y

// Myocardium thicknesses
double l_a_trunc = 1.4;
double l_t_equator_b = 0.6;
double l_t_equator_c = l_t_equator_b;
double l_t_apex = 0.3;

// Major epicardium axes
double l_a_epi = l_a_endo + l_t_apex;
double l_b_epi = l_b_endo + l_t_equator_b;
double l_c_epi = l_c_endo + l_t_equator_c;

//Right ventricle
// Major endocardium axes
double r_a_endo = 2.2;
double r_b_endo = 1.4;
double r_c_endo = r_b_endo;

// Myocardium thicknesses
double r_a_trunc = 0.75;
double r_t_equator_b = 0.25;
double r_t_equator_c = r_t_equator_b;
double r_t_apex = 0.25;

//Origin
double r_o_x = l_o_x - 1.0;
double r_o_y = l_o_y + 0.0;
double r_o_z = l_o_z + l_a_trunc - r_a_trunc;

// Major epicardium axes
double r_a_epi = r_a_endo + r_t_apex;
double r_b_epi = r_b_endo + r_t_equator_b;
double r_c_epi = r_c_endo + r_t_equator_c;

// Function
FT ellipsoid (const Point& p,
	      double o_x, double o_y, double o_z,
	      double a, double b, double c)
{
    return sqrt(b*b*c*c*(p.z() - o_z)*(p.z() - o_z) +
		c*c*a*a*(p.x() - o_x)*(p.x() - o_x) +
		a*a*b*b*(p.y() - o_y)*(p.y() - o_y)) - a*b*c;
}

FT left_ventricle (const Point& p)
{
    if((ellipsoid(p, l_o_x, l_o_y, l_o_z, l_a_epi,  l_b_epi,  l_c_epi)  <= 0.0 &&
	ellipsoid(p, l_o_x, l_o_y, l_o_z, l_a_endo, l_b_endo, l_c_endo) >= 0.0 &&
	p.z() <= l_a_trunc) ||
       (ellipsoid(p, r_o_x, r_o_y, r_o_z, r_a_epi,  r_b_epi,  r_c_epi)  <= 0.0 &&
	ellipsoid(p, r_o_x, r_o_y, r_o_z, r_a_endo, r_b_endo, r_c_endo) >= 0.0 &&
	ellipsoid(p, l_o_x, l_o_y, l_o_z, l_a_endo, l_b_endo, l_c_endo) >  0.0 &&
	p.z() <= r_a_trunc))
	return -1;
    else
	return 1;
}

#include <cmath>

int main()
{
  // Domain
  Mesh_domain domain(left_ventricle, K::Sphere_3(Point(0, 0, 0), 9.));

  // Mesh criteria
  double h = 0.2;
  Mesh_criteria criteria(edge_size = h,
                         facet_angle = 25, facet_size = h,
                         cell_radius_edge_ratio = 2, cell_size = h);

  // Mesh generation without feature preservation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_features());

  std::ofstream medit_file("biventricle.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();
  c3t3.clear();

  return 0;
}
