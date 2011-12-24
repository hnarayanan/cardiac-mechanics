#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>

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

// Declare parameters defining the ventricular geometry (defined in parameters.ini)

// General
double h;

// Left ventricle
// Center
double l_c_x, l_c_y, l_c_z;
// Major endocardium axes
double l_a_endo, l_a_trunc, l_b_endo, l_c_endo;

// Myocardium thicknesses
double l_t_equator_b, l_t_equator_c, l_t_apex;

// Major epicardium axes
double l_a_epi, l_b_epi, l_c_epi;

//Right ventricle
// Center
double r_c_x, r_c_y, r_c_z;

// Major endocardium axes
double r_a_endo, r_a_trunc, r_b_endo, r_c_endo;

// Myocardium thicknesses
double r_t_equator_b, r_t_equator_c, r_t_apex;

// Major epicardium axes
double r_a_epi, r_b_epi, r_c_epi;

void load_parameters()
{

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("parameters.ini", pt);

    h = boost::lexical_cast<double>(pt.get<std::string>("general.h"));

    l_c_x = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.center_x"));
    l_c_y = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.center_y"));
    l_c_z = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.center_z"));

    l_a_endo = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.a_endo"));
    l_a_trunc = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.a_trunc"));
    l_b_endo = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.b_endo"));
    l_c_endo = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.c_endo"));

    l_t_equator_b = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.t_equator_b"));
    l_t_equator_c = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.t_equator_c"));
    l_t_apex = boost::lexical_cast<double>(pt.get<std::string>("left_ventricle.t_apex"));

    l_a_epi = l_a_endo + l_t_apex;
    l_b_epi = l_b_endo + l_t_equator_b;
    l_c_epi = l_c_endo + l_t_equator_c;

    r_c_x = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.center_x"));
    r_c_y = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.center_y"));
    r_c_z = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.center_z"));

    r_a_endo = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.a_endo"));
    r_a_trunc = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.a_trunc"));
    r_b_endo = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.b_endo"));
    r_c_endo = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.c_endo"));

    r_t_equator_b = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.t_equator_b"));
    r_t_equator_c = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.t_equator_c"));
    r_t_apex = boost::lexical_cast<double>(pt.get<std::string>("right_ventricle.t_apex"));

    r_a_epi = r_a_endo + r_t_apex;
    r_b_epi = r_b_endo + r_t_equator_b;
    r_c_epi = r_c_endo + r_t_equator_c;

}

// Function
FT ellipsoid (const Point& p,
	      double o_x, double o_y, double o_z,
	      double a, double b, double c)
{
    return sqrt(b*b*c*c*(p.z() - o_z)*(p.z() - o_z) +
		c*c*a*a*(p.x() - o_x)*(p.x() - o_x) +
		a*a*b*b*(p.y() - o_y)*(p.y() - o_y)) - a*b*c;
}

FT ventricles (const Point& p)
{
    if((ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_epi,  l_b_epi,  l_c_epi)  <= 0.0 &&
	ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_endo, l_b_endo, l_c_endo) >= 0.0 &&
	p.z() <= l_a_trunc) ||
       (ellipsoid(p, r_c_x, r_c_y, r_c_z, r_a_epi,  r_b_epi,  r_c_epi)  <= 0.0 &&
	ellipsoid(p, r_c_x, r_c_y, r_c_z, r_a_endo, r_b_endo, r_c_endo) >= 0.0 &&
	ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_endo, l_b_endo, l_c_endo) >  0.0 &&
	p.z() <= r_a_trunc))
	return -1;
    else
	return 1;
}

#include <cmath>

int main()
{
    //Load paramters from file
    load_parameters();

    // Domain
    Mesh_domain domain(ventricles, K::Sphere_3(Point(0, 0, 0), 9.));

    // Mesh criteria
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
