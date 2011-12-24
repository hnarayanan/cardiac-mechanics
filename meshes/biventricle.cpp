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
int angle;

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
    angle = boost::lexical_cast<int>(pt.get<std::string>("general.angle"));

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
	      double c_x, double c_y, double c_z,
	      double a, double b, double c)
{
    return sqrt(b*b*c*c*(p.z() - c_z)*(p.z() - c_z) +
		c*c*a*a*(p.x() - c_x)*(p.x() - c_x) +
		a*a*b*b*(p.y() - c_y)*(p.y() - c_y)) - a*b*c;
}

FT ventricles (const Point& p)
{
    if((ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_epi,  l_b_epi,  l_c_epi)  <= 0.0 &&
	ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_endo, l_b_endo, l_c_endo) >= 0.0 &&
	p.z() <= l_a_trunc) ||
       (ellipsoid(p, r_c_x, r_c_y, r_c_z, r_a_epi,  r_b_epi,  r_c_epi)  <= 0.0 &&
	ellipsoid(p, r_c_x, r_c_y, r_c_z, r_a_endo, r_b_endo, r_c_endo) >= 0.0 &&
	ellipsoid(p, l_c_x, l_c_y, l_c_z, l_a_endo, l_b_endo, l_c_endo) >  0.0 &&
	p.z() <= r_a_trunc + r_c_z))
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
    Mesh_domain domain(ventricles, K::Sphere_3(Point(0, 0, 0), 10.0));

    // Mesh criteria
    Mesh_criteria criteria(edge_size = h,
    			   facet_angle = 25, facet_size = 2*h,
    			   cell_radius_edge_ratio = 2, cell_size = 2*h);

    // Create edge that we want to preserve
    Polylines l_base_endo (1);
    Polyline& l_base_endo_segment = l_base_endo.front();
    Polylines l_base_epi (1);
    Polyline& l_base_epi_segment = l_base_epi.front();

    Polylines r_base_endo (1);
    Polyline& r_base_endo_segment = r_base_endo.front();
    Polylines r_base_epi (1);
    Polyline& r_base_epi_segment = r_base_epi.front();

    double r_l_endo = l_b_endo*std::sqrt(1.0 - l_a_trunc*l_a_trunc/(l_a_endo*l_a_endo));
    double r_l_epi  = l_b_epi*std::sqrt(1.0 - l_a_trunc*l_a_trunc/(l_a_epi*l_a_epi));

    double r_r_endo = r_b_endo*std::sqrt(1.0 - r_a_trunc*r_a_trunc/(r_a_endo*r_a_endo));
    double r_r_epi  = r_b_epi*std::sqrt(1.0 - r_a_trunc*r_a_trunc/(r_a_epi*r_a_epi));

    for(int i = 0; i < 360; i += 5)
    {
	Point l_p_endo (l_c_x + r_l_endo*std::cos(i*CGAL_PI/180),
			l_c_y + r_l_endo*std::sin(i*CGAL_PI/180),
			l_c_z + l_a_trunc);
	l_base_endo_segment.push_back(l_p_endo);
	Point l_p_epi (l_c_x + r_l_epi*std::cos(i*CGAL_PI/180),
		       l_c_y + r_l_epi*std::sin(i*CGAL_PI/180),
		       l_c_z + l_a_trunc);
	l_base_epi_segment.push_back(l_p_epi);

	if (i >= 180 - angle/2 && i <= 180 + angle/2)
	{
	    Point r_p_endo (r_c_x + r_r_endo*std::cos(i*CGAL_PI/180),
			    r_c_y + r_r_endo*std::sin(i*CGAL_PI/180),
			    r_c_z + r_a_trunc);
	    r_base_endo_segment.push_back(r_p_endo);
	}

	if (i >= 180 - angle && i <= 180 + angle)
	{
	    Point r_p_epi (r_c_x + r_r_epi*std::cos(i*CGAL_PI/180),
			   r_c_y + r_r_epi*std::sin(i*CGAL_PI/180),
			   r_c_z + r_a_trunc);
	    r_base_epi_segment.push_back(r_p_epi);
	}
    }
    l_base_endo_segment.push_back(l_base_endo_segment.front());
    l_base_epi_segment.push_back(l_base_epi_segment.front());
    r_base_endo_segment.push_back(r_base_endo_segment.front());
    r_base_epi_segment.push_back(r_base_epi_segment.front());

    // Insert edge in domain
    domain.add_features(l_base_endo.begin(), l_base_endo.end());
    domain.add_features(l_base_epi.begin(), l_base_epi.end());
    domain.add_features(r_base_endo.begin(), r_base_endo.end());
    domain.add_features(r_base_epi.begin(), r_base_epi.end());

    // Mesh generation without feature preservation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    std::ofstream medit_file("biventricle.mesh");
    c3t3.output_to_medit(medit_file);
    medit_file.close();
    c3t3.clear();

    return 0;
}
