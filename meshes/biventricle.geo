////////////////////////////////////////////////////////////////////////

// Mesh density
h = 0.25;

// Parameters defining the left ventricular geometry
l_a_endo = 2.5;
l_b_endo = 1.1;
l_a_trunc = 1.4;
l_t_equator = 0.8;
l_t_apex = 0.4;

// Parameters defining the right ventricular geometry
r_a_endo = 2.1;
r_b_endo = 1.2;
r_a_trunc = 0.75;
r_t_equator = 0.25;
r_t_apex = 0.25;

// Parameters relating the two geometries
offset_x = -2.0;

////////////////////////////////////////////////////////////////////////

// Compute other useful dimensions for the left ventricle
l_a_epi = l_a_endo + l_t_apex;
l_b_epi = l_b_endo + l_t_equator;

l_b_endo_trunc = l_b_endo*Sqrt(1.0 - l_a_trunc^2/l_a_endo^2);
l_b_epi_trunc  = l_b_epi*Sqrt(1.0 - l_a_trunc^2/l_a_epi^2);

// Define control points for the left ventricle
Point(1)  = { l_b_endo,  0.0,      0.0, h};
Point(2)  = { 0.0     ,  l_a_endo, 0.0, h};
Point(3)  = {-l_b_endo,  0.0,      0.0, h};
Point(4)  = { 0.0     , -l_a_endo, 0.0, h};
Point(5)  = { l_b_epi, 0.0,     0.0, h};
Point(6)  = { 0.0,     l_a_epi, 0.0, h};
Point(7)  = {-l_b_epi, 0.0,     0.0, h};
Point(8)  = { 0.0   , -l_a_epi, 0.0, h};

Point(10) = { 0.0,     0.0,    0.0, h};

Point(11) = {l_b_endo_trunc, l_a_trunc, 0.0, h};
Point(12) = {l_b_epi_trunc, l_a_trunc, 0.0, h};

// Define control lines for the left ventricle
Ellipse(1) = {4, 10, 2, 1};
Ellipse(2) = {1, 10, 3, 11};
Line(3)    = {11, 12};
Ellipse(4) = {12, 10, 8, 5};
Ellipse(5) = {5, 10, 7, 8};
Line(6)    = {8, 4};

Line Loop(9)  = {1, 2, 3, 4, 5, 6};

// Define the control surface for the left ventricle
Plane Surface(10) = {9};

// Revolve the control surface to generate the volume of the left
// ventricle
vol01[] = Extrude{{0, 1, 0}, {0, 0, 0}, 2*Pi/3}{ Surface{10}; };
vol02[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol01[0]}; };
vol03[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol02[0]}; };

// Compute other useful dimensions for the right ventricle
r_a_epi = r_a_endo + r_t_apex;
r_b_epi = r_b_endo + r_t_equator;

r_b_endo_trunc = r_b_endo*Sqrt(1.0 - r_a_trunc^2/r_a_endo^2);
r_b_epi_trunc  = r_b_epi*Sqrt(1.0 - r_a_trunc^2/r_a_epi^2);

offset_y = 0.0 + l_a_trunc - r_a_trunc;

// Define control points for the right ventricle
Point(101)  = { r_b_endo,  0.0,      0.0, h};
Point(102)  = { 0.0     ,  r_a_endo, 0.0, h};
Point(103)  = {-r_b_endo,  0.0,      0.0, h};
Point(104)  = { 0.0     , -r_a_endo, 0.0, h};
Point(105)  = { r_b_epi, 0.0,     0.0, h};
Point(106)  = { 0.0,     r_a_epi, 0.0, h};
Point(107)  = {-r_b_epi, 0.0,     0.0, h};
Point(108)  = { 0.0   , -r_a_epi, 0.0, h};

Point(110) = { 0.0,     0.0,    0.0, h};

Point(111) = {r_b_endo_trunc, r_a_trunc, 0.0, h};
Point(112) = {r_b_epi_trunc, r_a_trunc, 0.0, h};

// Define control lines for the right ventricle
Ellipse(101) = {104, 110, 102, 101};
Ellipse(102) = {101, 110, 103, 111};
Line(103)    = {111, 112};
Ellipse(104) = {112, 110, 108, 105};
Ellipse(105) = {105, 110, 107, 108};
Line(106)    = {108, 104};

Line Loop(109) = {101, 102, 103, 104, 105, 106};

// Define the control surface for the right ventricle
Plane Surface(110) = {109};

// Translate the surface of the right ventricle by the offset and
// rotate it initially
Translate { offset_x, offset_y, 0 } { Surface{110}; }
Rotate{{0, 1, 0}, {0.0 + offset_x, 0.0  + offset_y, 0}, Pi/3}{ Surface{110}; }

// Revolve the control surface to generate the volume of the left
// ventricle
vol04[] = Extrude{{0, 1, 0}, {0.0 + offset_x, 0.0  + offset_y, 0}, 2*Pi/3}{ Surface{110}; };
vol05[] = Extrude{{0, 1 ,0}, {0.0 + offset_x, 0.0  + offset_y, 0}, 2*Pi/3}{ Surface{vol04[0]}; };
