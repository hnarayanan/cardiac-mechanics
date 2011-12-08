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
r_a_endo = 1.0;
r_b_endo = 1.0;
r_a_trunc = 0.75;
r_t_equator = 0.25;
r_t_apex = 0.25;

////////////////////////////////////////////////////////////////////////

/* // Compute other useful dimensions for the left ventricle */
/* l_a_epi = l_a_endo + l_t_apex; */
/* l_b_epi = l_b_endo + l_t_equator; */

/* l_b_endo_trunc = l_b_endo*Sqrt(1.0 - l_a_trunc^2/l_a_endo^2); */
/* l_b_epi_trunc  = l_b_epi*Sqrt(1.0 - l_a_trunc^2/l_a_epi^2); */

/* // Define control points for the left ventricle */
/* Point(1)  = { l_b_endo,  0.0,      0.0, h}; */
/* Point(2)  = { 0.0     ,  l_a_endo, 0.0, h}; */
/* Point(3)  = {-l_b_endo,  0.0,      0.0, h}; */
/* Point(4)  = { 0.0     , -l_a_endo, 0.0, h}; */
/* Point(5)  = { l_b_epi, 0.0,     0.0, h}; */
/* Point(6)  = { 0.0,     l_a_epi, 0.0, h}; */
/* Point(7)  = {-l_b_epi, 0.0,     0.0, h}; */
/* Point(8)  = { 0.0   , -l_a_epi, 0.0, h}; */

/* Point(10) = { 0.0,     0.0,    0.0, h}; */

/* Point(11) = {l_b_endo_trunc, l_a_trunc, 0.0, h}; */
/* Point(12) = {l_b_epi_trunc, l_a_trunc, 0.0, h}; */

/* // Define control lines for the left ventricle */
/* Ellipse(1) = {4, 10, 2, 1}; */
/* Ellipse(2) = {1, 10, 3, 11}; */
/* Line(3)    = {11, 12}; */
/* Ellipse(4) = {12, 10, 8, 5}; */
/* Ellipse(5) = {5, 10, 7, 8}; */
/* Line(6)    = {8, 4}; */

/* Line Loop(9)  = {1, 2, 3, 4, 5, 6}; */

/* // Define the control surface for the left ventricle */
/* Plane Surface(10) = {9}; */

/* // Revolve the control surface to generate the volume of the left */
/* // ventricle */
/* vol01[] = Extrude{{0, 1, 0}, {0, 0, 0}, 2*Pi/3}{ Surface{10}; }; */
/* vol02[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol01[0]}; }; */
/* vol03[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol02[0]}; }; */

// Compute other useful dimensions for the right ventricle
r_a_epi = r_a_endo + r_t_apex;
r_b_epi = r_b_endo + r_t_equator;

r_b_endo_trunc = r_b_endo*Sqrt(1.0 - r_a_trunc^2/r_a_endo^2);
r_b_epi_trunc  = r_b_epi*Sqrt(1.0 - r_a_trunc^2/r_a_epi^2);

// Define control points for the right ventricle
Point(21)  = { r_b_endo,  0.0,      0.0, h};
Point(22)  = { 0.0     ,  r_a_endo, 0.0, h};
Point(23)  = {-r_b_endo,  0.0,      0.0, h};
Point(24)  = { 0.0     , -r_a_endo, 0.0, h};
Point(25)  = { r_b_epi, 0.0,     0.0, h};
Point(26)  = { 0.0,     r_a_epi, 0.0, h};
Point(27)  = {-r_b_epi, 0.0,     0.0, h};
Point(28)  = { 0.0   , -r_a_epi, 0.0, h};

Point(30) = { 0.0,     0.0,    0.0, h};

Point(31) = {r_b_endo_trunc, r_a_trunc, 0.0, h};
Point(32) = {r_b_epi_trunc, r_a_trunc, 0.0, h};

// Define control lines for the right ventricle
Ellipse(21) = {24, 30, 22, 21};
Ellipse(22) = {21, 30, 23, 31};
Line(23)    = {31, 32};
Ellipse(24) = {32, 30, 28, 25};
Ellipse(25) = {25, 30, 27, 28};
Line(26)    = {28, 24};

Line Loop(29)  = {21, 22, 23, 24, 25, 26};

// Define the control surface for the right ventricle
Plane Surface(30) = {29};

// Revolve the control surface to generate the volume of the left
// ventricle
vol04[] = Extrude{{0, 1, 0}, {0, 0, 0}, 2*Pi/3}{ Surface{30}; };
vol05[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol04[0]}; };
vol06[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol04[0]}; };
