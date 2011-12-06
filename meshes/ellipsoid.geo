lc = 0.5;

a_endo = 2.5;
b_endo = 1.1;
t_equator = 1.0;
t_apex = 0.5;

a_epi = a_endo + t_apex;
b_epi = b_endo + t_equator;


Point(1)  = { b_endo,  0.0,    0.0, lc};
Point(2)  = { 0.0   ,  a_endo, 0.0, lc};
Point(3)  = {-b_endo,  0.0,    0.0, lc};
Point(4)  = { 0.0   , -a_endo, 0.0, lc};

Point(5)  = { b_epi,   0.0,    0.0, lc};
Point(6)  = { 0.0   ,  a_epi,  0.0, lc};
Point(7)  = {-b_epi,   0.0,    0.0, lc};
Point(8)  = { 0.0   , -a_epi,  0.0, lc};

Point(10) = { 0.0,     0.0,    0.0, lc};


Ellipse(1) = {4, 10, 2, 1};
Ellipse(2) = {1, 10, 3, 2};
Line(3)    = {2, 6};
Ellipse(4) = {6, 10, 8, 5};
Ellipse(5) = {5, 10, 7, 8};
Line(6)    = {8, 4};

Line Loop(9)  = {1, 2, 3, 4, 5, 6};
Plane Surface(10) = {9};

vol01[] = Extrude{{0, 1, 0}, {0, 0, 0}, 2*Pi/3}{ Surface{10}; };
vol02[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol01[0]}; };
vol03[] = Extrude{{0, 1 ,0}, {0, 0, 0}, 2*Pi/3}{ Surface{vol02[0]}; };
