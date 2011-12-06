lc = 0.1;

a_endo = 2.5;
b_endo = 1.0;
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


Ellipse(1) = {1, 10, 3, 2};
Ellipse(2) = {2, 10, 4, 3};
Ellipse(3) = {3, 10, 1, 4};
Ellipse(4) = {4, 10, 2, 1};

Ellipse(5) = {5, 10, 7, 6};
Ellipse(6) = {6, 10, 8, 7};
Ellipse(7) = {7, 10, 5, 8};
Ellipse(8) = {8, 10, 6, 5};

Line Loop(9)  = {1, 2, 3, 4};
Line Loop(10) = {5, 6, 7, 8};

Plane Surface(11) = {9,10};
