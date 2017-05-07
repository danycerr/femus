// Gmsh project created on Mon Oct 17 11:43:57 2016
Point(1) = {0, 0.0, 0, 0.0};
Point(2) = {1, 0, 0, 0.0};
Point(3) = {1, 1.0, 0, 0.0};
Point(4) = {0, 1.0, 0, 0.0};
Point(5) = {0.4, 0.0, 0, 0.0};
Point(6) = {0.6, 0.0, 0, 0.0};
Point(7) = {0.4, 0.6, 0, 0.0};
Point(8) = {0.6, 0.6, 0, 0.0};

Line(1) = {1, 5};
Line(2) = {5, 7};
Line(3) = {7, 8};
Line(4) = {8, 6};
Line(5) = {6, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line(9) = {5, 6};
Line(10) = {6, 8};
Line(11) = {8, 7};
Line(12) = {7, 5};
Line Loop(13) = {1, 2, 3, 4,5,6,7,8};
Line Loop(14) = {9,10,11,12};
Plane Surface(15) = {13};
Plane Surface(16) = {14};
Transfinite Line "*" = 20 Using Bump 1;
Transfinite Line {9,3} = 4 Using Bump 1;
Transfinite Line {2,4} = 12 Using Bump 1;
// Transfinite Surface {16} Alternated;
// Recombine Surface {15, 16};
Physical Line("12") = {1, 5, 7};
Physical Line("22") = {9};
Physical Line("1000") = {2, 3, 4};
Physical Line("10") = {8};
Physical Line("11") = {6};

Physical Surface("4") = {16};
Physical Surface("2") = {15};
