cl__1 = 1;
Point(1) = {0., 0., 0, 1};
Point(2) = {1., 0., 0, 1};
Point(3) = {1., 2., 0, 1};
Point(4) = {0., 2., 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Transfinite Surface {6} Alternated;
Transfinite Line {3, 1} = 3  Using Bump 1;
Transfinite Line {4, 2} = 6  Using Bump 1;
Recombine Surface {6};
Physical Surface("2") = {6};
Physical Line("10") = {1};
Physical Line("12") = {4, 2};
Physical Line("11") = {3};
