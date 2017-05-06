// Gmsh project created on Wed Oct 12 11:55:22 2016
cl__1 = 1;
Point(1) = {0., 0., 0, 1};
Point(2) = {1., 0., 0, 1};
Point(3) = {1., 1., 0, 1};
Point(4) = {0., 1., 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Recombine Surface {6};
Transfinite Surface {6} Alternated;
Transfinite Line  {1, 2, 3, 4} = 3  Using Bump 1;
Recombine Surface {6};
Extrude {0, 0, 2} {
  Surface{6};Layers{3};Recombine;
}
// Recombine Surface {15,19,28};
Physical Volume("2") = {1};
Physical Surface("12") = {15, 19, 23, 27};
Physical Surface("10") = {6};
Physical Surface("11") = {28};
