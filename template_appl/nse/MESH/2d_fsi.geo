cl__1 = 1;
nx =4;
ny = 6;
Point(1) = {0., 0., 0, 1};
Point(2) = {0.2, 0., 0, 1};
Point(3) = {0.3, 0., 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Transfinite Line {1, 2} = nx  Using Bump 1;
Extrude {0, 0.5, 0} {
  Line{1, 2};Layers{ny};Recombine;
}
Physical Surface("2") = {6};
Physical Line("10") = {1,2};
Physical Line("11") = {3,7};
Physical Line("12") = {4, 9};
Physical Line("1000") = {5};
Physical Surface("2") += {6};
Physical Surface("4") = {10};
