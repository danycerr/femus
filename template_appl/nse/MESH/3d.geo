Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(5) = {0, 1, 0, 1};
Line(1) = {1, 5};
Line(2) = {3, 5};
Line(3) = {3, 2};
Line(4) = {1, 2};
Line Loop(6) = {1, -2, 3, -4};
Transfinite Line {1, -2, 3, -4}= 6 Using Bump 1;
Plane Surface(6) = {6};
Transfinite Surface {6};
Recombine Surface {6};
Extrude {0, 0, 1} {
  Surface{6};
}
Transfinite Line {2, 18, 9, 14,11,8}= 6 Using Bump 1;
Recombine Surface {19,28,15};


