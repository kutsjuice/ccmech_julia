// Gmsh project created on Fri May 10 13:59:29 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, -0.02, 0, 0.005};
//+
Point(2) = {0, 0.02, 0, 0.005};
//+
Point(3) = {0.04, 0.02, 0, 0.005};
//+
Point(4) = {0.24, 0.02, 0, 0.005};
//+
Point(5) = {0.24, -0.02, 0, 0.005};
//+
Point(6) = {0.04, -0.02, 0, 0.005};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 6};
//+
Line(3) = {6, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 6};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Physical Surface("fixed", 8) = {1};
//+
Physical Curve("fixed", 9) = {4, 1, 2, 3};
//+
Physical Point("fixed", 10) = {1, 2, 3, 6};
//+
Extrude {0, 0, 0.025} {
  Surface{1}; Surface{2}; 
}
//+
Physical Volume("body", 21) = {1, 2};
