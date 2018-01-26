//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 2, 2, 0};
//+
Line Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {4, 1, 2, 3};
//+
Surface(3) = {3};
//+
Physical Surface("plate") = {3};
