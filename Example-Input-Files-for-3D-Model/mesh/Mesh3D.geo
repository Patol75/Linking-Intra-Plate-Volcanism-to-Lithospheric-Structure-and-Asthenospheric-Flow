height = 1e6;
width = 4e6;
thick = 4e6;
LAB = 9.1e5;
lowMtl = 3.4e5;
/*
ptcCenterX = 2.6e6;
ptcCenterY = 2e6;
ptcCenterZ = 0;
ptcOuterRadius = 2e5;
ptcInnerRadius = 5e4;
*/

// Inner Patch
/*
Point(1) = {ptcCenterX, ptcCenterY, ptcCenterZ};
Point(2) = {ptcCenterX + ptcInnerRadius, ptcCenterY, ptcCenterZ};
Point(3) = {ptcCenterX, ptcCenterY + ptcInnerRadius, ptcCenterZ};
Point(4) = {ptcCenterX - ptcInnerRadius, ptcCenterY, ptcCenterZ};
Point(5) = {ptcCenterX, ptcCenterY - ptcInnerRadius, ptcCenterZ};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
For i In { 1 : 4 }
    Physical Line(i) = {i};
EndFor
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
// Outer Patch
Point(6) = {ptcCenterX + ptcOuterRadius, ptcCenterY, ptcCenterZ};
Point(7) = {ptcCenterX, ptcCenterY + ptcOuterRadius, ptcCenterZ};
Point(8) = {ptcCenterX - ptcOuterRadius, ptcCenterY, ptcCenterZ};
Point(9) = {ptcCenterX, ptcCenterY - ptcOuterRadius, ptcCenterZ};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
For i In { 5 : 8 }
    Physical Line(i) = {i};
EndFor
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2, 1};
Physical Surface(2) = {2};
*/
// Front Face (Y == 0)
Point(10) = {0, 0, height};
Point(11) = {width, 0, height};
Point(12) = {width, 0, LAB};
Point(13) = {width, 0, lowMtl};
Point(14) = {width, 0, 0};
Point(15) = {0, 0, 0};
Point(16) = {0, 0, lowMtl};
Point(17) = {0, 0, LAB};
Line(9) = {10, 11};
Line(10) = {17, 12};
Line(11) = {16, 13};
Line(12) = {15, 14};
Line(13) = {11, 12};
Line(14) = {12, 13};
Line(15) = {13, 14};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 10};
For i In { 9 : 18 }
    Physical Line(i) = {i};
EndFor
Line Loop(3) = {9, 13, -10, 18};
Line Loop(4) = {10, 14, -11, 17};
Line Loop(5) = {11, 15, -12, 16};
For i In { 3 : 5 }
    Plane Surface(i) = {i};
    Physical Surface(i) = {i};
EndFor
// Back Face (Y == thick)
Point(18) = {0, thick, height};
Point(19) = {width, thick, height};
Point(20) = {width, thick, LAB};
Point(21) = {width, thick, lowMtl};
Point(22) = {width, thick, 0};
Point(23) = {0, thick, 0};
Point(24) = {0, thick, lowMtl};
Point(25) = {0, thick, LAB};
Line(19) = {18, 19};
Line(20) = {25, 20};
Line(21) = {24, 21};
Line(22) = {23, 22};
Line(23) = {19, 20};
Line(24) = {20, 21};
Line(25) = {21, 22};
Line(27) = {23, 24};
Line(28) = {24, 25};
Line(29) = {25, 18};
For i In { 19 : 29 }
    Physical Line(i) = {i};
EndFor
Line Loop(6) = {19, 23, -20, 29};
Line Loop(7) = {20, 24, -21, 28};
Line Loop(8) = {21, 25, -22, 27};
For i In { 6 : 8 }
    Plane Surface(i) = {i};
    Physical Surface(i) = {i};
EndFor
// Right Face (X == width)
Line(30) = {11, 19};
Line(31) = {12, 20};
Line(32) = {13, 21};
Line(33) = {14, 22};
For i In { 30 : 33 }
    Physical Line(i) = {i};
EndFor
Line Loop(9) = {30, 23, -31, -13};
Line Loop(10) = {31, 24, -32, -14};
Line Loop(11) = {32, 25, -33, -15};
For i In { 9 : 11 }
    Plane Surface(i) = {i};
    Physical Surface(i) = {i};
EndFor
// Left Face (X == 0)
Line(34) = {10, 18};
Line(35) = {17, 25};
Line(36) = {16, 24};
Line(37) = {15, 23};
For i In { 34 : 37 }
    Physical Line(i) = {i};
EndFor
Line Loop(12) = {34, -29, -35, 18};
Line Loop(13) = {35, -28, -36, 17};
Line Loop(14) = {36, -27, -37, 16};
For i In { 12 : 14 }
    Plane Surface(i) = {i};
    Physical Surface(i) = {i};
EndFor
// Top Face (Z == height)
Line Loop(15) = {9, 30, -19, -34};
Plane Surface(15) = {15};
Physical Surface(15) = {15};
// Bottom Face (Z == 0)
Line Loop(16) = {12, 33, -22, -37};

// Plane Surface(16) = {16, 2};
Plane Surface(16) = {16};
Physical Surface(16) = {16};

// Whole Domain

// Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
Surface Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
Volume(1) = {1};
Physical Volume(1) = {1};

Field[1] = Box;
Field[1].VIn = 15e3;
Field[1].VOut = 2e5;
Field[1].XMin = 0;
Field[1].XMax = width;
Field[1].YMin = 0;
Field[1].YMax = thick;
Field[1].ZMin = 8e5;
Field[1].ZMax = height;
/*
Field[2] = Box;
Field[2].VIn = 15e3;
Field[2].VOut = 2e5;
Field[2].XMin = ptcCenterX - ptcOuterRadius;
Field[2].XMax = ptcCenterX + ptcOuterRadius;
Field[2].YMin = ptcCenterY - ptcOuterRadius;
Field[2].YMax = ptcCenterY + ptcOuterRadius;
Field[2].ZMin = 0;
Field[2].ZMax = 5e4;

Field[3] = Min;
Field[3].FieldsList = {1, 2};

Background Field = 3;
*/
Background Field = 1;
