// Cube for ImpedancePredictionVIE

// lx ly and lz via gmsh command line
//lx = 1.0; // length x direction
//ly = 1.0; // length y direction
//lz = 1.0; // length z direction



// #### 1st cube ##########################################

// Upper points
Point(1) = {+lx/2, +ly/2, 0.0};
Point(2) = {-lx/2, +ly/2, 0.0};
Point(3) = {-lx/2, -ly/2, 0.0};
Point(4) = {+lx/2, -ly/2, 0.0};

// Lower points
Point(5) = {+lx/2, -ly/2, -lz/2};
Point(6) = {-lx/2, -ly/2, -lz/2};
Point(7) = {-lx/2, +ly/2, -lz/2};
Point(8) = {+lx/2, +ly/2, -lz/2};


// Top
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Top2Bottom???
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Bottom
Line(9) = {4, 5};
Line(10) = {1, 8};
Line(11) = {2, 7};
Line(12) = {3, 6};


// 6 cube faces
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-4, 9, -8, -10};
Curve Loop(3) = {12, -5, -9, -3};
Curve Loop(4) = {-7, -11, -1, 10};
Curve Loop(5) = {-6, -12, -2, 11};
Curve Loop(6) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Cube volume 
Surface Loop(1) = {4, 6, 3, 5, 1, 2};

Volume(1) = {1};







// #### 2nd cube ##########################################

// Upper points
Point(9) = {+lx/2, +ly/2, lz/2};
Point(10) = {-lx/2, +ly/2, lz/2};
Point(11) = {-lx/2, -ly/2, lz/2};
Point(12) = {+lx/2, -ly/2, lz/2};


// Top
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

// Top2Bottom
Line(17) = {9, 1};
Line(18) = {10, 2};
Line(19) = {11, 3};
Line(20) = {12, 4};


// 6 cube faces
// Curve Loop(1) = {1, 2, 3, 4}; see above....
Curve Loop(7) = {1, -18, -13, 17};
Curve Loop(8) = {-2, -18, 14, 19};
Curve Loop(9) = {-3, -19, 15, 20};
Curve Loop(10) = {4, -17, -16, 20};
Curve Loop(11) = {13, 14, 15, 16};

// Plane Surface(1) = {1};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};

// Cube volume 
Surface Loop(2) = {1, 7, 8, 9, 10, 11};

Volume(2) = {2};






// #### "Physical" ##########################################


Physical Surface("TopElectrode") = {11};
Physical Surface("BottomElectrode") = {6};
Physical Surface("InsulatingContact") = {2,3,4,5,7,8,9,10};

Physical Volume("BodyVolume") = {1,2};