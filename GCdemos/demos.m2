-----Demos for the paper: An algebraic framework for -----------
------------------------ geometrically continuous splines ------
-----------------
--- 2 patches ---
-----------------
restart
load "./GeometricContinuousSplines.m2"
load "./examples/twoPatches.m2"
d = 4;
hTtwopatches = gSplineBasis(E,{I}, d);
mBasis = hTtwopatches#"basis";
dimS = hTtwopatches#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));

-------------------------
--- a star of a vertex --
-------------------------
restart
load "./GeometricContinuousSplines.m2"
load "./examples/vertexStar.m2"
d = 4;
hTthreeSplit = gSplineBasis(E, for e in E list hTIdeals#e, d);
mBasis = hTthreeSplit#"basis";
dimS = hTthreeSplit#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));

------------------------
--- the cube example ---
------------------------
restart
load "./GeometricContinuousSplines.m2"
load "./examples/cube.m2"
hTIdeals; -- Hash table of ideals associated to each edge
d = 4;
hTcube = gSplineBasis(E, for e in E list hTIdeals#e, d);
mBasis = hTcube#"basis";
dimS = hTcube#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));
--interpolation
mP = mBasis*inverse sub(mBasis, for i from 0 to 11 list R_i=>0)*matrix{{0,0,1/1},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,-1}};
mP = transpose mP;
