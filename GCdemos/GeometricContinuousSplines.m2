newPackage(
    "GeometricContinuousSplines",
    Version => "0.1",
    Date => "4 May 2023",
    Headline => "This is a package for geometrically continuous splines.",
    Authors => {{ Name => "Beihui Yuan", Email => "by238@cornell.edu", HomePage => "https://sites.google.com/view/beihuiyuan/home"}},
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageExports => {"Cyclotomic"}
    )

export {
    "computeBaseField",
    "createStarVertexPatch",
    "exportXYZcoeff",
    "generateAmbientRing",
    "gSplineBasis",
    "interpolationAtMiddlePoints",--experimental
    "monomialBasisBiDegree",
    "monomialBasisT"
    }

-* Todo *-
-- Give G splines a Type (?)--
-- improve interpolation, which does not work very well--


-* Code section *-
---------------------------
computeBaseField = method()
---------------------------
---------------------------
--This method compute a field extension of QQ such that 
--the transition map can be defined over. 
---------------------------
--Inputs: 
---------------------------
--valences = list of positive integers.
---------------------------
--Outputs:
---------------------------
--A field of extension over QQ.
--usage: computeBaseField(valences)
---------------------------
computeBaseField(List):= Ring =>(valences) ->(
    n := lcm(valences);
    KK0 := cyclotomicField(n);
    ksi := KK0_0;
    e := symbol e;
    R := QQ[for i in valences list e_i];
    phi := map(KK0, R, matrix {for i in valences list ksi^(sub(n/i,ZZ))/2+ksi^(sub(n-n/i,ZZ))/2});
    KK := R/ker phi;
    KK
    )
---------------------------
createStarVertexPatch = method()
---------------------------
---------------------------
--This method compute a basis for G spline spaces over a star of a vertex, 
--assuming that transition maps are from symmetric gluing data. 
---------------------------
--Inputs: 
---------------------------
--valences = list of positive integers. May have repetition.
-------------input the valences at each BOUNDARY vertex clockwisely 
--deg = a non-negative integer, the degree bound for the spline space
---------------------------
--Outputs:
---------------------------
--A matrix. A basis for the spline space over base filed RR
--usage: createStarVertexPatch(valences, deg)
---------------------------
---------------------------
--Function dependence:
---------------------------
-- gSplineBasis
-- computeBaseField
-- generateAmbientRing
---------------------------
createStarVertexPatch(List, ZZ):= (valences, deg) ->(
    n := #valences;
    E := for i from 1 to n-1 list {i,i+1};
    E = {{n,1}}|E;
    uniValences := unique ({n}|valences);
    --uniValences = {n}|uniValences;
    posValences := for m in ({n}|valences) list position(uniValences,i->i==m);
    KK := computeBaseField(uniValences);
    S := generateAmbientRing(E,KK);
    SRR := generateAmbientRing(E,RR);
    gluea := (f, i)->(2*KK_(posValences_0)*(1-f)^2-2*KK_(posValences_i)*f^2); 
    glueb := -1;
    um := S_(2*n-2);
    vm := S_(2*n-1);
    up := S_0;
    vp := S_1;
    ideals := {ideal (um^2,vp^2,um+vp,vm-up-vp*gluea(up,n))};
    for i from 0 to n-2 do (
        um = S_(2*i);
        vm = S_(2*i+1);
        up = S_(2*i+2);
        vp = S_(2*i+3);
        ideals = ideals|{ideal(um^2,vp^2,um+vp,vm-up-vp*gluea(up,i+1))};);
    starPatchVertexBasis := (gSplineBasis(E,ideals,deg))#"basis";
    sub(starPatchVertexBasis, 
    (for i from 0 to #uniValences-1 list (KK_i => cos(2*pi/uniValences_i)))|
    for j from 0 to 2*n-1 list (S_j => SRR_j))
    )---complete on June 26, 2023.

---------------------------
exportXYZcoeff = method()
---------------------------
---------------------------
--This method creates .txt file containing coefficient matrices
---------------------------
--Inputs:
---------------------------
-- mP = a (3xn) matrix, 
--      giving the x,y,z coordinates, 
--      each entry is a polynomial.
-- uvrange = a list of two numbers.
-- fileName = a string, the name of output file.
---------------------------
--Outputs:
---------------------------
--A .txt file
--usage: createPyFile(mP, fileName)
---------------------------
exportXYZcoeff(Matrix, List, String) := (mP, uvrange, fileName) -> (
    S := ring mP;
    nump := numColumns mP;--number of patches
    f := concatenate{fileName, ".txt"} << ""; --initial a file
    polycoef := for k from 0 to (nump-1) list
      for polynf in entries mP_k list
      for i from 0 to degree(S_(2*k),polynf) list
      for j from 0 to degree(S_(2*k+1),polynf) list coefficient(S_(2*k)^i*S_(2*k+1)^j,polynf); 
    f << nump << endl;
    stringPolycoeff := replace("{","[",toString polycoef);
    stringPolycoeff = replace("}","]",stringPolycoeff);
    f << stringPolycoeff << endl;
    f << toString uvrange << endl << close;
    f
)

---------------------------
generateAmbientRing = method()
---------------------------
---------------------------
--This method creates an ambient ring for G spline space. 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
-- kk = base field
---------------------------
--Outputs:
---------------------------
--A polynomial ring kk[u_sigma,v_sigma:sigma in vert]
--usage: generateAmbientRing(E)
---------------------------
generateAmbientRing(List,Ring):= Ring => (E,kk) ->(
    vert := sort unique flatten E;
    u := symbol u;
    v := symbol v;
    S := kk[flatten for sigma in vert list {u_sigma,v_sigma}];
    S
    )

generateAmbientRing(List,InexactFieldFamily):= Ring => (E,kk) ->(
    vert := sort unique flatten E;
    u := symbol u;
    v := symbol v;
    S := kk[flatten for sigma in vert list {u_sigma,v_sigma}];
    S
    )

---------------------------
interpolationAtMiddlePoints = method()--not working well, needs tests
---------------------------
---------------------------
--This method computes subspaces of splines passing through given points 
---------------------------
--Inputs:
---------------------------
-- mBasis = a matrix of the basis of spline space, 
-- points3D = a matrix of points in 3D
---------------------------
--Outputs:
---------------------------
--A matrix of bases of splines passing through given points
--usage: interpolationAtMiddlePoints()
---------------------------
interpolationAtMiddlePoints(Matrix, Matrix) := (mBasis, points3D) -> (
    S:= ring mBasis_(0,0);
    varsS := flatten entries vars S;
    n := numRows mBasis;
    A := sub(mBasis, flatten for i from 0 to n-1 list {S_(2*i) => 0.0, S_(2*i+1) => 1.0});
    kerA := gens ker A;
    B := solve(A,points3D);
    (mBasis*kerA, mBasis*B)
)

---------------------------
monomialBasisBiDegree = method()
---------------------------
--This method creates a list of monomial basis (for bi-degree)
--for each coordinate ring with bi-degree no more than (d,d). 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--S = ambient ring
--deg = bi-degree
---------------------------
--Outputs:
---------------------------
--A hash table sigma => a list of monomials over vertex sigma
--usage: monomialBasisBiDegree(E,S,deg)
---------------------------
monomialBasisBiDegree(List,Ring,ZZ) := (E,S,deg) ->(
    vert := sort unique flatten E;
    monB:= hashTable for sigma in vert list 
    sigma => flatten entries monomials((S_(2*sigma-2)+1)^deg*(S_(2*sigma-1)+1)^deg);--flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    monB
    )

---------------------------
monomialBasisT = method()
---------------------------
--This method creates a list of monomial basis (for total degree)
--for each coordinate ring with total degree no more than d. 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--S = ambient ring
--deg = the total degree
---------------------------
--Outputs:
---------------------------
--A hash table sigma => a list of monomials over vertex sigma
--usage: monomialBasisT(E,S,deg)
---------------------------
monomialBasisT(List,Ring,ZZ) := (E,S,deg) ->(
    vert := sort unique flatten E;
    monB:= hashTable for sigma in vert list 
    sigma => flatten entries monomials((S_(2*sigma-2)+S_(2*sigma-1)+1)^deg);--flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    monB
    )

---------------------------
gSplineBasis = method()
---------------------------
---------------------------
--This method computes the geometrically continuous spline spaces 
--associated with a graph whose edges are labeled by ideals. 
--Also see generalizedSplines in AlgebraicSplines.m2 package.
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--ideals = list of ideals that labeled by the edges
--deg = an integer, which is the degree of the spline space
---------------------------
--Outputs:
---------------------------
--A Hash Table of dimension and basis in the given geometrically continuous spline spaces, 
--up to the given degree
--usage: gSplineBasis(E, L, d)#"dimension", gSplineBasis(E,L,d)#"basis"
---------------------------
gSplineBasis(List,List,ZZ) := HashTable => (E,ideals,deg) ->(
    vert := sort unique flatten E;
    S:= ring first ideals;
    --one has to input the underlying ring before using this function
    --make sure ideals all lie in the same ring
    ideals = apply(ideals, I->sub(I,S));
    --Possible future development: check if the generators of ideals are
    ------------------in variables of the two faces, using function support(f)
    --hashTable E=>ideals, label ideals with the correponding edge
    labelIdeals := hashTable for i from 0 to #E-1 list E#i=>ideals#i;
    --Boundary map from Edges to vertices, 
    --its rows are labeled by edges and columns are labeled by vertices
    boundaryEV:= matrix apply(E,
	e->apply(vert,
	    sigma->if(sigma===first e) then 1
	    else if(sigma===last e) then -1
	    else 0));
    boundaryEV = sub(boundaryEV,S);
    --Generate a monomial basis for each vertex sigma and degree no more than deg
    --monB is a hash table sigma => a list of monomial basis over sigma
    varsList := flatten entries vars S;
    monB:= hashTable for sigma in vert list 
    sigma => flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    --hashTable remainders, {e,sigma, m} => remainder
    hTrem:= hashTable flatten flatten for e in E list
    for sigma in vert list
    for m in monB#sigma list {e,sigma,m}=>
    (if (sigma===first e) then m%(labelIdeals#e)
    else if (sigma===last e) then -m%(labelIdeals#e)
    else 0);
    --generate a matrix corresponding to delta_2 from the hashTable hTEvm
    --of which rows labeled by edges, columns labeled by faces
    bdryM := matrix for e in E list {
	(coefficients matrix{flatten for sigma in vert list 
		for m in monB#sigma list hTrem#{e,sigma,m}})_1
	};
    --compute the kernel of bdryM, denoted by kerBdryM.
    --columns of which form a basis to the spline space.
    kerBdryM := ker bdryM;
    dimBasis := numColumns gens kerBdryM;
    sourceB := directSum(for sigma in vert list
	matrix {for m in monB#sigma list m});
    GSplineSpaceB := sourceB * gens kerBdryM;
    --output: hashTable of dimension and basis
    hTresults := new HashTable from { 
	"dimension" => dimBasis, 
	"basis" => GSplineSpaceB};
    hTresults
    --coding completes, yet to be tested. --Mar. 12, 2023
)


-* Documentation section *-
beginDocumentation()

doc ///
Key
  GeometricContinuousSplines
Headline
 A package for G splines
Description
  Text
   Still trying to figure this out.
  --Example
  --CannedExample
Acknowledgement
Contributors
References
Caveat
SeeAlso
Subnodes
///

doc ///
Key
 gSplineBasis
Headline
 Headline for gSplineBasis
Usage
 gSplineBasis(edges,ideals,deg)
Inputs
 edges:List
     list of edges
 ideals:List
     list of ideals
 deg:ZZ
     desired degree
Outputs
 LB:List
      a list of basis
--Consequences
--  Item
Description
  Text
      This method works for...
  --Example
  --CannedExample
  --Code
  --Pre
--ExampleFiles
--Contributors
--References
--Caveat
SeeAlso
///

-* Test section *-
-- test code and assertions here
-- may have as many TEST sections as needed
TEST /// -* The 2-patches case *-
E = {{1,2}};
R = generateAmbientRing(E,QQ);
r=1;
mathfraka=(f) -> 2*f-1;--(n_0,n_1)=(3,3)
mathfraka=(f)-> f^2; -- (n_0,n_1)=(4,3)
mathfraka=(f)->-f^2+2*f-1; --(n_0,n_1)=(3,4)
--mathfraka=(f)-> 0--(n_0,n_1)=(4,4)
gtm=(uminus,vminus,uplus,vplus)->{uminus+vplus, vminus-(uplus+vplus*mathfraka(uplus))};
J = ideal gtm(R_0,R_1,R_2,R_3);
I = J+ideal(R_0^(r+1),R_3^(r+1));
d = 4;
hTtwopatches = gSplineBasis(E,{I}, d);
monomialBasisT(E,R,3)
-- test exportXYZcoeff and createPyFile --
mBasis = hTtwopatches#"basis";
dimS = hTtwopatches#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));
///
end--

-* Development section *-
restart
path = append(path,"./")
debug needsPackage "GeometricContinuousSplines"
check "GeometricContinuousSplines"

uninstallPackage "GeometricContinuousSplines"
restart
installPackage "GeometricContinuousSplines"
viewHelp "GeometricContinuousSplines"