------cube case-----
---Step 1: generating the corresponding cell complex, in this case a cube--
needsPackage "SimplicialComplexes"
numfaces = 6;
positionfaces={{0,0,1},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,-1}};
T=ZZ[sigma_1..sigma_numfaces,Degrees=>positionfaces];
Octahedral=simplicialComplex monomialIdeal(sigma_1*sigma_6,sigma_2*sigma_4,sigma_3*sigma_5);
facesOfCube = (faces(Octahedral))#0;
edgesOfOct = (faces(Octahedral))#1;
verticesOfCube = (faces(Octahedral))#2;
E = for m in edgesOfOct list (indices m)+{1,1};
vert = unique flatten E;

--Step 2: generating the list of ideals--
r=1;
R=generateAmbientRing(E,QQ);
basisOnChart={{R_0,R_1,1},{1,R_2,R_3},{R_4,1,R_5},{-1,R_6,R_7},{R_8,-1,R_9},{R_10,R_11,-1}};

mathfraka = (f)->(2*f-1);--(n_0,n_1)=(3,3)
mathfrakb = (f)-> -1;
Lgtm = (uplus,vplus,uminus,vminus)->({uminus-mathfrakb(uplus)*vplus, vminus-uplus-vplus*mathfraka(uplus)});

endpoints = (indexface1,indexface2) -> (faces star(Octahedral,sigma_indexface1*sigma_indexface2))#2
theOtherEndPoint = (indexface1, indexface2, indexface3) -> (
    if member(sigma_indexface1*sigma_indexface2*sigma_indexface3,endpoints(indexface1,indexface2)) then 
    for x in support(sum endpoints(indexface1,indexface2)-sigma_indexface1*sigma_indexface2*sigma_indexface3) list index(x)+1)  

orientationOfFaces = (indexface1, indexface2, indexface3)->(det matrix{degree sigma_indexface1, degree sigma_indexface2, degree sigma_indexface3})
putIntoClockwise = (indexface1, indexface2, indexface3)->(
    if orientationOfFaces(indexface1,indexface2,indexface3) == 1 then {indexface2, indexface1, indexface3}
    else if orientationOfFaces(indexface1,indexface2,indexface3) == -1 then {indexface1,indexface2,indexface3}
    else "N/A"
    )
verticesOfFace = (indexface) -> (faces star (Octahedral,sigma_indexface))#2

coordPlus = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    coordp00=product for x in Lindexfaces list sigma_x;
    coordp10=product for x in theOtherEndPoint(Lindexfaces_0,Lindexfaces_1,Lindexfaces_2) list sigma_x;
    coordp01=product for x in theOtherEndPoint(Lindexfaces_0,Lindexfaces_2,Lindexfaces_1) list sigma_x;
    coordp11= (flatten delete(coordp00,delete(coordp01,delete(coordp10,verticesOfFace(Lindexfaces_0)))))_0;
    matrixcoordp =promote(matrix {degree coordp00, degree coordp01, degree coordp10, degree coordp11},QQ);
    coordup = solve(matrixcoordp,promote(matrix{{0},{0},{1},{1}},QQ));
    coordvp = solve(matrixcoordp,promote(matrix{{0},{1},{0},{1}},QQ));
    transpose matrix{flatten entries coordup, flatten entries coordvp}
    )
coordMinus = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    coordm00=product for x in Lindexfaces list sigma_x;
    coordm01=product for x in theOtherEndPoint(Lindexfaces_1,Lindexfaces_0,Lindexfaces_2) list sigma_x;
    coordm10=product for x in theOtherEndPoint(Lindexfaces_1,Lindexfaces_2,Lindexfaces_0) list sigma_x;
    coordm11= (flatten delete(coordm00,delete(coordm01,delete(coordm10,verticesOfFace(Lindexfaces_1)))))_0;
    matrixcoordm =promote(matrix {degree coordm00, degree coordm01, degree coordm10, degree coordm11},QQ);
    coordum = solve(matrixcoordm,promote(matrix{{0},{0},{1},{1}},QQ));
    coordvm = solve(matrixcoordm,promote(matrix{{0},{1},{0},{1}},QQ));
    transpose matrix{flatten entries coordum, flatten entries coordvm}
    )
affineTransformation = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    indexp = Lindexfaces_0;
    indexm = Lindexfaces_1;
    pmatrix = matrix {basisOnChart_(indexp-1)}*coordPlus(indexface1,indexface2,indexface3);
    mmatrix = matrix {basisOnChart_(indexm-1)}*coordMinus(indexface1,indexface2,indexface3);
    {flatten entries pmatrix, flatten entries mmatrix}--output {{up,vp},{um,vm}}
    )

generateItau = (indexface1,indexface2)->(
    endpointsL = endpoints(indexface1,indexface2)/(sigma_indexface1*sigma_indexface2);
    indexf3 = (index endpointsL_0)+1;
    indexf4 = (index endpointsL_1)+1;
    matrixOfuvAround123=matrix affineTransformation(indexface1,indexface2, indexf3);
    up123 = matrixOfuvAround123_(0,0);
    vp123 = matrixOfuvAround123_(0,1);
    um123 = matrixOfuvAround123_(1,0);
    vm123 = matrixOfuvAround123_(1,1);
    matrixOfuvAround124=matrix affineTransformation(indexface1,indexface2, indexf4);
    up124 = matrixOfuvAround124_(0,0);
    vp124 = matrixOfuvAround124_(0,1);
    um124 = matrixOfuvAround124_(1,0);
    vm124 = matrixOfuvAround124_(1,1);
    ideal Lgtm(up123,vp123,um123,vm123)+ideal Lgtm(up124,vp124,um124,vm124)+ideal(um123^(r+1),vp123^(r+1),um124^(r+1),vp124^(r+1))
    )
--hashTable edge => ideals
hTIdeals = hashTable for e in E list e=>generateItau(first e, last e);