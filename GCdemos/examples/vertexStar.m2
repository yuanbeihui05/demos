---- A star of a vertex (3 faces) example ---
E = {{1,2},{2,3},{3,1}};
R = generateAmbientRing(E,QQ);
r=1;
mathfraka = (f)->-f^2+2*f-1; 
mathfrakb = (f) -> -1;
Lgtm = (uplus,vplus,uminus,vminus)->({uminus-mathfrakb(uplus)*vplus, vminus-uplus-vplus*mathfraka(uplus)});
generateItau = (indexface1,indexface2)->(
    ideal Lgtm(R_(indexface1*2-2),R_(indexface1*2-1),R_(indexface2*2-2),R_(indexface2*2-1))+ideal(R_(indexface1*2-1)^(r+1),R_(indexface2*2-2)^(r+1))
    );
hTIdeals = hashTable for e in E list e=>generateItau(first e, last e);