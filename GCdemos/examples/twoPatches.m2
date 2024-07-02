---- two patches example ----
E = {{1,2}};
R = generateAmbientRing(E,QQ)
r=1;
--mathfraka=(f) -> 2*f-1;--(n_0,n_1)=(3,3)
--mathfraka=(f)-> f^2; -- (n_0,n_1)=(4,3)
mathfraka=(f)->-f^2+2*f-1; --(n_0,n_1)=(3,4)
--mathfraka=(f)-> 0--(n_0,n_1)=(4,4)
gtm=(uminus,vminus,uplus,vplus)->{uminus+vplus, vminus-(uplus+vplus*mathfraka(uplus))};
J = ideal gtm(R_0,R_1,R_2,R_3);
I = J+ideal(R_0^(r+1),R_3^(r+1));