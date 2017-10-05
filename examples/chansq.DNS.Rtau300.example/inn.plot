c: infile inn.plotU 1        U1 contours U2 U3 vectors  at x=0,2,4

c: if imagebarU > 0 find HAVEBARU continue
c: bar imagebarU 200 50 b
c: image outgif imagebarU out/Ubar end
HAVEBARU

c: infile inn.plotUxa 1       x/oct ave at current time
c: infile inn.plotUave 1      time/x ave then oct ave of that
c: infile inn.plotconv 1      convergence plots U2 Rd drho/dt (time)

