function u = solveEDO (massa_bc, Cd, Rhoa,S, g)

syms y(t);
eqn = diff(y,t) == g - (1/(2*massa_bc) *Rhoa*Cd*S*y^2);
cond = y(0) == 0;
u = dsolve (eqn, cond);
end