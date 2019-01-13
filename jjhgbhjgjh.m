js = 1366;
s = 5.67e-8;
syms T1 T2

alfa1 = 0.3;
rho1 = 0.3;
t2 = 0.4;
eqn1 = js*alfa1-2*s*T1^4+s*T2^4;
eqn2 = t2*js-s*T2^4+s*T1^4;

sol = solve([eqn1, eqn2],[T1,T2]);
sol.T1
sol.T2
temp1 = sol.T1(2);
temp2 = round(sol.T2)




