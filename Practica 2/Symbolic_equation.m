%%How to use symbolic equations

syms u t ;
u = 3*t^5;
t=2;
vpa(subs(u)); %vpa approximates the result instead of returning the exact symbolic equation
