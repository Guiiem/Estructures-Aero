function [Iz_B, S_B] = computeInertia_B(b, h, t)
% We get Iz_B and section_B
Iz_B = 2*((1/12)*b*t^3+b*t*(h/2)^2)+t*(h-t)^3*(1/12);
S_B = (h-t)*t+2*b*t;
end