
function sigcr = sigvinc(Nelements, x, Ic, Ib, mat, Tmat, Tn)

% Compute critical stress for buckling (PANDEO)
sigcr = zeros(Nelements,1); %sigma crÌtica
for e=1:Nelements
K=1; %beta dels apunts de pandeo
x1=x(Tn(e,1),1); x2=x(Tn(e,2),1);
y1=x(Tn(e,1),2); y2=x(Tn(e,2),2);
z1=x(Tn(e,1),3); z2=x(Tn(e,2),3);
le=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
if Tmat(e)==1
E=mat(Tmat(e),1);
Area=mat(Tmat(e),2);
I=Ic;
else
E=mat(Tmat(e),1);
Area=mat(Tmat(e),2);
I=Ib;
Area = 4.34e-7;
end
Pcr=(-pi*2*E*I)/(K*le)^2;
sigcr(e)=Pcr/Area;
end
end