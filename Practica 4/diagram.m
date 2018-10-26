function [N, Ty, Mz] = diagram(f, Nelements, NnodesXelement, Td, R, theta)

Ty = zeros(Nelements, NnodesXelement);
N =  zeros(Nelements, NnodesXelement);
Mz =  zeros(Nelements, NnodesXelement);
R = R';

for i = 1:Nelements-1
    Ty(i,1) = f(Td(i,2))+R(Td(i,2));
    Ty(i,2) = f(Td(i,5))+R(Td(i,5));
    N(i,1) = f(Td(i,1))+R(Td(i,1));
    N(i,2) = f(Td(i,4))+R(Td(i,4));
    %Mz(i,1) = f(Td(i,3))+R(Td(i,3));
    %Mz(i,2) = Mz(i,1)+f(Td(i,6)) + R(Td(i,6));
end

i = 3;
    Ty(i,1) = f(Td(i,2))*cos(theta)+f(Td(i,1))*sin(theta)+R(Td(i,2))*cos(theta)+R(Td(i,5))*sin(theta);
    Ty(i,2) = f(Td(i,1))*cos(theta)+f(Td(i,4))*sin(theta)+R(Td(i,5))*cos(theta)+R(Td(i,2))*sin(theta)+Ty(i,1);
    N(i,1) = f(Td(i,1))*cos(theta)+f(Td(i,2))*sin(theta)+R(Td(i,1))*cos(theta)+R(Td(i,4))*sin(theta);
    N(i,2) = f(Td(i,4))+R(Td(i,4))*cos(theta)+R(Td(i,1))*sin(theta)+N(i,1);
    %Mz(i,1) = f(Td(i,3))+R(Td(i,3));
    %Mz(i,2) = Mz(i,1)+f(Td(i,6)) + R(Td(i,6));

