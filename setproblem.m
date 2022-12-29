function [ ]= setproblem
% - d/dx [k(x) du/dx] + c(x) du/dx + b(x) * u = f(x) %
%boundary conditions alfa1 du/dx + beta1 u = gama1 alfa2 du/dx + beta2 u = gama2 % 
%input page for coefficients and boundary conditions%
global k1 k2 kx c b f f0 f0x nel deg h alfa1 alfa2 beta1 beta2 gama1 gama2 beta0 xmin xmax nint;
k1='-50'; %first k term%
k2='-50'; %second k term%
kx=0; %the location of change of k occures%
c='0';
b='0';
f='-100';
f0=0; %point load magnitude%
f0x=2;  %point load location in x coordinate%
alfa1=0;
alfa2=-5;
beta1=1;
beta2=0;
gama1=0;
gama2=15;
beta0=10000000; %penalty number%
xmin=2; %domain starting point%
xmax=8; %domain ending point%
nint=4;%number of gauss quadrature points range: 1-4%
deg=1; %degree of shape function 1=linear 2=quadratic 3=cubic%
nel=3; %number of elements%
h=ones(1,nel)*(xmax-xmin)/(nel); % equally spaced%
%if  not equally spaced elements enter each element length%
