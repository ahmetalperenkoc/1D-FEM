function [K,F,Kr,Fr,u,ee,een,l2n,uh,uph,ue,upe]=buildKF
global k nel deg h f0 f0x alfa1 alfa2 beta1 beta2 gama1 gama2 beta0 xmin xmax nint gaus gausw S dS nnodes xe xb ntot;
setproblem();
nnodes=deg+1; %number of nodes per element%
ntot=nel*(nnodes-1)+1; %number of total nodes%
xb=xmin; %initial element start point%
if nint==1 %gauss quadrature points and corresponding gauss quadrature weights for numeric integration %
gaus=0;
gausw=2;    
elseif nint==2
gaus=[-sqrt(1/3) sqrt(1/3)];
gausw=[1 1];    
elseif nint==3
gaus=[-sqrt(3/5) 0 sqrt(3/5)];
gausw=[5/9 8/9 5/9]; 
elseif nint==4
gaus=[-sqrt((15+2*sqrt(30))/35) -sqrt((15-2*sqrt(30))/35) sqrt((15-2*sqrt(30))/35) sqrt((15+2*sqrt(30))/35)];
gausw=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36 ];    
end
[S,dS]=shape(deg);%call shape functions 
K=zeros(ntot); % create an empty stiffness matrix%
F=zeros(ntot,1);% create an empty load vector%
F0=zeros(ntot,1); %create an empty vector for point loads%
z=f0x/h(1)*deg+1; %point load node number%   
for i=1:ntot
    if i==z
F0(i)=f0;
    end
end
for j=1:nel
xe=xb+h(j); %update element end point%
    for l=1:nnodes
    nrow=(j-1)*(nnodes-1)+l;
    [Ke,Fe]=KeFe;
    F(nrow,1)=F(nrow,1)+sum(Fe(l,:)); %assembly of local load vector to global%
    for z=1:nnodes
    ncol=(j-1)*(nnodes-1)+z;
    K(nrow,ncol)=K(nrow,ncol)+Ke(l,z); %assembly of local stiffness matrix to global%
    end
    end
xb=xb+h(j); %update element starting point%
end
F=F+F0;
if (alfa1~=0) && beta1 ==0 && beta2==0 %if two neumann boundary conditions%
    Kr=K;
    Kr(1,:)=zeros(1,ntot);
    Kr(1,1)=1;
    Fr=F;
    Fr(1,1)=0;    
elseif (alfa1~=0) %if general boundary condition at xmin%
 x=xmin;
 Kr=K;
 Kr(1,1)=K(1,1)-eval(k)*beta1/alfa1;
 Fr=F;
 Fr(1,1)=F(1,1)-eval(k)*gama1/alfa1;
else %if esential boundary condition at xmin use penalty method%
x=xmin;
 Kr=K;
 Kr(1,1)=K(1,1)-eval(k)*beta0;
 Fr=F;
 Fr(1,1)=F(1,1)-eval(k)*beta0*gama1;
end
if (alfa1~=0) && beta1 ==0 && beta2==0 %if two neumann boundary conditions%
x=xmax;
Kr=Kr;
Fr(end,1)=F(end,1)+eval(k)*gama2/alfa2;
elseif (alfa2~=0) %if general boundary condition at xmin%
 x=xmax;
 Kr(ntot,ntot)=K(ntot,ntot)+eval(k)*beta2/alfa2;
 Fr(ntot,1)=F(ntot,1)+eval(k)*gama2/alfa2;
else %if essential boundary condition at xmax use penalty method%
x=xmax;
 Kr(ntot,ntot)=K(ntot,ntot)+eval(k)*beta0;
 Fr(ntot,1)=F(ntot,1)+eval(k)*beta0*gama2;
end
u=inv(Kr)*Fr;
[ee,een,l2n]=norm123(u);
[uh,uph,ue,upe]=interp(u);