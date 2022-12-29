function [Ke,Fe]=KeFe 
global k k1 k2 kx xe xb gaus gausw S dS f c b nnodes nint J;
ksi=gaus; %gauss quadrature points%
Fe=zeros(nnodes,1); %create an empty matrix for local load vector%
Ke=zeros(nnodes,nnodes); %create an empty matrix for local stiffness vector%
J=(xe-xb)/2; %Jacobian which is element length/2)
for j=1:nnodes
for i=1:nint
    x=(xe-xb)*0.5*ksi(i)+(xe+xb)*0.5; %local to global coordinate transfoormation for function evaluation%
    Fe(j,1)=Fe(j,1)+eval(f,x)*S(j,i).*J*gausw(i); %calculate local load vector%
end
end
for j=1:nnodes
    for l=1:nnodes
        for i=1:nint
            x=(xe-xb)*0.5*ksi(i)+(xe+xb)*0.5; % local to global coordinate change if k,b and c are function of x, we should evaluate them in global c%
            if x<kx
                k=k1;
            else
                k=k2;
            end
            Ke(j,l)=Ke(j,l)+(eval(k,x)*dS(j,i)/J*dS(l,i)/J+eval(c,x)*S(j,i)*dS(l,i)/J+...
eval(b,x)*S(j,i)*S(l,i))*J*gausw(i); %calculate local stiffness matrix%
        end 
    end
end
end