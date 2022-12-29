function [S,dS]=shape(deg)
global gaus nint
S=zeros(deg,nint); %deg by nint because each shape function evaluated at nint gauss quadrature points%
dS=zeros(deg,nint); %the derivatives of shape functions are also calculated at nint points%
ksi=gaus; %gauss quadrature points%
if deg==1 %linear shape functions%
    S(1,:)=(1-ksi)./2;
    S(2,:)=(1+ksi)./2;
    dS(1,:)=-0.5;
    dS(2,:)=0.5;
elseif deg==2 %quadratic shape functions%
    S(1,:)=-ksi.*(1-ksi)./2;
    S(2,:)=1-ksi.^2;
    S(3,:)=ksi.*(1+ksi)./2;
    dS(1,:)=ksi-0.5;
    dS(2,:)=-2*ksi;
    dS(3,:)=ksi+0.5;
elseif deg==3 %cubic shape functions%
    S(1,:)=-9/16.*(ksi+1/3).*(ksi-1/3).*(ksi-1);
    S(2,:)=27/16.*(ksi+1).*(ksi-1/3).*(ksi-1);
    S(3,:)=-27/16.*(ksi+1).*(ksi+1/3).*(ksi-1);
    S(4,:)=9/16.*(ksi+1).*(ksi+1/3).*(ksi-1/3);
    dS(1,:)=1/16.*(-27*ksi.^2+18*ksi+1);
    dS(2,:)=9/16.*(9*ksi.^2-2*ksi-3);
    dS(3,:)=-9/16.*(9*ksi.^2+2*ksi-3);
    dS(4,:)=1/16.*(27*ksi.^2+18*ksi-1);
end
end