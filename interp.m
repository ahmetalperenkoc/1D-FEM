function [uh,uph,ue,upe]=interp(u)
global deg nel 
Ksi=linspace(-1,1,11);
uph=zeros(nel,length(Ksi));
uh=zeros(nel,length(Ksi));
kk=zeros(nel,length(Ksi));
xb=0;
for i=1:nel
    xe(i)=xb+1/(nel); 
    xb=xb+1/(nel);
end
xb=[0 xe(1:nel)];
xx=zeros(nel,length(Ksi));
for i=1:nel
    for j=1:length(Ksi)
    ksi=Ksi(j);
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
        for k=1:deg+1
            uph(i,j)=uph(i,j)+dS(k)*u((i-1)*deg+k)*(2/(xe(1)-xb(1))); %interpolate derivative of u between nodes%
            uh(i,j)=uh(i,j)+S(k)*u((i-1)*deg+k);    %interpolate u between nodes%
        end
    end
end
for i=1:nel
xx(i,:)=(xe(i)-xb(i)).*0.5.*Ksi+(xe(i)+xb(i))*0.5; %calculate gauss quadrature points in global coordinate system%
kk(i,:)=1./5+5*(xx(i,:)-0.5).^2; %calculate k at gauss quadrature points%
upe(i,:)=-atan(5*(xx(i,:)-0.5))-atan(2.5)+((1-xx(i,:)).*5./(1+25*(xx(i,:)-0.5).^2)); %calculate exact solution of secondary variable at gauss quadrature points%
ue(i,:)=(1-xx(i,:)).*(atan(5*(xx(i,:)-0.5))+atan(2.5)); %calculate exact solution of primary variable at gauss quadrature points%
end
ud=abs(ue-uh); %absolute different of primary variable %
upd=abs(upe-uph); %absolute different of secondary variable %
%the below codes are not output of the function I used them to create my plot%
[a b]=min(ud,[],2); %index number of primary variable error in each element and its magnitude%
[c d]=min(upd,[],2); %index number of secondary variable error in each element and its magnitude%
b=-1+0.2*(b-1); %convert matrix index to Ksi coordinate%
d=-1+0.2*(d-1); %convert matrix index to Ksi coordinate%
