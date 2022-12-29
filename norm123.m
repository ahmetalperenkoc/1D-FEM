function [ee,een,l2n]=norm123(u)
global deg nel
Ksi=[-0.9324695142 -0.6612093865 -0.2386191861 0.2386191861 0.6612093865 0.9324695142]; %6-point GQ%
w=[0.1713244924 0.3607615730 0.4679139346 0.4679139346 0.3607615730 0.1713244924]; %GQ points%
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
upe(i,:)=-atan(5*(xx(i,:)-0.5))-atan(2.5)+((1-xx(i,:)).*5./(1+25*(xx(i,:)-0.5).^2)); %calculte exact solution of secondary variable at gauss quadrature points%
ue(i,:)=(1-xx(i,:)).*(atan(5*(xx(i,:)-0.5))+atan(2.5)); %calculte exact solution of primary variable at gauss quadrature points%
end
ee=zeros(nel,1);
een=zeros(nel,1);
l2n=zeros(nel,1);
eene=zeros(nel,1);
l2ne=zeros(nel,1);
for i=1:nel
    for j=1:length(Ksi)
        ee(i)=ee(i)+0.5*kk(i,j)*(upe(i,j)-uph(i,j))^2*(xe(1)-xb(1))/2*w(j); %calculate energy norm with GQ&
        een(i)=een(i)+(kk(i,j)*(upe(i,j)-uph(i,j))^2*(xe(1)-xb(1))/2*w(j)); %calculate normalized energy norm numerator with GQ%
        eene(i)=eene(i)+(kk(i,j)*(upe(i,j))^2*(xe(1)-xb(1))/2*w(j)); %calculate normalized energy norm denominator with GQ%
        l2n(i)=l2n(i)+((ue(i,j)-uh(i,j))^2*(xe(1)-xb(1))/2*w(j)); %calculate normalized L2 norm numerator with GQ%
        l2ne(i)=l2ne(i)+((ue(i,j))^2*(xe(1)-xb(1))/2*w(j)); %calculate normalized L2 norm denominator with GQ%
    end
    
end
ee=sqrt(sum(ee)); %sum elemental error to get total error%
een=sqrt(sum(een))/sqrt(sum(eene)); %sum elemental error to get total error%
l2n=sqrt(sum(l2n))/sqrt(sum(l2ne)); %sum elemental error to get total error%