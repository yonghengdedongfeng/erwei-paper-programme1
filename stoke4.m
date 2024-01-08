%%
%偏振奇点计算部分
clear;
clc;

nm = 1e-9;
lamda = 627*nm;
%高度要注意调节，越高则倏逝波越小
d = 40*nm;
k0 = 2*pi/lamda;
n = 1;
NA = n*sin(pi/2);

kx = linspace(-n*k0,n*k0,2001);
ky = linspace(-n*k0,n*k0,2001);
[KX,KY] = meshgrid(kx,ky);
for i=1:2001
    for j=1:2001
        if KX(i,j)^2+KY(i,j)^2 > 0.9991*k0^2
            KX(i,j) = 0;
            KY(i,j) = 0;
        end
    end
end
K0 = ones(2001)*k0;
KZ = sqrt(K0.^2-KX.^2-KY.^2);
KZ2 = sqrt(n^2*K0.^2-KX.^2-KY.^2);
C = exp(1j*KZ*d).*sqrt(n^2*K0.^2-KX.^2-KY.^2)./KZ;
TP = 2*n*KZ./(KZ2 + n^2*KZ);
TS = 2*KZ./(KZ + KZ2);

for ampz = (0:0.025:0.5)*2*2
    %这是颗粒散射偶极矩
    p = [1,1j,ampz*1j]; %ampz = (-0.5:0.1:0.5)*2
    p = [2.58,-2*1j,1*1j];
    %p = [1,1j,0.2*exp(1j*ampz)]; %ampz = (-0.5:0.1:0.5)*2*pi
    %p = [0.76,0,1*1j];
    
    %上半区域
    pEfp = (p(1)*C.*KX.*KZ)./(sqrt(KX.^2+KY.^2).*K0) + (p(2)*C.*KY.*KZ)./(sqrt(KX.^2+KY.^2).*K0) - p(3)*C.*sqrt(KX.^2+KY.^2)./K0;
    pEfs = (-1*p(1)*C.*KY./sqrt(KX.^2+KY.^2)) + p(2)*C.*KX./sqrt(KX.^2+KY.^2);

    pI = abs(pEfs).^2 + abs(pEfp).^2;

    sinphi = KY./sqrt(KX.^2+KY.^2);
    cosphi = KX./sqrt(KX.^2+KY.^2);
%     %上半区域
%     pEx = -1*pEfs.*sinphi + pEfp.*cosphi;
%     pEy =    pEfs.*cosphi + pEfp.*sinphi;
%     %这里是圆偏态的共轭转置处理之后
%     pElp = (pEx + 1j*pEy)/sqrt(2);
%     pErp = (pEx - 1j*pEy)/sqrt(2);

    pl = (pEfp + 1j*pEfs).*exp(-1j*atan2(KX,KY))/sqrt(2);
    pr = (pEfp - 1j*pEfs).*exp(1j*atan2(KX,KY))/sqrt(2);
    pIrp =abs(pr).^2;
    pIlp =abs(pl).^2;
     %画图
    pf = figure(1);
    pf.Position(1:2) = [1100 200];
    pf.Position(3:4) = [920 800];
   
    sgtitle(num2str(ampz))

    subplot(3,2,1)
    surf(KX/k0,KY/k0,KZ/k0,pI)
    axis equal
    colorbar
    shading interp;
    xlabel('kx/k')
    ylabel('ky/k')
    zlabel('kz/k')

    subplot(3,2,2)
    imagesc(kx/k0,ky/k0,pI);title('z>0')
    colormap("jet")
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,3)
    imagesc(kx/k0,ky/k0,abs(pr).^2);title('Irp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,4)
    imagesc(kx/k0,ky/k0,abs(pl).^2);title('Ilp')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,5)
%     imagesc(kx/k0,ky/k0,angle(pElp));title('phase(left)')
    imagesc(kx/k0,ky/k0,angle(pr));title('phase(right)')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,6)
    imagesc(kx/k0,ky/k0,angle(pl));title('phase(left)')
%     imagesc(kx/k0,ky/k0,angle(pErp));title('phase(right)')
    colorbar
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    pause(3)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%偏振椭圆反演部分
%%
clear;
nm = 1e-9;
lamda = 633*nm;
k0 = 2*pi/lamda;

% k1x = k0*(0.499);
% k1y = k0*(0.388);
% k2x = k0*(-0.499);
% k2y = k0*(0.388);
k1x = k0*(-0.3386);
k1y = k0*(0.3298);
k2x = k0*(-0.1102);
k2y = k0*(0.0867);
p1 = [k1x;k1y;sqrt(k0^2-k1x^2-k1y^2)];
p2 = [k2x;k2y;sqrt(k0^2-k2x^2-k2y^2)];
figure(1)
quiver3(0,0,0,p1(1),p1(2),p1(3),0,'Color','r','LineWidth',3);
hold on
quiver3(0,0,0,p2(1),p2(2),p2(3),0,'Color','b','LineWidth',5);
hold off
axis xy
xlabel('x')
ylabel('y')

phi1 = atan2(p1(2),p1(1));
%costhe = sqrt(1-(p1(1)^2+p1(2)^2)/k0^2);
costhe = p1(3)/k0;
Rzp = [cos(phi1),sin(phi1),0;-sin(phi1),cos(phi1),0;0,0,1];
Ryt = [costhe,0,-sqrt(1-costhe^2);0,1,0;sqrt(1-costhe^2),0,costhe];

p11 = mtimes(Rzp,p1);
p12 = mtimes(Ryt,p11);
p21 = mtimes(Rzp,p2);
p22 = mtimes(Ryt,p21);
figure(2)
quiver3(0,0,0,p12(1),p12(2),p12(3),0,'Color','r','LineWidth',3);
hold on
quiver3(0,0,0,p22(1),p22(2),p22(3),0,'Color','b','LineWidth',5);
axis xy
xlabel('x')
ylabel('y')

%防止计算出现的近似零的小数两两出现，导致不能够按照作为零来处理
for i = 1:1:3
    commax = max(p22);
    if abs(p22(i)/commax)<=1e-10
        p22(i) = 0;
    end
end
the2 = acos(p22(3)/k0);
phi2 = atan2(p22(2),p22(1));
%这里问题很大
% Ex2 = cos(phi2+pi)*(1j) + cos(phi2+3*pi/2)*(1);
% Ey2 = cos(phi2+pi/2)*(1j) + cos(phi2+pi)*(1);
% Ez2 = tan(the2/2)*(1j);
Ex2 = cos(phi2+pi)*(-1j) + cos(phi2+3*pi/2)*(1);
Ey2 = cos(phi2+pi/2)*(-1j) + cos(phi2+pi)*(1);
Ez2 = tan(the2/2)*(-1j);

E2 = [Ex2;Ey2;Ez2];
E = transpose(Rzp)*transpose(Ryt)*E2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%新归一化归一化过程
% a1 = 2382-3086j;
% a2 = 3166+2275j;
% a3 = -168-1151j;
a1 = -0.7241+0.6612j;
a2 = -0.6894-0.7037j;
a3 = -0.0202+0.3152j;
pp1 = 0.5*angle(a1*a1+a2*a2+a3*a3);
AA = sqrt(abs(a1)^2 + abs(a2)^2 +abs(a3)^2);
ee = [a1,a2,a3]*exp(-1j*pp1)/AA

