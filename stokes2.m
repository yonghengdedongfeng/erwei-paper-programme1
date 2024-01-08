%这其实是对下半部分辐射的一个 特例推导，另一个是通用型推导，特例推导主要是验证性的，包括一个画偏振椭圆信息的模块需要被用到
clear;
clc;
nm = 1e-9;
lamda = 633*nm;
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
%p = [-0.9797-0.2003j,0,0.0352-0.2234j];

for ampz = (0:0.5:0.5)*2*2
    %这是颗粒散射偶极矩
    %这个改变的是pz相位，很有意思
    %p = [1,1j,1*exp(1j*ampz*pi/2)]; %ampz = (-0.5:0.1:0.5)*2
    %p = [1,1j,0.2*exp(1j*ampz)]; %ampz = (-0.5:0.1:0.5)*2*pi
    p = [2.58,-2.41*1j,1*1j];
    p = [2.58,-2*1j,1*1j];
     %下半区域
    nEfp = (-1*p(1)*C.*KX.*KZ)./(sqrt(KX.^2+KY.^2).*K0) + (-1*p(2)*C.*KY.*KZ)./(sqrt(KX.^2+KY.^2).*K0) - p(3)*C.*sqrt(KX.^2+KY.^2)./K0;
    nEfs = (-1*p(1)*C.*KY./sqrt(KX.^2+KY.^2)) + p(2)*C.*KX./sqrt(KX.^2+KY.^2);
    nI = abs(nEfs).^2 + abs(nEfp).^2;

    sinphi = KY./sqrt(KX.^2+KY.^2);
    cosphi = KX./sqrt(KX.^2+KY.^2);
    
    %下半区域%%%%%%%%%%%%%%%%%%%%
%     nEx = -1*nEfs.*sinphi - nEfp.*cosphi;
%     nEy =    nEfs.*cosphi - nEfp.*sinphi;
%     nElp = (nEx - 1j*nEy)/sqrt(2);
%     nErp = (nEx + 1j*nEy)/sqrt(2);
    %这个并不是投影计算,这个在补充位相
    pl = (nEfp - 1j*nEfs).*exp(-1j*atan2(KX,KY))/sqrt(2);
    pr = (nEfp + 1j*nEfs).*exp(1j*atan2(KX,KY))/sqrt(2);
    nIrp =abs(pr).^2;
    nIlp =abs(pl).^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%下半空间的
    nf = figure(2);
    nf.Position(1:2) = [1100 200];
    nf.Position(3:4) = [920 800];
    
    sgtitle(num2str(0.3))
    
    subplot(3,2,1)
    surf(KX/k0,KY/k0,-KZ/k0,nI)
    axis equal
    colorbar
    shading interp;
    xlabel('kx/k')
    ylabel('ky/k')
    zlabel('kz/k')

    subplot(3,2,2)
    imagesc(kx/k0,ky/k0,nI);title('z<0')
    colormap("jet")
    colorbar
    xlabel('kx/k')
    ylabel('ky/k')
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,3)
    %imagesc(kx/k0,ky/k0,log(abs(nErp).^2));title('Irp')
    imagesc(kx/k0,ky/k0,abs(pr).^2);title('Irp')
    colorbar
    xlabel('kx/k')
    ylabel('ky/k')
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off
    
    subplot(3,2,4)
    imagesc(kx/k0,ky/k0,abs(pl).^2);title('Ilp')
    colorbar
    xlabel('kx/k')
    ylabel('ky/k')
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,5)
    imagesc(kx/k0,ky/k0,angle(pl));title('phase(left)')
%     imagesc(kx/k0,ky/k0,angle(nErp));title('phase(right)')
    colorbar
    xlabel('kx/k')
    ylabel('ky/k')
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    subplot(3,2,6)
%     imagesc(kx/k0,ky/k0,angle(nElp));title('phase(left)')
    imagesc(kx/k0,ky/k0,angle(pr));title('phase(right)')
    colorbar
    xlabel('kx/k')
    ylabel('ky/k')
    axis xy
    hold on
    k01 = ones(1,100)*k0;
    x=linspace(-k0,k0,100);
    plot(x/k0,sqrt(k01.^2-x.^2)/k0,'g--',x/k0,-sqrt(k01.^2-x.^2)/k0,'g--')
    hold off

    pause(3)
end

%椭圆反演模块
%%
clear;
nm = 1e-9;
lamda = 633*nm;
k0 = 2*pi/lamda;

k1x = k0*(-0.499);
k1y = k0*(-0.388);
k2x = k0*(0.499);
k2y = k0*(-0.388);
p1 = [k1x;k1y;-sqrt(k0^2-k1x^2-k1y^2)];
p2 = [k2x;k2y;-sqrt(k0^2-k2x^2-k2y^2)];
figure(1)
quiver3(0,0,0,p1(1),p1(2),p1(3),0,'Color','r','LineWidth',3);
hold on
quiver3(0,0,0,p2(1),p2(2),p2(3),0,'Color','b','LineWidth',5);
hold off
axis xy
xlabel('x')
ylabel('y')

phi1 = atan2(p1(2),p1(1));
costhe = sqrt(1-(p1(1)^2+p1(2)^2)/k0^2);
Rzp = [cos(phi1),sin(phi1),0;-sin(phi1),cos(phi1),0;0,0,1];
Ryt = [costhe,0,sqrt(1-costhe^2);0,1,0;-sqrt(1-costhe^2),0,costhe];

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
the2 = acos(abs(p22(3))/k0);
phi2 = atan2(p22(2),p22(1));
%这里问题很大
Ex2 = cos(phi2)*(1j) - sin(phi2)*(1);
Ey2 = sin(phi2)*(1j) + cos(phi2)*(1);
Ez2 = tan(the2/2)*(1j);

E2 = [Ex2;Ey2;Ez2];
E = transpose(Rzp)*transpose(Ryt)*E2

%画偏振椭圆模块
%%
%对上面求出来的东西进行检验
h=1;

for phi=0.1:0.1:2*pi
    Ex1 = 2.58;
    Ey1 = -2j;
    Ez1 = 1j;

    t=linspace(0,2*pi,50);
    x = real(Ex1*exp(1j*t));
    y = real(Ey1*exp(1j*t));
    z = real(Ez1*exp(1j*t));

    mochang = abs(sqrt(Ex1^2+Ez1^2+Ey1^2));
    changshu = sqrt(Ex1^2+Ey1^2+Ez1^2);
    majorx = real(conj(Ex1)*changshu)/mochang;
    majory = real(conj(Ey1)*changshu)/mochang;
    majorz = real(conj(Ez1)*changshu)/mochang;
    sqrt(majorx^2+majory^2+majorz^2)
    minx = imag(conj(Ex1)*changshu)/mochang;
    miny = imag(conj(Ey1)*changshu)/mochang;
    minz = imag(conj(Ez1)*changshu)/mochang;
    sqrt(minx^2+miny^2+minz^2)
    E1 = [Ex1,Ey1,Ez1];
    normal = imag(cross(conj(E),E));
    n = sqrt(normal(1)^2+normal(2)^2+abs(normal(3))^2);
    normal = normal/n;
    
end

for phi=0.1:0.1:2*pi
    Ex2 = 1.1539j;
    Ey2 = 0.8942;
    Ez2 =-0.4477;

    t=linspace(0,2*pi,50);
    x2 = real(Ex2*exp(1j*t));
    y2 = real(Ey2*exp(1j*t));
    z2 = real(Ez2*exp(1j*t));

    mochang = abs(sqrt(Ex2^2+Ez2^2+Ey2^2));
    changshu = sqrt(Ex2^2+Ey2^2+Ez2^2);
    majorx = real(conj(Ex2)*changshu)/mochang;
    majory = real(conj(Ey2)*changshu)/mochang;
    majorz = real(conj(Ez2)*changshu)/mochang;
    sqrt(majorx^2+majory^2+majorz^2)
    minx = imag(conj(Ex2)*changshu)/mochang;
    miny = imag(conj(Ey2)*changshu)/mochang;
    minz = imag(conj(Ez2)*changshu)/mochang;
    sqrt(minx^2+miny^2+minz^2)
    E2 = [Ex2,Ey2,Ez2];
    normal = imag(cross(conj(E),E));
    n = sqrt(normal(1)^2+normal(2)^2+abs(normal(3))^2);
    normal = normal/n;
    
end

plot3(x,y,z)
hold on
plot3(x2,y2,z2)
hold on 
plot(x,y)
hold on
q1=quiver3(0,0,0,majorx,majory,majorz,0,'Color','r','LineWidth',3);
q2=quiver3(0,0,0,minx,miny,minz,0,'Color','b','LineWidth',5);
q3=quiver3(0,0,0,normal(1),normal(2),normal(3),0,'Color','g','LineWidth',5);
xlabel('x')
ylabel('y')
zlabel('z')
axis xy
