clear
clear workspace;
clc

%% Inputs del sistema: Taper ratio(lambda), aspect ratio A, 
% quarter chord sweep angle, b=1, c=1, wing tip twistn, alphal0 (root and
% tip, flap position ( y (incio-fin)), alpha_g. 


% We suppose symmetrical profile with no twist nor sweep.
aspect_ratio=10;
delta25=0;
epsilon_ct=0;
alphal0t=0;
alphal0r=0;
taper_ratio=0.2;
S=1/(aspect_ratio);
N=100;
Cly_graf=zeros(5,N);
z=1;
%To impose CL=0.25, we get two values of CL at random angles of attack
%(between  and 10º) and we calculate the slope, so we can know the angle of
%attack to have CL=0.25. We will use, for example 2º and 5º
while(taper_ratio<=1);
cr=(2/aspect_ratio)/(1+taper_ratio);
ct=cr*taper_ratio;
CLx=zeros(1,2);
x=1;
alpha=2;
alpha = deg2rad(alpha);  
beta = deg2rad(0);
for i=1:2
Vinf = [cos(alpha)*cos(beta) cos(alpha)*sin(beta) sin(alpha)];

%% Wing plot 

%%Nos hace falta la tangenete del angulo que forma la linea de 3c/4 con el eje y 
%%para el computo de los control points
tandelta75=2*((0.75*cr)-(0.25*cr+0.5*tan(delta25)+0.5*ct));
y=-0.5:(1/N):0.5;
xvhm=zeros(1,N); %% X Vortex head midpoints
yvhm=zeros(1,N); %% Y Vortex head midpoints

%% coordenas HV panels
xa=zeros(1,N); 
xb=zeros(1,N);
xc=zeros(1,N);
xd=zeros(1,N);
ya=zeros(1,N);
yb=zeros(1,N);
yc=zeros(1,N);
yd=zeros(1,N);
s=zeros(1,N);
chord=zeros(1,N);
%% coordenas HV control points

xcp=zeros(1,N);
ycp=zeros(1,N);
zcp=zeros(1,N);
epsilon=zeros(1,N);
alphal0=zeros(1,N);
alphal0tot=zeros(1,N);
nx=zeros(1,N);
ny=zeros(1,N);
nz=zeros(1,N);

for i=1:N
   yvhm(i)=y(i)+((1/N)/2);
   xvhm(i)=(cr/4)+(abs(yvhm(i))*tan(delta25));
   xcp(i)=((3*cr)/4)-(abs(yvhm(i))*tandelta75);
   ycp(i)= yvhm(i);
   epsilon(i)=epsilon_ct*2*abs(ycp(i));
   alphal0(i)=(2*abs(ycp(i))*(alphal0t-alphal0r))+alphal0r;
end

for i=1:N
   yb(i)=yvhm(i)-((1/N)*0.5);
   xb(i)=xvhm(i);
   yc(i)=yvhm(i)+((1/N)*0.5);
   xc(i)=xvhm(i);
   ya(i)=yb(i);
   xa(i)=20;
   yd(i)=yc(i);
   xd(i)=20;
   chord1(i)=cr-((cr-ct)*abs(yb(i))*2);
   chord2(i)=cr-((cr-ct)*abs(yc(i))*2);
   s(i)=(chord1(i)+chord2(i))*(1/N)*0.5;
end

 %%Calculo del vector n para cada panel (sin tener en cuenta flap)
 for i=1:N
 alphal0tot(i)=epsilon(i)+alphal0(i);
 nx(i)=sin(alphal0tot(i));
 ny(i)=0;
 nz(i)=cos(alphal0tot(i));
 end

 %% Calculo de las velocidades inducidas   
 
%AB
v1=zeros(N,N);
u1=zeros(N,N);
w1=zeros(N,N);
v1m=zeros(N,N);
u1m=zeros(N,N);
w1m=zeros(N,N);
%BC
v2=zeros(N,N);
u2=zeros(N,N);
w2=zeros(N,N);
%CD
v3=zeros(N,N);
u3=zeros(N,N);
w3=zeros(N,N);
v3m=zeros(N,N);
u3m=zeros(N,N);
w3m=zeros(N,N);


%% Induced Velocity (Control Point)
for i=1:N
    for j=1:N
    % ab vortex
        a=0;
        b=0;
        c=((xcp(i)-xa(j))*(ycp(i)-yb(j)))-((ycp(i)-ya(j))*(xcp(i)-xb(j)));
        d=a*a+b*b+c*c;

        r1= sqrt(((xcp(i)-xa(j))*(xcp(i)-xa(j)))+((ycp(i)-ya(j))*(ycp(i)-ya(j))));
        r2= sqrt(((xcp(i)-xb(j))*(xcp(i)-xb(j)))+((ycp(i)-yb(j))*(ycp(i)-yb(j))));

        ror1=(xb(j)-xa(j))*(xcp(i)-xa(j))+(yb(j)-ya(j))*(ycp(i)-ya(j));
        ror2=(xb(j)-xa(j))*(xcp(i)-xb(j))+(yb(j)-ya(j))*(ycp(i)-yb(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u1(i,j)=a*com;
        v1(i,j)=b*com;
        w1(i,j)=c*com;

    end
end

% bc vortex
for i=1:N
    for j=1:N
        a=0;
        b=0;
        c=((xcp(i)-xb(j))*(ycp(i)-yc(j)))-((ycp(i)-yb(j))*(xcp(i)-xc(j)));
        d=a*a+b*b+c*c;

        r1= sqrt(((xcp(i)-xb(j))*(xcp(i)-xb(j)))+((ycp(i)-yb(j))*(ycp(i)-yb(j))));
        r2= sqrt(((xcp(i)-xc(j))*(xcp(i)-xc(j)))+((ycp(i)-yc(j))*(ycp(i)-yc(j))));
        
        ror1=((xc(j)-xb(j))*(xcp(i)-xb(j)))+((yc(j)-yb(j))*(ycp(i)-yb(j)));
        ror2=((xc(j)-xb(j))*(xcp(i)-xc(j)))+((yc(j)-yb(j))*(ycp(i)-yc(j)));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u2(i,j)=a*com;
        v2(i,j)=b*com;
        w2(i,j)=c*com;

    end
end

%cd vortex
for i=1:N
    for j=1:N
        a=0;
        b=0;
        c=((xcp(i)-xc(j))*(ycp(i)-yd(j)))-((ycp(i)-yc(j))*(xcp(i)-xd(j)));
        d=a*a+b*b+c*c;

        r1= sqrt(((xcp(i)-xc(j))*(xcp(i)-xc(j)))+((ycp(i)-yc(j))*(ycp(i)-yc(j))));
        r2= sqrt(((xcp(i)-xd(j))*(xcp(i)-xd(j)))+((ycp(i)-yd(j))*(ycp(i)-yd(j))));
        
        ror1=((xd(j)-xc(j))*(xcp(i)-xc(j)))+((yd(j)-yc(j))*(ycp(i)-yc(j)));
        ror2=((xd(j)-xc(j))*(xcp(i)-xd(j)))+((yd(j)-yc(j))*(ycp(i)-yd(j)));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u3(i,j)=a*com;
        v3(i,j)=b*com;
        w3(i,j)=c*com;

    end
end

% Velocidades inducidas totales en el control point
U=u1+u2+u3;
V=v1+v2+v3;
W=w1+w2+w3;

%% Induced Velocity (Vortex Head Midpoint)
for i=1:N
    for j=1:N
    % ab vortex
        a=0;
        b=0;
        c=((xvhm(i)-xa(j))*(yvhm(i)-yb(j)))-((yvhm(i)-ya(j))*(xvhm(i)-xb(j)));
        d=a*a+b*b+c*c;

        r1= sqrt(((xvhm(i)-xa(j))*(xvhm(i)-xa(j)))+((yvhm(i)-ya(j))*(yvhm(i)-ya(j))));
        r2= sqrt(((xvhm(i)-xb(j))*(xvhm(i)-xb(j)))+((yvhm(i)-yb(j))*(yvhm(i)-yb(j))));

        ror1=(xb(j)-xa(j))*(xvhm(i)-xa(j))+(yb(j)-ya(j))*(yvhm(i)-ya(j));
        ror2=(xb(j)-xa(j))*(xvhm(i)-xb(j))+(yb(j)-ya(j))*(yvhm(i)-yb(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u1m(i,j)=a*com;
        v1m(i,j)=b*com;
        w1m(i,j)=c*com;

    end
end


%cd vortex
for i=1:N
    for j=1:N
        a=0;
        b=0;
        c=((xvhm(i)-xc(j))*(yvhm(i)-yd(j)))-((yvhm(i)-yc(j))*(xvhm(i)-xd(j)));
        d=a*a+b*b+c*c;

        r1= sqrt(((xvhm(i)-xc(j))*(xvhm(i)-xc(j)))+((yvhm(i)-yc(j))*(yvhm(i)-yc(j))));
        r2= sqrt(((xvhm(i)-xd(j))*(xvhm(i)-xd(j)))+((yvhm(i)-yd(j))*(yvhm(i)-yd(j))));
        
        ror1=((xd(j)-xc(j))*(xvhm(i)-xc(j)))+((yd(j)-yc(j))*(yvhm(i)-yc(j)));
        ror2=((xd(j)-xc(j))*(xvhm(i)-xd(j)))+((yd(j)-yc(j))*(yvhm(i)-yd(j)));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u3m(i,j)=a*com;
        v3m(i,j)=b*com;
        w3m(i,j)=c*com;

    end
end

% Velocidades inducidas totales en el control point
Um=u1m+u3m;
Vm=v1m+v3m;
Wm=w1m+w3m;


%% Matriz de coeficientes
B=zeros(N,1);
A=zeros(N,N);

for i=1:N
    for j=1:N
    A(i,j)=U(i,j)*nx(i)+V(i,j)*ny(i)+W(i,j)*nz(i);
    end
    B(i)=-Vinf(1,1)*nx(i)-Vinf(1,2)*ny(i)-Vinf(1,3)*nz(i);
    
end

%% Fuerzas
Cly=zeros(N,1);
Circulacion=linsolve(A,B);
CL=0; 
sumaCL=0;
for i=1:N
    CL=(2/S)*Circulacion(i)*(1/N)+sumaCL;
    sumaCL=CL;
end
CLx(x) = CL;
alpha = 5;
alpha = deg2rad(alpha);  
x=x+1;
end
m=(CLx(2)-CLx(1))/3;
%imponemos que el CL=0.25, y con el angulo de ataque que obtendremos
%aseguraremos un CL=0.25.
alpha_def=(0.25/m);
alpha_def=deg2rad(alpha_def);
Vinf = [cos(alpha_def)*cos(beta) cos(alpha_def)*sin(beta) sin(alpha_def)];
B=zeros(N,1);
A=zeros(N,N);
for i=1:N
    for j=1:N
    A(i,j)=U(i,j)*nx(i)+V(i,j)*ny(i)+W(i,j)*nz(i);
    end
    B(i)=-Vinf(1,1)*nx(i)-Vinf(1,2)*ny(i)-Vinf(1,3)*nz(i);
    
end

%% Fuerzas
Cly=zeros(N,1);
Circulacion=linsolve(A,B);
CL=0; 
sumaCL=0;
for i=1:N
    CL=(2/S)*Circulacion(i)*(1/N)+sumaCL;
    sumaCL=CL;
    Cly_graf(z,i)=(2/s(i))*Circulacion(i)*(1/N);
end
taper_ratio=taper_ratio+0.2;
z=z+1;
end
Ni=zeros(1,100);
for i=1:100;
    Ni(i)=i;
end
figure
plot(yvhm, Cly_graf);
title('Cl distribution along the panels of the wing');
xlabel('y');
ylabel('Cl');

grid on;
legend show;
legend ('{\lambda}=0.2','{\lambda}=0.4','{\lambda}=0.6','{\lambda}=0.8','{\lambda}=1','Location','southeast');
hold on;

