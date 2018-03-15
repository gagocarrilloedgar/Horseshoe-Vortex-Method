%% Inputs del sistema: Taper ratio(lambda), aspect ratio A, 
% quarter chord sweep angle, b=1, c=1, wing tip twistn, alphal0 (root and
% tip, flap position ( y (incio-fin)), alpha_g. 

disp('WING PARAMETERS');
% Aspect Ratio
inputwing =0;
while (inputwing==0)
    
    aspect_ratio=input('Introduce a valid aspect ratio A: ');
    if ( aspect_ratio<0)
        disp ('Invalid Aspect Ratio, introduce another one ');
        fprintf('\n');
        inputwing =0;
    elseif (aspect_ratio<5)
        warning('Horshoe Vortex Method may not provide accurate results for such a low aspect ratio');
        inputwing=1;
    else 
     S=1/aspect_ratio;
     inputwing=1;   
    end
    
end

%Taper Ratio
inputwing =0;
while (inputwing==0)
    
   taper_ratio=input('Taper ratio \lamda: ');
    if ( taper_ratio==0 || taper_ratio>1)
        disp ('Invalid Taper Ratio, set another one ');
        fprintf('\n');
        inputwing =0;
    else 
     cr=(2/aspect_ratio)/(1+taper_ratio);
     ct=cr*taper_ratio;
     inputwing=1;   
    end
    
end



disp('GEOMETRIC TWIST');
epsilon_ct = input('Tip airfoil epsilon (º): ');
epsilon_ct= deg2rad(epsilon_ct);


disp(' AERODYNAMIC CONDITIONS ');

alphal0t = input('Zerolift tip alpha (º): ');
alphal0t=deg2rad(alphal0t);
alphal0r = input('Zerolift root alpha(º): ');
alphal0r= deg2rad(alphal0r);

N=100;
CLA=zeros(61,1);
CL1=zeros(61,1);
CL2=zeros(61,1);
Cly1=zeros(N,1);
Cly2=zeros(N,1);
CL_alpha=zeros(61,1);
Cla=zeros(N,1);
Clb=zeros(N,1);
contador=1;
XacT=zeros(61,1);

for t=0:60
%Quarter chord swep angle (delta25)
tandeltaLE=tan(deg2rad(t));
delta25=atan(tandeltaLE-2*0.25*(cr-ct));

for h=5:5:10
alpha = deg2rad(h); 
beta=deg2rad(0);
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
chord1=zeros(1,N);
chord2=zeros(1,N);
chord3=zeros(1,N);

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
   chord3(i)=cr-((cr-ct)*abs(yvhm(i))*2);
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
CL=0; CMle=0; CDi=0;
sumaCL=0; sumaCMle=0; sumaCDi=0;
cmg=(2/3)*cr*((1+taper_ratio+(taper_ratio)^2)/(1+taper_ratio));

for i=1:N
    CL=(2/S)*Circulacion(i)*(1/N)+sumaCL;
    sumaCL=CL;
    CMle=-(2/(S*cmg))*Circulacion(i)*xvhm(i)*(1/N)+sumaCMle;
    sumaCMle=CMle;
    Cly(i)=(2/s(i))*Circulacion(i)*(1/N);
for j=1:N
    CDi=-(2/S)*Circulacion(i)*(1/N)*Circulacion(j)*Wm(i,j)+sumaCDi; 
    sumaCDi=CDi;
end
end

if h==5
CL1(contador)=CL;
for i=1:N
Cly1(i)=Cly(i);
end

else
CL2(contador)=CL;
for i=1:N
Cly2(i)=Cly(i);
end
end
end


CL_alpha(contador)= (CL1(contador))/(deg2rad(5));


for i=1:N
Cla(i)=(Cly1(i)-Cly2(i))/(CL1(contador)-CL2(contador));
end

for i=1:N
Clb(i)=Cly1(i)-(Cla(i)*CL1(contador));
end

sumaCla=0;
for i=1:N
suma=(1/S)*Cla(i)*chord3(i)*(1/N)+sumaCla;
sumaCla=suma;

end
CLA(contador)=suma;

sumaXac=0;
for i=1:N
Xac=((1/(S*cmg))*(Cla(i))*chord3(i)*abs(yvhm(i))*tan(delta25)*(1/N))+sumaXac;
sumaXac=Xac;
end

XacT(contador,1)=Xac+(cr*0.25)/cmg;
contador=contador+1;
end

X=0:1:60;
X2=[0 30 45 60];
Y=[4.2 3.9 3.4 2.6];
figure
plot(X,CL_alpha,'DisplayName','CL_{\alpha}');
hold on
plot(X,XacT,'r','DisplayName','x_{ac/c}');
plot(X2,Y,'o','DisplayName','Katz and Plotkin CL_{\alpha}');
xlabel('delta','interpreter', 'latex' );
ylabel('$CL_{\alpha}$','interpreter','latex');
hold off;
legend show;
legend('Location','southeast');
figure;
plot(yvhm, Cla,'o');
