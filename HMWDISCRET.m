clear
clear workspace;
clc



%% Wing plot 
N=8; %%Numero de HV panels
delta25=deg2rad(0);
A=4;
estrechamiento=1;
cr=(2/A)/(estrechamiento+1);
ct=estrechamiento*cr;
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
%% coordenas HV control points
xcp=zeros(1,N);
ycp=zeros(1,N);
epsilont=deg2rad(-3);
alphal0t=deg2rad(-2);
alphal0r=deg2rad(-1);
epsilon=zeros(1,N);
alphal0=zeros(1,N);
alphal0tot=zeros(1,N);
nx=zeros(1,N);
ny=zeros(1,N);
nz=zeros(1,N);

for i=1:(N/2)
   yvhm(i)=y(i)+((1/N)/2);
   xvhm(i)=(cr/4)+(abs(yvhm(i))*tan(delta25));
   xcp(i)=((3*cr)/4)-(abs(yvhm(i))*tandelta75);
   ycp(i)= yvhm(i);
   epsilon(i)=epsilont*2*abs(ycp(i));
   alphal0(i)=(2*abs(ycp(i))*(alphal0t-alphal0r))+alphal0r;
end


for i=((N/2)+2):(N+1)
   yvhm(i-1)=y(i)-((1/N)/2);
   xvhm(i-1)= (cr/4)+(yvhm(i-1)*tan(delta25));
   xcp(i-1)=((3*cr)/4)-(yvhm(i-1)*tandelta75);
   ycp(i-1)=yvhm(i-1);
   epsilon(i-1)=epsilont*2*abs(ycp(i-1));
   alphal0(i-1)=(2*abs(ycp(i-1))*(alphal0t-alphal0r))+alphal0r;
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
end

 %%Calculo del vector n para cada panel (sin tener en cuenta flap)
 for i=1:N
 alphal0tot(i)=epsilon(i)+alphal0(i);
 nx(i)=sin(alphal0tot(i));
 ny(i)=0;
 nz(i)=cos(alphal0tot(i));
 end
 
 


       
figure
plot(xvhm,yvhm,'c*');
hold on;
plot(xa,ya,'c*');
plot(xb,yb,'c*');
plot(xc,yc,'c*');
plot(xd,yd,'c*');
plot(xcp,ycp,'c*');
hold off;

