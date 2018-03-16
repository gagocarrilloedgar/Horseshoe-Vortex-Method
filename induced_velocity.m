clear
clear workspace;
clc
%AB
v1=zeros(M,M);
u1=zeros(M,M);
w1=zeros(M,M);
%BC
v2=zeros(M,M);
u2=zeros(M,M);
w2=zeros(M,M);
%CD
v3=zeros(M,M);
u3=zeros(M,M);
w3=zeros(M,M);

pi=3.141592;

%% Induced Velocity vector 
for i=1:M
    for j=i:M
    % ab vortex
        a=0;
        b=0;
        c=(xcp(i)-xa(j))*(ycp(i)-yb(j))-(ycp(i)-ya(j))*(xcp(i)-xb(j));
        d=a*a+b*b+c*c;

        r1= sqrt((xcp(i)-xa(j))*(xcp(i)-xa(j))+(ycp(i)-ya(j))*(ycp(i)-ya(j)));
        r2= sqrt((xcp-xb)*(xcp-xb)+(ycp-yb)*(ycp-yb));

        ror1=(xb(j)-xa(j))*(xcp(i)-xa(j))+(yb(i)-ya(j))*(ycp(i)-ya(j));
        ror2=(xb(j)-xa(j))*(xcp(i)-xb(j))+(yb(i)-ya(j))*(ycp(i)-yb(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u1(i,j)=a*com;
        v1(i,j)=b*com;
        w1(i,j)=c*com;

    end
end

% bc vortex
for i=1:M
    for j=i:M
        a=0;
        b=0;
        c=(xcp(i)-xb(j))*(ycp(i)-yc(j))-(ycp(i)-yb(j))*(xcp(i)-xc(j));
        d=a*a+b*b+c*c;

        r1= sqrt((xcp(i)-xb(j))*(xcp(i)-xb(j))+(ycp(i)-yb(j))*(ycp(i)-yb(j)));
        r2= sqrt((xcp(i)-xc(j))*(xcp(i)-xc(j))+(ycp(i)-yc(j))*(ycp(i)-yc(j)));
        
        ror1=(xc(j)-xb(j))*(xcp(i)-xb(j))+(yc(i)-yb(j))*(ycp(i)-yb(j));
        ror2=(xc(j)-xb(j))*(xcp(i)-xc(j))+(yc(i)-yb(j))*(ycp(i)-yc(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u2(i,j)=a*com;
        v2(i,j)=b*com;
        w2(i,j)=c*com;

    end
end

%cd vortex
for i=1:M
    for j=i:M
        a=0;
        b=0;
        c=(xcp(i)-xc(j))*(ycp(i)-yd(j))-(ycp(i)-yc(j))*(xcp(i)-xd(j));
        d=a*a+b*b+c*c;

        r1= sqrt((xcp(i)-xc(j))*(xcp(i)-xc(j))+(ycp(i)-yc(j))*(ycp(i)-yc(j)));
        r2= sqrt((xcp(i)-xd(j))*(xcp(i)-xd(j))+(ycp(i)-yd(j))*(ycp(i)-yd(j)));
        
        ror1=(xd(j)-xc(j))*(xcp(i)-xc(j))+(yd(i)-yc(j))*(ycp(i)-yc(j));
        ror2=(xd(j)-xc(j))*(xcp(i)-xd(j))+(yd(i)-yc(j))*(ycp(i)-yd(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u3(i,j)=a*com;
        v3(i,j)=b*com;
        w3(i,j)=c*com;

    end
end

% Velocidades inducidadas totales en todo el panel
U=u1+u2+u3;
V=v1+v2+v3;
W=w1+w2+w3;