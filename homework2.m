clear
clear workspace;
clc

v=zeros(1,N);
u=zeros(1,N);
w=zeros(1,N);

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

        u(i,j)=a*com;
        v(i,j)=b*com;
        w(i,j)=c*com;

    end
end

% bc vortex

for i=1:M
    for j=i:M
    % ab vortex
        a=0;
        b=0;
        c=(xcp(i)-xb(j))*(ycp(i)-yc(j))-(ycp(i)-yb(j))*(xcp(i)-xc(j));
        d=a*a+b*b+c*c;

        r1= sqrt((xcp(i)-xb(j))*(xcp(i)-xb(j))+(ycp(i)-yb(j))*(ycp(i)-yb(j)));
        r2= sqrt((xcp(i)-x(j))*(xcp(i)-xc(j))+(ycp(i)-yc(j))*(ycp(i)-yb(j)));
        
        ror1=(xb(j)-xa(j))*(xcp(i)-xa(j))+(yb(i)-ya(j))*(ycp(i)-ya(j));
        ror2=(xb(j)-xa(j))*(xcp(i)-xb(j))+(yb(i)-ya(j))*(ycp(i)-yb(j));

        com=(1/(4*pi*d))*((ror1/r1)-(ror2/r2));

        u(i,j)=a*com;
        v(i,j)=b*com;
        w(i,j)=c*com;

    end
end