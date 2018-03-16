%% calculo de la polar 

CDP=zeros(1,N);
CDy=zeros(1,N);
CDiy=zeros(1,N);
CLy=zeros(1,N);

for i=1:N
    CDP(i)=(1/S)*(1/N)*(0.00585*(cly(i))^2 -0.01867*cly(i)+0.0085)*chord3(i); 
    
        for j=1:N
        CDiy(i)=-(2/S)*Circulacion(i)*(1/N)*Circulacion(j)*Wm(i,j);
        end
    
    CDy(i)= CD(i)+CDP(i);
    CLy=(2/S)*Circulacion(i)*(1/N);
end

figure;
plot(CLy,CD);
xlim([-0.5,0,5]);
title('Polar wing curve CD as function od CL','interpreter','latex');
xlabel('CL','interpreter','latex');
ylabel('CD''interpreter','latex');
grid on;
