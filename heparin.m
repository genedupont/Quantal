%ROD kinetic scheme with addition of heparin
clear
clc
for iii=0:1
    for jjj=1:1
%Kinetic constants of the model
k1=0.4;
k2=0.2;
k3=0.004;
km1=k1/100;
km2=k2/100;

Rt=1.966*10^8; %Number of IP3 receptors in the well divided by 100 to keep short computing times
%Rt=1.966*10^7; %Number of IP3 receptors in the well divided by 1000 to keep short computing times
%Rt=1.966*10^6; %Number of IP3 receptors in the well divided by 10000 to keep short computing times


N=[0,0,30786,0,0]; %Number of tetraliganded receptors at 100nM IP3 divided by 100, as for Rtot
%N=[0,0,1.112*10^5,0,0]; %Number of tetraliganded receptors at 300 nM IP3 divided by 1000, as for Rtot
%N=[0,0,1.9*10^5,0,0]; %Number of tetraliganded receptors at 1 microM IP3 divided by 10000, as for Rtot


hep=iii;
NN=jjj;
tau=100;
krel=40; %40/Number of receptors per cell(5e^4) = 0.0008
b=0.0000001;

if hep == 1
    couleur='b';
else
    couleur='r';
end

C=0.22;
Cer(1,1)=500;

dt=1/0.4/100 %Time step (0.4 being the largest kinetic constant)
tf=300; %End time
dtsave=1/0.2/200 %Time step for saving data. 
tsave(:,1)=0:dtsave:tf; % Time vector for saving data
%Matrix for saving data 
Osave=zeros(1/dtsave*tf+1,1);
Rs=ones(size(Osave,1),1);
Os=zeros(size(Osave,1),1);
Ds=zeros(size(Osave,1),1);

for i=1:3 %Number of pulses
    for j=N(i)+1:N(i+1) %Select each IP3R bound to IP3 
        t=(i-1)*tf/4;
        %These parameters allows to consider the reaction in a
        %non-stochastic manner if needed
        to=0;   
        tr=0;
        td=0;
        nsave=(size(Osave,1)-1)/4*(i-1)+1;
        R=1;
        O=0;
        D=0;
        
        
        for step=1/dt*tf/4*(i-1)+1:1/dt*tf %For each time step        
            if R==1 
                PDR=k1*dt;
                Pnothing=1-PDR;
                if rand < PDR
                    R=0;
                    O=1;
                    to=t;
                end
            elseif O==1
                PDR1=k2*dt;
                PDR2=km1*dt;
                summa=PDR1+PDR2;
                Pnothing=1-PDR1-PDR2;
                rr=rand;
                if rr < PDR1
                    O=0;
                    D=1;
                    td=t;
                elseif (rr > PDR1) && (rr < summa)
                    O=0;
                    R=1;
                    td=t; 
                end
            elseif D==1 
                 PDR1=k3*dt;
                PDR2=km2*dt;
                summa=PDR1+PDR2;
                Pnothing=1-PDR1-PDR2;
                rr=rand;
                if rr < PDR1
                    D=0;
                    R=1;
                    td=t;
                elseif (rr > PDR1) && (rr < summa)
                    O=1;
                    D=0;
                    td=t;
                end
            end
            t=t+dt;
            
            if (t-tsave(nsave,1))>dtsave
                nsave=nsave+1;
                Rs(nsave,1)=Rs(nsave,1)+R;
                Os(nsave,1)=Os(nsave,1)+O;
                Ds(nsave,1)=Ds(nsave,1)+D;
            end
        end
    end
end



size(tsave)
size(Os)
% 
% figure(1)
% plot(tsave, sum(Os,2))
% hold on
% xlabel('Time (s)')
% ylabel('Number of open IP_{3}Rs')
% figure(2)
% plot(tsave, sum(Ds,2))
% hold on
% xlabel('Time (s)')
% ylabel('Number of desensitised IP_{3}Rs')

 

for i=2:size(tsave,1)
    dt=tsave(i,1)-tsave(i-1,1);
    o4=Os(i,1);
    i;
%     if Cer(i-1,1)>500
%         Vp=0;
%     else
%         Vp=VP0;
%     end
    if tsave(i,1)<100
        krel=krel;
    else
        krel=krel*hep;
    end
    Cer(i,1)=rg4dcer(dt,tau,krel, Cer(i-1,1),o4,Rt,b) ;


end
figure(jjj)
hold on
%plot(tsave,(100-15)*(Cer./5/100)+15,couleur)
plot(tsave,(100-20)*(Cer./5/100)+20,couleur)
Cerbis=(100-20)*(Cer./5/100)+20;
xlabel('Time (s)')
ylabel('Ca^{2+} content (%)')
ylim([20 100]);

%pour sauver les données
if (iii == 0)
    T1 = table(tsave,Cerbis,'VariableNames',{ 'time', 'calcium'});
    writetable(T1, 'calcium.txt');
else
    T1bis = table(tsave,Cerbis,'VariableNames',{ 'time', 'calcium'});
    writetable(T1bis, 'calcium_hep.txt');
end

    end
end

    function y=dcer(Cer,tau,krel,o4,Rtot,b) 
    %y=Vp*(Kcer)^6/((Cer)^6+Kcer^6)*0.22^2/(0.22^2+0.35^2)-krel*(Cer-0.22)*(b+o4/Rtot);
    y=(500-Cer)/tau-krel*(Cer-0.22)*(b+o4/Rtot);
    end
    function y=rg4dcer(dt,tau,krel,Cer,o4,Rtot,b) 
    k1=dcer(Cer,tau,krel,o4,Rtot,b) ;
    k2=dcer(Cer+dt/2*k1,tau,krel,o4,Rtot,b) ;
    k3=dcer(Cer+dt/2*k2,tau,krel,o4,Rtot,b) ;
    k4=dcer(Cer+dt*k3,tau,krel,o4,Rtot,b) ;
    y=Cer+dt/6*(k1+2*k2+2*k3+k4);
    end


