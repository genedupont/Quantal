%Minimalist model for quantal release (ROD). R->O and D->R are stochastic.
clear
clc
for iii=0:1
    for jjj=1:1
%Parameters of the model
K=[0.4,0.2,0.004] 
k1=K(1,1);
k2=K(1,2);
k3=K(1,3);

%Rt=1.966*10^10;
Rt=1.966*10^7;

heparin=iii;
NN=jjj;
tau=100;

%300nM IP3
%N=[0,0,1.1117*10^8,0,0];
N=[0,0,1.1117*10^5,0,0];

krel=40;
b=0.0000045;

if heparin == 1
    couleur='b';
else
    couleur='r';
end

C=0.22;
Cer(1,1)=500;

dt=1/(max(K))/100 %Time step
tf=300; %End time
dtsave=1/K(2)/200; %Time step for saving data. I took it twice smaller than k2 to be sure to record each opening.
tsave(:,1)=0:dtsave:tf; % Time vector for saving data
%Matrix for saving data 
Osave=zeros(1/dtsave*tf+1,1);
Rs=ones(size(Osave,1),1);
Os=zeros(size(Osave,1),1);
Ds=zeros(size(Osave,1),1);

for i=1:3 %Number of pulses
    for j=N(i)+1:N(i+1) %Select each IP3R bound to IP3 
        t=(i-1)*tf/4;
        %These parameters allow to consider the reaction in a
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
                
                PDR=k2*dt;
                Pnothing=1-PDR;
                if rand < PDR
                    O=0;
                    D=1;
                    td=t;
                end
            elseif D==1 
                PDR=k3*dt;
                Pnothing=1-PDR;
                if rand < PDR
                    R=1;
                    D=0;
                    tr=t;
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


for i=2:size(tsave,1)
    dt=tsave(i,1)-tsave(i-1,1);
    o4=Os(i,1);
    i;

    if tsave(i,1)<100
        krel=krel;
    else
        krel=krel*heparin;
    end
    Cer(i,1)=rg4dcer(dt,tau,krel, Cer(i-1,1),o4,Rt,b) ;


end
figure(jjj)
hold on
plot(tsave,(100-20)*(Cer./5/100)+20,couleur)
Cerbis=(100-20)*(Cer./5/100)+20;
xlabel('Time (s)')
ylabel('Ca^{2+} content (%)')
ylim([20 100]);

%to save data
if (iii == 0)
    T1 = table(tsave,Cerbis,'VariableNames',{ 'time', 'calcium'});
    writetable(T1, 'calcium_hep.txt');
else
    T1bis = table(tsave,Cerbis,'VariableNames',{ 'time', 'calcium'});
    writetable(T1bis, 'calcium.txt');
end

    end
end

    function y=dcer(Cer,tau,krel,o4,Rtot,b) 
    y=(500-Cer)/tau-krel*(Cer-0.22)*(b+o4/Rtot);
    end
    function y=rg4dcer(dt,tau,krel,Cer,o4,Rtot,b) 
    k1=dcer(Cer,tau,krel,o4,Rtot,b) ;
    k2=dcer(Cer+dt/2*k1,tau,krel,o4,Rtot,b) ;
    k3=dcer(Cer+dt/2*k2,tau,krel,o4,Rtot,b) ;
    k4=dcer(Cer+dt*k3,tau,krel,o4,Rtot,b) ;
    y=Cer+dt/6*(k1+2*k2+2*k3+k4);
    end


