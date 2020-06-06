%Minimalist model for quantal release (ROD). R->O, O->D and D->R are stochastic.
clear
clc

%Kinetic constants of ROD transitions
K=[0.4,0.2,0.004] ;
k1=K(1,1);
k2=K(1,2);
k3=K(1,3);
%Parameters of calcium evolution
krel=40;
b=0.0000045;
%Values for one cell (to get fast simulations)
Rt=49164;

%Number of tetraliganded receptors for 60, 120 and 220 nM IP3
N1=Rt*(60/(60+794))^4;
N2=Rt*(120/(120+794))^4;
N3=Rt*(220/(220+794))^4;

N=[0, 0, N1, N2, N3];

%Numeric parameters
dt=1/(max(K))/200 %Time step
tf=200; %End time
dtsave=1/K(2)/20;
tsave(:,1)=0:dtsave:tf; % Time vector for saving data
%Matrix for saving data 
Osave=zeros(1/dtsave*tf+1,1);
Rs=ones(size(Osave,1),1);
Os=zeros(size(Osave,1),1);
Ds=zeros(size(Osave,1),1);

for i=1:4 %Number of pulses
 for j=N(i)+1:N(i+1) %Select each IP3R bound to IP3
        t=(i-1)*tf/4;
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
            
            if (t-tsave(nsave))>dtsave
                nsave=nsave+1;
                Rs(nsave,1)=Rs(nsave,1)+R;
                Os(nsave,1)=Os(nsave,1)+O;
                Ds(nsave,1)=Ds(nsave,1)+D;
            end
        end
    end
end

figure(1)
plot(tsave, sum(Os,2))
hold on
xlabel('Time (s)')
ylabel('Number of open IP_{3}Rs')
figure(2)
plot(tsave, sum(Ds,2))
hold on
xlabel('Time (s)')
ylabel('Number of inactivated IP_{3}Rs')

R4=sum(Os,2);
t=tsave;
figure(3)
hold on

cf=0.1;   
cf(1)=500;

 for i=2:size(R4)
     tt(i)=t(i);
       cf(i)=exp(-krel*(b+R4(i)/Rt)*(t(i)-t(i-1))+log(cf(i-1)-0.22))+0.22   ;
 end
  
    cf=(100-20)*(cf./5/100)+20;
    
    plot(tt, cf)
    axis ([0 200 0 100])
    xlabel('Time (s)')
    ylabel('Ca^{2+} content in the ER (%)')
    
    %to save data
    tt2=tt.'; 
    ccf=cf.';
  
    T1 = table(tt2,ccf,'VariableNames',{ 'time', 'calcium'});
    writetable(T1, 'calcium.txt');
   
    T2=table(t,sum(Os,2),'VariableNames',{ 'time', 'open_rec'});
    writetable(T2, 'open.txt');
    
    T3=table(t,sum(Ds,2),'VariableNames',{ 'time', 'inact_rec'});
    writetable(T3, 'inact.txt');
    
 
  


 

    



   


