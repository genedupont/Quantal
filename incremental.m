%ROD scheme. R-O and O-D are reversible. D-R is irreversible
clear
clc 

%Kinetic constants of state transitions
k1=0.4;
k2=0.2;
k3=0.004;
km1=k1/100;
km2=k2/100;
K=[k1,k2,k3] ;
%
%kinetic constant for Ca2+
krel=8;   %krel/Rt~0.0008
b=0.0000001;

%Total number of receptors/cell
Rt=9101;

%Number of tetraliganded IP3 receptors after each IP3 addition
%N=[0,0,5,22,61];
%N=[0,0,10,44,122];
N=[0,0,20,88,244];

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

for i=1:4 
 for j=N(i)+1:N(i+1) 
        t=(i-1)*tf/4;
        nsave=(size(Osave,1)-1)/4*(i-1)+1;
        R=1;
        O=0;
        D=0;
        
        for step=1/dt*tf/4*(i-1)+1:1/dt*tf      
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
plot(tsave, sum(Os,2),'LineWidth',2)
hold on
xlabel('Time (s)','FontSize',14)
ylabel('Number of open IP_{3}Rs','FontSize',14)
figure(2)
plot(tsave, sum(Ds,2),'LineWidth',2)
hold on
xlabel('Time (s)','FontSize',14)
ylabel('Number of desensitised IP_{3}Rs','FontSize',14)

R4=sum(Os,2);
t=tsave;
figure(3)
hold on

cf=0.1;   
cf(1)=500;

 for i=2:size(R4)
     tt(i)=t(i);
       cf(i)=exp(-krel*(b+R4(i)/(Rt))*(t(i)-t(i-1))+log(cf(i-1)-0.22))+0.22   ;
 end
  
    cf=(100-20)*(cf./5/100)+20;
    
    plot(tt, cf,'LineWidth',2)
    axis ([0 200 0 100])
    xlabel('Time (s)','FontSize',14)
    ylabel('Ca^{2+} content in the ER (%)','FontSize',14)
    
    %To save the date in text files
    tt2=tt.'; 
    ccf=cf.';
  
    T1 = table(tt2,ccf,'VariableNames',{ 'time', 'calcium'});
    writetable(T1, 'calcium.txt');
   
    T2=table(t,sum(Os,2),'VariableNames',{ 'time', 'open_rec'});
    writetable(T2, 'open.txt');
    
    T3=table(t,sum(Ds,2),'VariableNames',{ 'time', 'inact_rec'});
    writetable(T3, 'inact.txt');
    
 
  


 

    



   


