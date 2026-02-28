close all;
clear
clc;
% long term aging experiment
% purposely distribute energy throughput on healthier cells to prolong pack
% lifetime
%% Parameters setup
celltype=1; %1: LFP; 2:LMO
para.Doubled_calendar_aging=0;%change to 1 to double calendar aging
para.Doubled_cyclic_aging=0;%change to 1 to double cyclic aging
knee_point_SOH=0.75;
n_cell=20; %number of cells
para.R0=0.01;%Internal resistance
more_diverse_cells=0;%change to 1 to to increase the varaince in the aging speed of different cells
more_fast_charging=0;%change to 1 to to increase the proportion of fast charging sessions
more_diverse_temperature=0;%change to 1 to to increase the variance of operating temperature
higher_operating_temperature=0;%change to 1 to increase the average operating temperature
end_of_life_SOH=0.7;%Battery cell and battery pack end-of-life SOH
%% Load Data
DATA = csvread("data_tesla.csv",1,0);
original_time=3600*24*100;
stepnum=length(DATA);
stepcSOC=[100;DATA(:,6)]/100;
stepdSOC=DATA(:,3)/100;
chargetime=DATA(:,5);
chargestarttime=DATA(:,1)+original_time;
chargeendtime=DATA(:,2)+original_time;
dischargestarttime=DATA(:,2)+original_time;
DCcurrent=4*ones(stepnum,1);
dischargecurrent=2*ones(stepnum,1);
stepT=25*ones(stepnum,1);%relaxation temperature
clear DATA
% clean the data
for i=1:stepnum
    if(stepcSOC(i)<stepdSOC(i)+0.001)
        stepcSOC(i)=stepdSOC(i)+0.001;
    end
end
%maximum charging time is 10h
for i=1:stepnum
    if(chargetime(i)>3600*10)
        chargetime(i)=3600*10;
        chargeendtime(i)=chargetime(i)+chargestarttime(i);
    end
end
for i=1:stepnum
    if(stepcSOC(i+1)>0.995)
        stepcSOC(i+1)=0.995;
    end
end
%% DC charging setup
averageCrate=zeros(stepnum,1);
for i=1:stepnum
    averageCrate(i)=(stepcSOC(i+1)-stepdSOC(i))/chargetime(i)*3600;
end
[~,idx]=sort(averageCrate,"descend");
if more_fast_charging
    DCpercentage=0.5;
else
    DCpercentage=0.318;
end

num=DCpercentage*stepnum;
whether_DC=zeros(stepnum,1);
for i=1:num
    whether_DC(idx(i))=1;
end
%% Other Parameters
if(celltype==2)
    multiplier=4;%repeat the aging profile for 4 times
else
    multiplier=16;%repeat the aging profile for 16 times
end
Uphase=2.5*n_cell;%phase voltage
m=6;%SOC division number during CV stage
%% battery initial states generation
rng(100)%set random seed 
if more_diverse_cells
    gamma=1+0.15*randn(n_cell,1);
else
    gamma=1+0.1*randn(n_cell,1);
end
std=0.001;
Qmax=2.3;
Qmax1=Qmax.*(1+std*randn(n_cell,1));
if more_diverse_temperature
    T_std=4;
else
    T_std=2;
end
if higher_operating_temperature
    T_mean_add=20; %compared with the relaxation temperature
else
    T_mean_add=10;
end
deltaT=T_mean_add+T_std*randn(n_cell,1);
Ahrecord1=zeros(1,stepnum);
Ahrecord2=zeros(1,stepnum);
%% SOC balancing
Iremember=zeros(stepnum,1);
[SOH_algorithm1, SOH_algorithm1_order, SOH_loss_cyc1, SOH_loss_cal1, Lifetime1, Iremember, iplot]=simulation_main(celltype,n_cell,1,m,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,stepcSOC,stepdSOC,Qmax1,Qmax,knee_point_SOH,end_of_life_SOH,0,Iremember,0,gamma,multiplier,para);
%% Proposed method (no noise)
[SOH_algorithm2, SOH_algorithm2_order, SOH_loss_cyc2, SOH_loss_cal2, Lifetime2, ~, ~]=simulation_main(celltype,n_cell,n_cell,m,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,stepcSOC,stepdSOC,Qmax1,Qmax,knee_point_SOH,end_of_life_SOH,1,Iremember,0,gamma,multiplier,para);
%% Proposed method (2% SOH noise)
SOH_noise=0.02;
[SOH_algorithm3, SOH_algorithm3_order, SOH_loss_cyc3, SOH_loss_cal3, Lifetime3, ~, ~]=simulation_main(celltype,n_cell,n_cell,m,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,stepcSOC,stepdSOC,Qmax1,Qmax,knee_point_SOH,end_of_life_SOH,1,Iremember,SOH_noise,gamma,multiplier,para);
%% Figures
limit=ceil((multiplier*Lifetime2*1.01)/365/3600/24/5)*5;
f1=figure;
tiledlayout(1,3,"TileSpacing","tight","Padding","compact");
f1.Position = [100 100 1200 540];
ax = nexttile;
hold(ax,'on')
nLines  = 15;
C = parula(nLines);
%C = turbo(nLines);
colororder(ax, C);       % set the axes' color cycle
ax.ColorOrderIndex = 1;  % start from first color

%x=(0:stepnum)*multiplier;
x = (chargeendtime(1:stepnum)-chargeendtime(1))*multiplier;
y = SOH_algorithm1_order(1,:);
%plot(x(1:iplot),y(1:iplot))
for j=1:10  
    y=SOH_algorithm1_order(n_cell/10*j,:);
    plot(x(1:iplot)/365/24/3600,y(1:iplot)*100,DisplayName=append(int2str(n_cell/10*j),"th highest SOH"))
end
xlabel('Time (years)','FontSize',14)
ylabel('SOH (%)','FontSize',14)
set(gca,'Fontsize',14)
legend(Location="best",FontSize=12)
ylim([0.6 1]*100)
xlim([0 limit])
grid on
title(append('(i) Baseline method (SOC balancing)'),'FontSize',13.5,'FontWeight', 'normal');
%disp(append('Strategy 1 mean = ', num2str(mean(Ahtotal1))))
%disp(append('Strategy 1 min = ', num2str(min(Ahtotal1))))
%disp(append('Strategy 1 max = ', num2str(max(Ahtotal1))))
%disp(append('Strategy 1 variance = ', num2str(var(Ahtotal1))))
ax = nexttile;
hold(ax,'on')
nLines  = 15;
C = parula(nLines);
%C = turbo(nLines);
colororder(ax, C);       % set the axes' color cycle
ax.ColorOrderIndex = 1;  % start from first color

%plot((0:stepnum)*multiplier,SOH_algorithm2_order(1,:))
for j=1:10
    plot([x;0]/365/24/3600,SOH_algorithm2_order(n_cell/10*j,:)*100,DisplayName=append(int2str(n_cell/10*j),"th highest SOH"))
end
xlabel('Time (years)','FontSize',14)
%ylabel('SOH [%]','FontSize',12)
set(gca,'Fontsize',14)
legend(Location="best",FontSize=12)
ylim([0.6 1]*100)
xlim([0 limit])
grid on
%disp(append('Strategy 2 mean = ', num2str(mean(Ahtotal2))))
%disp(append('Strategy 2 min = ', num2str(min(Ahtotal2))))
%disp(append('Strategy 2 max = ', num2str(max(Ahtotal2))))
%disp(append('Strategy 2 variance = ', num2str(var(Ahtotal2))))
title(append('(ii) Our method with no noise (+',num2str(100*(Lifetime2-Lifetime1)/Lifetime1,'%.0f '),'% RUL)'),'FontSize',13.5,'FontWeight', 'normal');
disp(append("Improvement in percentage: ",num2str(100*(Lifetime2-Lifetime1)/Lifetime1)))
ax = nexttile;
hold(ax,'on')
nLines  = 15;
C = parula(nLines);
%C = turbo(nLines);
colororder(ax, C);       % set the axes' color cycle
ax.ColorOrderIndex = 1;  % start from first color
%x=(0:stepnum)*multiplier;
for j=1:10
    plot([x;0]/365/24/3600,SOH_algorithm3_order(n_cell/10*j,:)*100,DisplayName=append(int2str(n_cell/10*j),"th highest SOH"))
end
xlabel('Time (years)','FontSize',14)
set(gca,'Fontsize',14)
legend(Location="best",FontSize=12)
ylim([0.6 1]*100)
xlim([0 limit])
grid on
title(append('(iii) Our method with 2% noise (+',num2str(100*(Lifetime3-Lifetime1)/Lifetime1,'%.0f '),'% RUL)'),'FontSize',13.5,'FontWeight', 'normal');
disp(append("Algorithm 3 improvement in percentage: ",num2str(100*(Lifetime3-Lifetime1)/Lifetime1)))

%exportgraphics(f1,append('comparison',int2str(celltype),'.png'),'Resolution',900)
function OCV_o = OCV1(SOC)
    OCV_o=OCV(SOC, 1);
end
function OCV_o = OCV2(SOC)
    OCV_o=OCV(SOC, 2);
end


function OCV = OCV(SOC, celltype)
if(celltype==1)
    a1 = -5.863e-1;
    a2 = 21.9;
    a3 = 3.414;
    a4 = 1.102e-1;
    a5 = -1.718e-1;
    a6 = 8e-3;
    OCV = a1 * exp(-a2 * SOC) + a3 + a4 * SOC + a5 * exp( -a6 / (1-SOC));
else
    s = SOC;
    if(s<=0)
        s=0.01;
    end
    a=3.875;
    b=-0.335;
    c=-0.5332;
    d=0.8315;
    m=0.653;
    n=0.6;
    OCV = a + b .* (-log(s)).^m + c .* s + d .* exp(n .* (s - 1));
end
end
function SOC = I_CV_inverse(I_CV)
SOC=(2.6963-I_CV)/2.58;
if(SOC>1)
    SOC=1;
end
end
function Cmax = I_CV(SOC) %%%%
Cmax = -2.58*SOC + 2.6963;
end

function [fc_sum_new, Fadedelta]=agingmodelLMO(t,deltat,T,SOC,DOD,fc_sum,SOHnow,SOHknee,increase,multiplier)
if(abs(DOD)>0.01)
    kd1=1.4e5;
    kd2=-0.501;
    kd3=-1.23e5;
    ks=1.04;
    sigma=0.5;
    kT=0.0693;
    Tref=25;
    ST=exp(kT*(T-Tref)*(Tref+273.15)/(T+273.15));
    Ssig=exp(ks*(SOC-sigma));
    Sdelta=1/(kd1*DOD^kd2+kd3);
    fc_sum_new=fc_sum+multiplier*Sdelta*Ssig*ST/2;
else
    fc_sum_new=fc_sum;
end
Fadedelta=agingLMO(T,SOC,(t+deltat)*multiplier,fc_sum_new)-agingLMO(T,SOC,t*multiplier,fc_sum);
Fadedelta=Fadedelta*100;
if(SOHnow<SOHknee)
    Fadedelta=Fadedelta*(1+(SOHknee-SOHnow)/(0.01)*increase);
end
end

function Fade=agingLMO(T,SOC,t,fc_sum)
alpha=0.0575;
beta=121;
ks=1.04;
sigma=0.5;
kT=0.0693;
Tref=25;
kt=4.14e-10;
ST=exp(kT*(T-Tref)*(Tref+273.15)/(T+273.15));
Ssig=exp(ks*(SOC-sigma));
St=kt*t;
ft=St*Ssig*ST;
fd=fc_sum+ft;
Fade=1-(alpha*exp(-beta*fd)+(1-alpha)*exp(-fd));
end
function [fc_new, Fadecyc]=agingcyc(t,SOC,DOD,Qmax,deltaAh,T,crate,celltype,SOHnow,SOHknee,increase,multiplier,fc_sum,para)
if(celltype==1)
    Tem=T+273.15;
    a = 2.0916e-8;
    b = -1.2179e-5;
    c = 0.0018;
    d = -1.7082e-6;
    e = 0.0556;
    k_crate=(a*Tem^2+b*Tem+c)*exp((d*Tem+e)*crate);
    Fadecyc=k_crate*deltaAh;
    if(SOHnow<SOHknee)
        Fadecyc=Fadecyc*(1+(SOHknee-SOHnow)/(0.01)*increase);
    end
    Fadecyc=Fadecyc*multiplier/Qmax*100/2;
    fc_new=0;
else
    [fc_new, Fadecyc]=agingmodelLMO(t,0,T,SOC,abs(DOD),fc_sum,SOHnow,SOHknee,increase,multiplier);
    Fadecyc=Fadecyc/2;
end
if para.Doubled_cyclic_aging
    Fadecyc=Fadecyc*2;
end
end

function Fadecal=agingcal(t,deltat,T,SOC,celltype,Qmax,SOHnow,SOHknee,increase,multiplier,fc_sum,para)
if(celltype==1)
    Tem=T+273.15;
    f = 5.9808e6;
    g = 0.6898;
    h = -6.4647e3;
    day=t/3600/24;
    deltaday=deltat/3600/24;
    k_T=exp(h/Tem);
    k_SOC=f*exp(g*SOC);
    Fadecal=k_T*k_SOC*deltaday/2/sqrt(day*multiplier)*multiplier;
    if(SOHnow<SOHknee)
        Fadecal=Fadecal*(1+(SOHknee-SOHnow)/(0.01)*increase);
    end
    Fadecal=Fadecal*100/Qmax;
else
    DOD=0;
    [~,Fadecal]=agingmodelLMO(t*multiplier,deltat*multiplier,T,SOC,DOD,fc_sum,SOHnow,SOHknee,0.5,1);
end
if para.Doubled_calendar_aging
    Fadecal=Fadecal*2;
end
end


function [SOH_algorithm3, SOH_algorithm3_order, SOH_loss_cyc, SOH_loss_cal, Lifetime3, Irememberout, iplot]=simulation_main(celltype,n1,n2,m,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,stepcSOC,stepdSOC,Qmax1,Qmax,SOHknee,Pack_SOH_end,known_current,Iremember,SOH_noise,gamma,multiplier,para)
rng(100)%set random seed 
ttt=[];
fc_sum=zeros(n1,1);
increase=0.5;
Ahtotal3=100*ones(n1,1);
R0=para.R0;%Internal resistance
Lifetime3=0;
Qnow=Qmax1;
Ahrecord3=zeros(1,stepnum);
SOH_loss_cyc=zeros(n1,stepnum);
SOH_loss_cal=zeros(n1,stepnum);
SOH_algorithm3=NaN(n1,stepnum+1);
SOH_algorithm3_order=NaN(n1,stepnum+1);
SOH_algorithm3(:,1)=Qmax1/Qmax;
SOH_algorithm3_order(:,1)=sort(SOH_algorithm3(:,1),'descend');
Irememberout=Iremember;
iplot=stepnum;
timesremember=1;
for i=1:stepnum
    SOH_algorithm3(:,i+1)=SOH_algorithm3(:,i);
    if(i>1)
        SOH_loss_cal(:,i)=SOH_loss_cal(:,i-1);
        SOH_loss_cyc(:,i)=SOH_loss_cyc(:,i-1);
    end
    %discharge aging
    d=real(acos(OCV(stepdSOC(i),celltype)*(0:1:n1)/Uphase))/pi*2;
    %optimization
    duty2=zeros(1,n1);
    for ii=1:n1
        duty2(ii)=0.5*(d(ii)+d(ii+1));
    end
    totalAh=sum(Qnow)-sum(SOH_algorithm3(:,i+1)*Qmax)*stepdSOC(i);
    deltaAh=zeros(n1,1);
    num=2000;
    if(n2==1)
        deltaAh=totalAh*SOH_algorithm3(:,i+1)/sum(SOH_algorithm3(:,i+1));
    else
        for ii=1:num
            [~,idx]=sort(Qnow-deltaAh,'descend');
            for jj=1:n1
                deltaAh(idx(jj))=deltaAh(idx(jj))+totalAh/num*duty2(jj)/sum(duty2);
            end
        end
    end
    Ahtotal3=Ahtotal3+deltaAh;%*multiplier;
    Ahrecord3(i)=totalAh;
    %aging
    for ii=1:n1
        DOD=deltaAh(ii)/SOH_algorithm3(ii,i+1)/Qmax;
        SOCtemp=(Qnow(ii)*2-deltaAh(ii))/SOH_algorithm3(ii,i+1)/Qmax/2;

        [fc_sum(ii), cyc_loss]=agingcyc(dischargestarttime(i),SOCtemp,DOD,Qmax,deltaAh(ii),stepT(i)+deltaT(ii),dischargecurrent(i)/Qmax,celltype,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para);
        cyc_loss=cyc_loss/100*gamma(ii);
        SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
        SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;
        if(i>1)
            %assume that the temperature is always 25 deg C for calendar
            cal_loss=agingcal(dischargestarttime(i),chargestarttime(i)-dischargestarttime(i-1),stepT(i),SOCtemp,celltype,Qmax,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para)/100*gamma(ii);
            SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
            SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
        end
    end
    %Qnow=Qnow-deltaAh;
    %SOCnow=Qnow./Qmax./SOH_algorithm3(:,i);
    SOCnow=ones(n1,1)*stepdSOC(i);
    Qnow=SOCnow.*Qmax.*SOH_algorithm3(:,i+1);
    %charge aging   
    if(whether_DC(i)==1)
        %DC
        deltaAh=stepcSOC(i+1)*SOH_algorithm3(:,i+1)*Qmax-Qnow;
        DCrate=max(deltaAh)/(chargeendtime(i)-chargestarttime(i))/Qmax;
        for ii=1:n1
           SOCtemp=(stepcSOC(i+1)+Qnow(ii)/SOH_algorithm3(ii,i+1)/Qmax)/2;
           DOD=stepcSOC(i+1)-Qnow(ii)/SOH_algorithm3(ii,i+1)/Qmax;
           [fc_sum(ii), cyc_loss]=agingcyc(chargestarttime(i),SOCtemp,DOD,Qmax,deltaAh(ii),stepT(i)+deltaT(ii),DCrate,celltype,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para);
           cyc_loss=cyc_loss/100*gamma(ii);
           SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
           SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;
           %calendar
           cal_loss=agingcal(chargestarttime(i),chargeendtime(i)-chargestarttime(i),stepT(i)+deltaT(ii),SOCtemp,celltype,Qmax,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para)/100*gamma(ii);
           SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
           SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
           if(i<stepnum)
               SOCtemp=(stepcSOC(i+1)+Qnow(ii)/SOH_algorithm3(ii,i+1)/Qmax)/2;
               cal_loss=agingcal(chargeendtime(i),dischargestarttime(i+1)-chargeendtime(i),stepT(i),SOCtemp,celltype,Qmax,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para)/100*gamma(ii);
               SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
               SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
           end
        end
        Qnow=stepcSOC(i+1)*SOH_algorithm3(:,i+1)*Qmax;

    else
        %AC
        %calculate Qnow by optimization
        %record also the c-rates
        %calculate the duty cycles
        SOH_algorithm3(:,i+1)=max(0.1,SOH_algorithm3(:,i+1));
        SOH_now_noisy=SOH_algorithm3(:,i+1)+SOH_noise*randn(n1,1);

        deltaSOH=NaN(n1,1);
        for iii=1:n1
            deltaSOH(iii)=max(SOH_now_noisy(iii)-Pack_SOH_end,0.001);
        end
        % Use maximum current, find the maximum voltage that makes charging
        % possible
        Tlimit_modified=chargetime(i)/3600;
        Unow = n1*OCV((stepdSOC(i)+stepcSOC(i+1))/2,celltype);
        Q_opt=[];
        Iavg=min(I_CV(0.1), I_CV(stepdSOC(i)));
        params.SOH_EOL=Pack_SOH_end;
        while(known_current==1 && isempty(Q_opt) && Unow>4)
            SOC_CV=I_CV_inverse(Iavg);
            if(SOC_CV>=1)
                SOCgrid=ones(1,m+1);
            else
                SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
            end
            SOCgrid=max(SOCgrid,stepdSOC(i));
            SOCgrid2=[stepdSOC(i) SOCgrid];
            SOCs=zeros(1,m+1);
            for ii=1:m+1
                SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
            end
            Is=I_CV(SOCs);
            Is(1)=Iavg;
            if celltype==1
                [Q_opt, ~, ~] = solve_charge_allocation_lp( ...
                    SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), ...
                    stepcSOC(i+1)*sum(SOH_algorithm3(:,i+1))*Qmax, Is, SOCgrid, ...
                    Unow, 'sinusoidal', @OCV1, R0, Qmax, Tlimit_modified, ...
                    Uphase, dischargecurrent(i)/Qmax, params);
            else
                [Q_opt, ~, ~] = solve_charge_allocation_lp( ...
                    SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), ...
                    stepcSOC(i+1)*sum(SOH_algorithm3(:,i+1))*Qmax, Is, SOCgrid, ...
                    Unow, 'sinusoidal', @OCV2, R0, Qmax, Tlimit_modified, ...
                    Uphase, dischargecurrent(i)/Qmax, params);
            end
            Unow = Unow*0.95;
        end

        % Fix the voltage, find best current current
        %Find minimum current (C rate)
        if known_current && Iremember(i)>0
            Iavgmin=Iremember(i);
        else
            d=real(acos(OCV(0.5*stepdSOC(i)+0.5*stepcSOC(i+1),celltype)*(0:1:n1)/Unow))/pi*2;
            %optimization
            duty=zeros(1,n1);
            for ii=1:n1
                duty(ii)=0.5*(d(ii)+d(ii+1));
            end
            Iavgmin=sum(SOH_algorithm3(:,i+1))*(stepcSOC(i+1)-stepdSOC(i))/sum(duty)/Tlimit_modified;
        end
        % search for optimal current
        best=inf;
        Ahbest=[];
        if celltype==1
            while Iavg>Iavgmin*0.7
                SOC_CV=I_CV_inverse(Iavg);
                if(SOC_CV>=1)
                    SOCgrid=ones(1,m+1);
                else
                    SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
                end
                SOCgrid=max(SOCgrid,stepdSOC(i));
                SOCgrid2=[stepdSOC(i) SOCgrid];
                SOCs=zeros(1,m+1);
                for ii=1:m+1
                    SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
                end
                Is=I_CV(SOCs);
                Is(1)=Iavg;
                if n2==1 % SOC balancing
                    [Q_opt2, sol, ~] = solve_charge_allocation_equalSOC_lp( ...
                        SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), stepcSOC(i+1), ...
                        Is, SOCgrid, ...
                        Unow, 'sinusoidal', @OCV1, R0, Qmax, Tlimit_modified, params);
                else % SOH balancing
                    [Q_opt2, sol, ~] = solve_charge_allocation_lp( ...
                        SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), ...
                        stepcSOC(i+1)*sum(SOH_algorithm3(:,i+1))*Qmax, Is, SOCgrid, ...
                        Unow, 'sinusoidal', @OCV1, R0, Qmax, Tlimit_modified, ...
                        Uphase, dischargecurrent(i)/Qmax, params);
                end
                if sol.fval<best
                    best=sol.fval;
                    Ahbest=Q_opt2;
                    Iavgbest=Iavg;
                    SOCgridbest=SOCgrid;
                end
                Iavg = Iavg*0.95;
            end
        else
            while Iavg>Iavgmin*0.7
                SOC_CV=I_CV_inverse(Iavg);
                if(SOC_CV>=1)
                    SOCgrid=ones(1,m+1);
                else
                    SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
                end
                SOCgrid=max(SOCgrid,stepdSOC(i));
                SOCgrid2=[stepdSOC(i) SOCgrid];
                SOCs=zeros(1,m+1);
                for ii=1:m+1
                    SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
                end
                Is=I_CV(SOCs);
                Is(1)=Iavg;
                if n2==1 % SOC balancing
                    [Q_opt2, sol, ~] = solve_charge_allocation_equalSOC_lp( ...
                        SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), stepcSOC(i+1), ...
                        Is, SOCgrid, ...
                        Unow, 'sinusoidal', @OCV2, R0, Qmax, Tlimit_modified, params);
                else % SOH balancing
                    [Q_opt2, sol, ~] = solve_charge_allocation_lp( ...
                        SOH_now_noisy, Qnow*0.99, Qmax.*SOH_algorithm3(:,i+1), ...
                        stepcSOC(i+1)*sum(SOH_algorithm3(:,i+1))*Qmax, Is, SOCgrid, ...
                        Unow, 'sinusoidal', @OCV2, R0, Qmax, Tlimit_modified, ...
                        Uphase, dischargecurrent(i)/Qmax, params);
                end
                if sol.fval<best
                    best=sol.fval;
                    Ahbest=Q_opt2;
                    Iavgbest=Iavg;
                    SOCgridbest=SOCgrid;
                end
                Iavg = Iavg*0.95;
            end
        end
        if(isempty(Ahbest))
            disp(append('Failed: ',int2str(i)))
        else
            if(known_current==0)
                Irememberout(i)=Iavgbest;
            end
            timesremember=Iavgbest/Iavgmin;
        end
        SOC_copy=SOCnow;
        %AC charging cyclic aging
        SOH_memory=SOH_algorithm3(:,i+1);
        for ii=1:n1
            for iii=1:m+1
                %cyclic
                deltaAh=Ahbest(ii,iii);
                SOCtemp=(2*SOCnow(ii)+deltaAh/Qmax/SOH_memory(ii))/2;
                SOCnow(ii)=SOCnow(ii)+deltaAh/Qmax/SOH_memory(ii);
                DOD=deltaAh/Qmax/SOH_memory(ii);
                [fc_sum(ii), cyc_loss]=agingcyc(chargestarttime(i),SOCtemp,DOD,Qmax,deltaAh,stepT(i)+deltaT(ii),dischargecurrent(i)/Qmax,celltype,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para);
                cyc_loss=cyc_loss/100*gamma(ii);
                SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
                SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;
            end
        end
        %AC charging calendar aging
        for ii=1:n1
            deltat=chargeendtime(i)-chargestarttime(i);
            avgSOC=SOC_copy(ii)+sum(Ahbest(ii,:))/(SOH_memory(ii)*Qmax)/2;
            cal_loss=agingcal(chargestarttime(i),deltat,stepT(i)+deltaT(ii),avgSOC,celltype,Qmax,SOH_algorithm3(ii,i+1),SOHknee,increase,multiplier,fc_sum(ii),para)/100*gamma(ii);
            SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
            SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
            if(i<stepnum)
                avgSOC=SOC_copy(ii)+sum(Ahbest(ii,:))/(SOH_memory(ii)*Qmax);
                deltat=dischargestarttime(i+1)-chargeendtime(i);
                cal_loss=agingcal(chargeendtime(i),deltat,stepT(i),avgSOC,celltype,Qmax,SOH_memory(ii),SOHknee,increase,multiplier,fc_sum(ii),para)/100*gamma(ii);
                SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
                SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
            end
        end
        Qnow = SOCnow.*SOH_memory.*Qmax;
        clear Ahbest SOC_copy SOH_memory
    end
    ttt=[ttt std(SOCnow) / mean(SOCnow)];

    SOH_algorithm3_order(:,i+1)=sort(SOH_algorithm3(:,i+1),'descend');
    SOHtemp=SOH_algorithm3(:,i+1);
    for ii=1:n1
        if(SOHtemp(ii)<Pack_SOH_end)
            SOHtemp(ii)=0;
        end
    end
    disp(append(num2str(timesremember)," (",int2str(i)," cycles completed)"))
    if (mean(SOHtemp)<Pack_SOH_end && Lifetime3==0)
        %Lifetime3 = sum(Ahrecord3);%%%
        Lifetime3 = chargeendtime(i)-dischargestarttime(1);%%%
        if(known_current==0)
            iplot=i+1;
%             for ii =1:n1
%                 SOH_algorithm3(ii,i+1) = SOH_algorithm3(ii,1);
%             end
        else
            iplot=1;
            break;
        end
    end
    if (min(SOH_algorithm3(:,i+1))<0 && Lifetime3~=0)
        break;
    end
end
end

function [Q_opt, sol, prob] = solve_charge_allocation_lp( ...
    SOH, Q_init, Q_max, Q_final_sum, I_no_avg, SOC_max, ...
    U_phase, modulationType, fOCV, R0, Q_no, t_total, ...
    U_dis_phase, I_no_dis_avg, params)
%SOLVE_CHARGE_ALLOCATION_LP  Solve the LP in (cf) with constraints 1)–8).
%
% Decision variables:
%   Q(i,j)   = stage-wise added charge (Ah) for cell i in stage j
%              i = 1..n1, j = 0..m   (MATLAB uses j = 1..(m+1))
%   t(k,j)   = epigraph auxiliary (free) for top-k sum in stage j
%   s(i,k,j) = epigraph slack (>=0) for top-k sum in stage j
%
% ------------------------- Inputs (required) -------------------------
% SOH            [n1x1]  cell SOH values (either 0–1 or 0–100; auto-detect)
% Q_init         [n1x1]  initial charge contents (Ah) of each cell
% Q_max          [n1x1]  current max capacity (Ah) of each cell
% Q_final_sum    [1x1]   target total charge content at end (Ah), sum_i Q_final,i
% I_no_avg       [1x(nStages)] normalized avg C-rate per stage (C-rate), nStages=m+1
% SOC_max        [1x(nStages)] per-stage SOC upper bound at end of stage k (0–1)
% U_phase        [1x1]   fixed phase-voltage magnitude (V)
% modulationType 'sinusoidal' or 'dc'
% fOCV           function handle, U = fOCV(SOC) where SOC in [0,1]
% R0             [1x1]   average source-side resistance (Ohm)
% Q_no           [1x1]   nominal capacity (Ah) used to map C-rate -> current
% t_total        [1x1]   charging time limit (hours)
% U_dis_phase    [1x1]   representative discharge phase-voltage magnitude (V)
% I_no_dis_avg   [1x1]   representative normalized discharge C-rate (C-rate)
%
% params (optional struct) fields:
%   .SOH_EOL  (default 0.70) EOL threshold
%   .kappa    (default 0.1)
%   .epsSOH   (default 1e-3)
%   .Mbig     (default 1e6)
%   .linprogOptions (default dual-simplex, Display=none)
%
% ------------------------- Outputs -------------------------
% Q_opt   [n1 x nStages] optimal stage-wise allocations Q(i,j) in Ah
% sol     struct with fields:
%   .exitflag, .fval, .output, .lambda, .x
%   .Q_final   [n1x1] final stored charge contents (Ah) = Q_init + sum_j Q_opt(:,j)
%   .weights   [n1 x nStages] w(i,j)
% prob    struct with useful intermediate quantities:
%   .d          [n1 x nStages] stage duty cycles d(i,j)
%   .dprime     [n1x1] discharge duty pattern d'(i) (level-indexed)
%   .U_ter      [1 x nStages] terminal-voltage proxy per stage
%   .Aineq,.bineq,.Aeq,.beq,.lb,.ub,.f  (LP matrices)
%
% Requires: Optimization Toolbox (linprog)

% -------------------- defaults & checks --------------------
if nargin < 16 || isempty(params), params = struct(); end
if ~isfield(params,'SOH_EOL'),  params.SOH_EOL = 0.70; end
if ~isfield(params,'kappa'),    params.kappa   = 0.1; end
if ~isfield(params,'epsSOH'),   params.epsSOH  = 1e-3; end
if ~isfield(params,'Mbig'),     params.Mbig    = 1e6; end
if ~isfield(params,'linprogOptions')
    params.linprogOptions = optimoptions('linprog', ...
        'Algorithm','dual-simplex', 'Display','none');
end

SOH      = SOH(:);
Q_init   = Q_init(:);
Q_max    = Q_max(:);

n1 = numel(SOH);
nStages = numel(I_no_avg);  % = m+1

assert(numel(Q_init)==n1 && numel(Q_max)==n1, 'SOH, Q_init, Q_max must have same length n1.');
assert(numel(SOC_max)==nStages, 'SOC_max must have length equal to numel(I_no_avg) = m+1.');
assert(isa(fOCV,'function_handle'), 'fOCV must be a function handle.');

% Auto-detect SOH scale (0–100 vs 0–1)
if max(SOH) > 1.5
    SOH = SOH/100;
end
SOH_EOL = params.SOH_EOL;
if SOH_EOL > 1.5
    SOH_EOL = SOH_EOL/100;
end

% -------------------- compute weights w(i,j) --------------------
kappa  = params.kappa;
epsSOH = params.epsSOH;
Mbig   = params.Mbig;

w = zeros(n1, nStages);
for j = 1:nStages
    factor = (1 + kappa*I_no_avg(j));
    for i = 1:n1
        if SOH(i) > SOH_EOL + epsSOH
            w(i,j) = factor / (SOH(i) - SOH_EOL)^2;
        else
            w(i,j) = factor * Mbig;
        end
    end
end

% -------------------- compute stage duty cycles d(i,j) --------------------
Q_init_sum = sum(Q_init);
SOC_avg = (Q_final_sum + Q_init_sum) / (2*sum(Q_max));
SOC_avg = min(max(SOC_avg, 0), 1);

U_ocv = fOCV(SOC_avg);
U_ter = zeros(1,nStages);
for j = 1:nStages
    U_ter(j) = U_ocv + (Q_no * I_no_avg(j)) * R0; % V
end

d = zeros(n1, nStages);
switch lower(string(modulationType))
    case "sinusoidal"
        for j = 1:nStages
            for i = 1:n1
                ratio = ((2*i-1) * U_ter(j)) / (2*U_phase);
                ratio = min(ratio, 1); % per paper
                d(i,j) = (2/pi) * acos(ratio);
            end
        end
    case "dc"
        for j = 1:nStages
            for i = 1:n1
                if (i-1)*U_ter(j) > U_phase
                    d(i,j) = 0;
                elseif i*U_ter(j) < U_phase
                    d(i,j) = 1;
                else
                    d(i,j) = (U_phase - (i-1)*U_ter(j)) / U_ter(j);
                end
            end
        end
    otherwise
        error('modulationType must be ''sinusoidal'' or ''dc''.');
end

sumd = sum(d,1);
if any(sumd <= 0)
    error('Some stages have sum_i d(i,j) = 0, making the time constraint ill-defined.');
end

% gamma(k,j) = (sum_{l=1}^k d(l,j)) / (sum_{l=1}^{n1} d(l,j))
gamma = zeros(n1, nStages);
for j = 1:nStages
    csum = cumsum(d(:,j));
    gamma(:,j) = csum ./ sumd(j);
end

% -------------------- discharge duty pattern d'(i) --------------------
SOC_dis = Q_final_sum / (2*sum(Q_max));
SOC_dis = min(max(SOC_dis, 0), 1);
U_ocv_dis = fOCV(SOC_dis);
U_dis_ter = U_ocv_dis - (Q_no * I_no_dis_avg) * R0;

dprime = zeros(n1,1);
for i = 1:n1
    ratio = ((2*i-1) * U_dis_ter) / (2*U_dis_phase);
    ratio = min(ratio, 1);
    dprime(i) = (2/pi) * acos(ratio);
end
sumdprime = sum(dprime);
if sumdprime <= 0
    error('Discharge duty pattern has sum(dprime)=0; check U_dis_phase and U_dis_ter.');
end

% -------------------- variable indexing --------------------
% Q(i,j): i=1..n1, j=1..nStages  -> index (j-1)*n1 + i
% t(k,j): k=1..n1, j=1..nStages  -> index NQ + (j-1)*n1 + k
% s(i,k,j): i=1..n1,k=1..n1,j=1..nStages
%          -> index NQ+Nt + (j-1)*n1*n1 + (k-1)*n1 + i
NQ = n1*nStages;
Nt = n1*nStages;
Ns = n1*n1*nStages;
N  = NQ + Nt + Ns;

idxQ = @(i,j) (j-1)*n1 + i;
idxT = @(k,j) NQ + (j-1)*n1 + k;
idxS = @(i,k,j) NQ + Nt + (j-1)*n1*n1 + (k-1)*n1 + i;

% Objective
f = zeros(N,1);
f(1:NQ) = w(:);

% Bounds
lb = -inf(N,1);
ub = inf(N,1);

% Q >= 0
lb(1:NQ) = 0;

% s >= 0
lb(NQ+Nt+1:end) = 0;

% t free (lb stays -inf)

% -------------------- equality: total added charge --------------------
deltaQ_total = Q_final_sum - Q_init_sum;
Aeq = sparse(1, N);
Aeq(1,1:NQ) = 1;
beq = deltaQ_total;

% -------------------- build inequalities Aineq*x <= bineq --------------------
% Row blocks:
%  (a) time limit: 1 row
%  (b) per-stage SOC upper bounds: n1*nStages rows
%  (c) epigraph constraints: n1*n1*nStages rows
%  (d) strong_lp constraints: n1*nStages rows
%  (e) SOH ordering (adjacent in sorted SOH): (n1-1) rows
%  (f) discharge utilization (k=1..n1-1): (n1-1) rows

nTime   = 1;
nSOC    = n1*nStages;
nEpi    = n1*n1*nStages;
nStrong = n1*nStages;
nOrder  = max(n1-1,0);
nDis    = max(n1-1,0);

nIneq = nTime + nSOC + nEpi + nStrong + nOrder + nDis;

% Estimate nnz for preallocation
nnz_time   = NQ;
nnz_soc    = n1 * (nStages*(nStages+1)/2);
nnz_epi    = 3 * nEpi;
nnz_strong = (2*n1 + 1) * nStrong;
nnz_order  = 2 * nStages * nOrder;
nnz_dis    = NQ * nDis;
nnz_est    = nnz_time + nnz_soc + nnz_epi + nnz_strong + nnz_order + nnz_dis;

rows = zeros(nnz_est,1);
cols = zeros(nnz_est,1);
vals = zeros(nnz_est,1);
bineq = zeros(nIneq,1);
p = 0;   % pointer into (rows,cols,vals)
r = 0;   % row counter

% (a) time limit
r = r + 1;
for j = 1:nStages
    cj = 1 / (Q_no * I_no_avg(j) * sumd(j));  % hours per Ah allocated in stage j
    base = (j-1)*n1;
    idxs = base + (1:n1);
    rows(p+(1:n1)) = r;
    cols(p+(1:n1)) = idxs;
    vals(p+(1:n1)) = cj;
    p = p + n1;
end
bineq(r) = t_total;

% (b) per-stage SOC upper bounds
% Q_init(i) + sum_{j=1..k} Q(i,j) <= Q_max(i)*SOC_max(k)
for i = 1:n1
    for k = 1:nStages
        r = r + 1;
        % coefficients for Q(i,1..k)
        stageIdxs = (0:(k-1))*n1 + i;  % indices in Q block
        nn = numel(stageIdxs);
        rows(p+(1:nn)) = r;
        cols(p+(1:nn)) = stageIdxs;
        vals(p+(1:nn)) = 1;
        p = p + nn;

        bineq(r) = Q_max(i)*SOC_max(k) - Q_init(i);
    end
end

% (c) epigraph constraints: Q(i,j) - t(k,j) - s(i,k,j) <= 0
for j = 1:nStages
    for k = 1:n1
        for i = 1:n1
            r = r + 1;
            rows(p+1) = r; cols(p+1) = idxQ(i,j); vals(p+1) = 1;   p = p+1;
            rows(p+1) = r; cols(p+1) = idxT(k,j); vals(p+1) = -1;  p = p+1;
            rows(p+1) = r; cols(p+1) = idxS(i,k,j); vals(p+1) = -1; p = p+1;
            bineq(r) = 0;
        end
    end
end

% (d) strong_lp constraints:
%   k*t(k,j) + sum_i s(i,k,j) <= gamma(k,j)*sum_i Q(i,j)
% => k*t + sum_i s - gamma*sum_i Q <= 0
for j = 1:nStages
    for k = 1:n1
        r = r + 1;

        % t term
        rows(p+1) = r; cols(p+1) = idxT(k,j); vals(p+1) = k; p = p+1;

        % sum_i s(i,k,j)
        for i = 1:n1
            rows(p+1) = r; cols(p+1) = idxS(i,k,j); vals(p+1) = 1; p = p+1;
        end

        % -gamma(k,j) * sum_i Q(i,j)
        gkj = gamma(k,j);
        base = (j-1)*n1;
        for i = 1:n1
            rows(p+1) = r; cols(p+1) = base + i; vals(p+1) = -gkj; p = p+1;
        end

        bineq(r) = 0;
    end
end

% (e) SOH-consistent ordering of final stored charge:
% Sort SOH ascending; enforce Q_final(sorted(r)) <= Q_final(sorted(r+1))
SOH_true = Q_max./Q_no;
[~, idxAsc] = sort(SOH_true, 'ascend');
for rr = 1:(n1-1)
    a = idxAsc(rr);
    b = idxAsc(rr+1);

    r = r + 1;
    % sum_j Q(a,j) - sum_j Q(b,j) <= Q_init(b) - Q_init(a)
    for j = 1:nStages
        rows(p+1) = r; cols(p+1) = idxQ(a,j); vals(p+1) = 1;  p = p+1;
        rows(p+1) = r; cols(p+1) = idxQ(b,j); vals(p+1) = -1; p = p+1;
    end
    bineq(r) = Q_init(b) - Q_init(a);
end

% (f) discharge utilization sufficient condition (k=1..n1-1),
% using ordering by SOH descending ("(i)" = i-th highest SOH).
[~, idxDesc] = sort(SOH_true, 'descend');
for k = 1:(n1-1)
    beta = sum(dprime(1:k)) / sumdprime;  % RHS ratio
    topSet = idxDesc(1:k);

    % Build per-cell coefficient: top -> (1-beta), others -> (-beta)
    coeffCell = -beta * ones(n1,1);
    coeffCell(topSet) = 1 - beta;

    r = r + 1;
    % Apply same per-cell coefficient across all stages
    for j = 1:nStages
        base = (j-1)*n1;
        idxs = base + (1:n1);
        rows(p+(1:n1)) = r;
        cols(p+(1:n1)) = idxs;
        vals(p+(1:n1)) = coeffCell;
        p = p + n1;
    end

    bineq(r) = beta*sum(Q_init) - sum(Q_init(topSet));
end

% Trim unused preallocation (in case of small n1)
rows = rows(1:p); cols = cols(1:p); vals = vals(1:p);

Aineq = sparse(rows, cols, vals, nIneq, N);

% -------------------- solve LP --------------------
try
    [x,fval,exitflag,output,lambda] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, params.linprogOptions);

catch ME
    % Save the failing case for debugging/support
    save("linprog_crash_case.mat","f","Aineq","bineq","Aeq","beq","lb","ub","ME");

    % Graceful fallback so your batch code can continue
    x = [];
    fval = NaN;
    exitflag = -999;  % custom "solver error" flag
    output = struct( ...
        "message", ME.message, ...
        "identifier", ME.identifier, ...
        "solver", "linprog", ...
        "note", "Caught exception; returned empty solution." );
    lambda = [];
end

if isempty(x)
    Q_opt = [];
    Q_final = Q_init;
else
    Q_opt = reshape(x(1:NQ), [n1, nStages]);
    Q_final = Q_init + sum(Q_opt, 2);
end
% -------------------- pack outputs --------------------
sol = struct();
sol.x = x;
sol.fval = fval;
sol.exitflag = exitflag;
sol.output = output;
sol.lambda = lambda;
sol.Q_final = Q_final;
sol.weights = w;

prob = struct();
prob.d = d;
prob.dprime = dprime;
prob.U_ter = U_ter;
prob.Aineq = Aineq; prob.bineq = bineq;
prob.Aeq = Aeq;     prob.beq = beq;
prob.lb = lb;       prob.ub = ub;
prob.f = f;
prob.index = struct('NQ',NQ,'Nt',Nt,'Ns',Ns,'N',N);
end

function [Q_opt, sol, prob] = solve_charge_allocation_equalSOC_lp( ...
    SOH, Q_init, Q_max, SOC_final, ...
    I_no_avg, SOC_max, ...
    U_phase, modulationType, fOCV, R0, Q_no, t_total, params)
%SOLVE_CHARGE_ALLOCATION_EQUALSOC_LP
% LP for the special case: all cells end at the SAME SOC after charging.
%
% Compared to the full LP, this version enforces:
%   (Equal-final-SOC)  (Q_init(i) + sum_j Q(i,j)) / Q_max(i) = SOC_final  for all i
% and checks ONLY:
%   (a) Charging-time limit
%   (b) Per-stage SOC upper bounds
% (Nonnegativity is always enforced.)
%
% ------------------------- Decision variables -------------------------
% Q(i,j): added charge (Ah) for cell i in stage j
%   i = 1..n1, j = 1..nStages   (your paper's j=0..m maps to MATLAB j=1..m+1)
%
% ------------------------- Inputs -------------------------
% SOH           [n1x1] cell SOH (0–1 or 0–100; auto-detect)
% Q_init        [n1x1] initial charge contents (Ah)
% Q_max         [n1x1] current max capacities (Ah)
% SOC_final     [1x1]  required common final SOC in [0,1]
%
% I_no_avg      [1x nStages] normalized avg C-rate per stage (C-rate), nStages=m+1
% SOC_max       [1x nStages] SOC upper bound at end of each stage k (0–1)
%
% U_phase       [1x1] phase-voltage magnitude (V)
% modulationType 'sinusoidal' or 'dc'
% fOCV          function handle U = fOCV(SOC), SOC in [0,1]
% R0            [1x1] avg source-side resistance (Ohm)
% Q_no          [1x1] nominal capacity (Ah) to map C-rate -> current
% t_total       [1x1] charging time limit (hours)
%
% params (optional struct) fields:
%   .SOH_EOL   (default 0.70)  EOL SOH threshold
%   .kappa     (default 0.1)
%   .epsSOH    (default 1e-3)
%   .Mbig      (default 1e6)
%   .linprogOptions (default dual-simplex, Display=none)
%
% ------------------------- Outputs -------------------------
% Q_opt  [n1 x nStages] optimal allocations (Ah)
% sol    struct:
%   .exitflag, .fval, .output, .lambda, .x
%   .Q_final   [n1x1] final charge contents (Ah)
%   .SOC_final [n1x1] final SOC (should be all SOC_final)
%   .time_used [1x1]  LHS of time constraint (hours)
%   .weights   [n1 x nStages] w(i,j)
% prob   struct:
%   .d        [n1 x nStages] duty cycles d(i,j)
%   .U_ter    [1 x nStages]  stage terminal-voltage proxies
%   .Aineq,.bineq,.Aeq,.beq,.lb,.ub,.f
%
% Requires: Optimization Toolbox (linprog)

% -------------------- defaults & checks --------------------
if nargin < 13 || isempty(params), params = struct(); end
if ~isfield(params,'SOH_EOL'),  params.SOH_EOL = 0.70; end
if ~isfield(params,'kappa'),    params.kappa   = 0.1; end
if ~isfield(params,'epsSOH'),   params.epsSOH  = 1e-3; end
if ~isfield(params,'Mbig'),     params.Mbig    = 1e6; end
if ~isfield(params,'linprogOptions')
    params.linprogOptions = optimoptions('linprog', ...
        'Algorithm','dual-simplex','Display','none');
end

SOH    = SOH(:);
Q_init = Q_init(:);
Q_max  = Q_max(:);

n1 = numel(SOH);
nStages = numel(I_no_avg);

assert(numel(Q_init)==n1 && numel(Q_max)==n1, 'SOH, Q_init, Q_max must have same length n1.');
assert(numel(SOC_max)==nStages, 'SOC_max must have same length as I_no_avg.');
assert(isa(fOCV,'function_handle'), 'fOCV must be a function handle.');
assert(SOC_final >= 0 && SOC_final <= 1, 'SOC_final must be in [0,1].');

% Auto-detect SOH scale
if max(SOH) > 1.5
    SOH = SOH/100;
end
SOH_EOL = params.SOH_EOL;
if SOH_EOL > 1.5
    SOH_EOL = SOH_EOL/100;
end

% -------------------- equal-final-SOC implies per-cell total added charge --------------------
DeltaQ = Q_max .* SOC_final - Q_init;         % required total added Ah per cell
if any(DeltaQ < -1e-12)
    error('Infeasible: some cells would require negative added charge to reach SOC_final.');
end
DeltaQ = max(DeltaQ, 0); % tiny numerical cleanup

% Also check last-stage SOC cap feasibility (necessary condition)
cap_end = Q_max .* SOC_max(end) - Q_init;     % max total addable by end
if any(DeltaQ > cap_end + 1e-12)
    error('Infeasible: some cells cannot reach SOC_final under SOC_max(end).');
end

% -------------------- compute weights w(i,j) as in your paper --------------------
kappa  = params.kappa;
epsSOH = params.epsSOH;
Mbig   = params.Mbig;

w = zeros(n1, nStages);
for j = 1:nStages
    factor = 1 + kappa * I_no_avg(j);
    for i = 1:n1
        if SOH(i) > SOH_EOL + epsSOH
            w(i,j) = factor / (SOH(i) - SOH_EOL)^2;
        else
            w(i,j) = factor * Mbig;
        end
    end
end

% -------------------- compute duty cycles d(i,j) and time coefficients --------------------
% Use the same proxy structure as your original model:
% SOC_avg = (SOC_init_avg + SOC_final)/2
SOC_init_avg = sum(Q_init) / sum(Q_max);
SOC_avg = min(max(0.5*(SOC_init_avg + SOC_final), 0), 1);

U_ocv = fOCV(SOC_avg);

U_ter = zeros(1,nStages);
for j = 1:nStages
    U_ter(j) = U_ocv + (Q_no * I_no_avg(j)) * R0;
end

d = zeros(n1,nStages);
switch lower(string(modulationType))
    case "sinusoidal"
        for j = 1:nStages
            for i = 1:n1
                ratio = ((2*i-1) * U_ter(j)) / (2*U_phase);
                ratio = min(ratio, 1);
                d(i,j) = (2/pi)*acos(ratio);
            end
        end
    case "dc"
        for j = 1:nStages
            for i = 1:n1
                if (i-1)*U_ter(j) > U_phase
                    d(i,j) = 0;
                elseif i*U_ter(j) < U_phase
                    d(i,j) = 1;
                else
                    d(i,j) = (U_phase - (i-1)*U_ter(j)) / U_ter(j);
                end
            end
        end
    otherwise
        error('modulationType must be ''sinusoidal'' or ''dc''.');
end

sumd = sum(d,1);
if any(sumd <= 0)
    error('Some stages have sum_i d(i,j)=0; time constraint becomes ill-defined.');
end
if any(I_no_avg <= 0)
    error('All I_no_avg entries must be > 0 to form the time constraint.');
end

% -------------------- build LP: variables x = Q(:) --------------------
% Q(i,j) is stored in x at index (j-1)*n1 + i
NQ = n1*nStages;

f = w(:);
lb = zeros(NQ,1);
ub = inf(NQ,1);

% Equality: for each cell i, sum_j Q(i,j) = DeltaQ(i)
Aeq = sparse(n1, NQ);
beq = DeltaQ;
for i = 1:n1
    cols = (0:(nStages-1))*n1 + i;
    Aeq(i, cols) = 1;
end

% Inequalities:
% (a) time limit:
%   sum_j [ (sum_i Q(i,j)) / (Q_no*I_no_avg(j)*sum_i d(i,j)) ] <= t_total
% (b) per-stage SOC upper bounds:
%   sum_{j=1..k} Q(i,j) <= Q_max(i)*SOC_max(k) - Q_init(i)
nTime = 1;
nSOC  = n1*nStages;
nIneq = nTime + nSOC;

% Preallocate triplets for sparse Aineq
nnz_est = NQ + sum(1:nStages)*n1;
rows = zeros(nnz_est,1);
cols = zeros(nnz_est,1);
vals = zeros(nnz_est,1);
bineq = zeros(nIneq,1);
p = 0; r = 0;

% (a) time row
r = r + 1;
for j = 1:nStages
    cj = 1 / (Q_no * I_no_avg(j) * sumd(j));  % hours per Ah in stage j
    base = (j-1)*n1;
    idxs = base + (1:n1);
    rows(p+(1:n1)) = r;
    cols(p+(1:n1)) = idxs;
    vals(p+(1:n1)) = cj;
    p = p + n1;
end
bineq(r) = t_total;

% (b) stage SOC caps
for i = 1:n1
    for k = 1:nStages
        r = r + 1;
        stageIdxs = (0:(k-1))*n1 + i;
        nn = numel(stageIdxs);
        rows(p+(1:nn)) = r;
        cols(p+(1:nn)) = stageIdxs;
        vals(p+(1:nn)) = 1;
        p = p + nn;

        bineq(r) = Q_max(i)*SOC_max(k) - Q_init(i);
    end
end

rows = rows(1:p); cols = cols(1:p); vals = vals(1:p);
Aineq = sparse(rows, cols, vals, nIneq, NQ);

% -------------------- solve --------------------
[x, fval, exitflag, output, lambda] = linprog( ...
    f, Aineq, bineq, Aeq, beq, lb, ub, params.linprogOptions);
if isempty(x)
    Q_opt = [];
    Q_final = Q_init;
else
    Q_opt = reshape(x, [n1, nStages]);
    Q_final = Q_init + sum(Q_opt, 2);
end
SOC_final_vec = Q_final ./ Q_max;

% compute time used
%cj = 1 ./ (Q_no .* I_no_avg(:)' .* sumd); % 1xnStages
%time_used = sum( sum(Q_opt,1) .* cj );

% -------------------- pack outputs --------------------
sol = struct();
sol.x = x;
sol.fval = fval;
sol.exitflag = exitflag;
sol.output = output;
sol.lambda = lambda;
sol.Q_final = Q_final;
sol.SOC_final = SOC_final_vec;
%sol.time_used = time_used;
sol.weights = w;

prob = struct();
prob.d = d;
prob.U_ter = U_ter;
prob.Aineq = Aineq; prob.bineq = bineq;
prob.Aeq = Aeq;     prob.beq = beq;
prob.lb = lb;       prob.ub = ub;
prob.f = f;
end