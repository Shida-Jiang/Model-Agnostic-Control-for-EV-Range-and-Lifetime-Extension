% long term aging experiment
% Very complex aging model
% purposely distribute aging on a few cells
close all;
clear
clc;
celltype=1;
%% Safety guarantees and reference-case configuration
% All SOH values are fractions (e.g., 0.70 means 70% SOH).
cfg.initialAgingOffsetDays = 100;
cfg.lmoDODCutoff = 0.001;

% Aging-model sensitivity switches.
cfg.calendarAgingScale = 1.0;  % Table I: set to 2 for calendar aging x2
cfg.cyclicAgingScale = 1.0;    % Table I: set to 2 for cyclic aging x2

% EOL, heterogeneity, and SOH-estimation settings.
cfg.randomSeed = 100;
cfg.packSOHEOL = 0.7;          % pack EOL based on usable pack capacity
cfg.cellSOHRetire = 0.6;       % individual-cell retirement threshold
cfg.kneeSOH = 0.75;            % onset of late-life acceleration
cfg.kneeAlpha = 25;            % multiplier = 1 + alpha*(SOH_knee - SOH)
cfg.gammaStd = 0.1;           % persistent cell-to-cell aging-rate variation
cfg.sohNoiseStd = 0.02;        % 2 percentage-point SOH estimation noise
cfg.sohBias = 0;            % additive SOH-estimation bias

% Controller objective. The SOC-balancing baseline keeps its own fixed,
% cell-independent stage weight so that these settings change only the
% proposed controller.
cfg.weightExponent = 2;        % exponent p in the SOH weighting law
cfg.kappa = 0.1;              % C-rate weight for the proposed controller
cfg.socBaselineKappa = 0.10;   % fixed stage weight for the SOC baseline
cfg.epsSOH = 1e-3;
cfg.Mbig = 1e6;                % absolute floor for the EOL penalty
cfg.MbigFactor = 10;           % M >= MbigFactor/epsSOH^p

% Proposed-controller LP formulation:
%   "full"    - original epigraph formulation; exact stage-wise
%               majorization representation, but O(n1^2) variables.
%   "ordered" - require each stage allocation to follow the slowly varying
%               capacity-based SOH order, then impose majorization directly.
%               This is more restrictive but uses only Q(i,j) variables.
cfg.lpFormulation = "full";   % set to "ordered" for the reduced formulation

% Controller-input mismatch.
% For the SOC-bias case, the controller shifts the initial and target SOC
% by the same event-wise amount. This preserves the requested session Ah
% and therefore the true pack SOC trajectory used for strategy comparison.
cfg.socBias = 0;            % Table I: set to 0.05 for +5 pp bias
cfg.socEstimateMax = 0.995;    % upper bound for the biased SOC estimate

% Capacity-error support is retained but disabled for the SOC-bias case.
% Capacity error is multiplicative, zero-mean across cells, and persistent.
cfg.capacityErrorStd = 0.00;
cfg.capacityErrorBias = 0.00;
cfg.capacitySafetySOC = 1;     % true-SOC saturation used by online correction
cfg.ocvBiasV = 0;          % Table I: set to -0.020 V

% Representative resistance used by the optimization. Exact cell-specific
% R0 values would make duty cycles depend on the cell-to-level assignment.
cfg.R0Physical = 0.01;         % true proxy used for diagnostics
cfg.R0Design = 0.01;           % resistance used by the optimization
% Table I scalar-resistance case: set BOTH values to 0.05 Ohm. Setting only
% R0Design to 0.05 instead tests a conservative design margin.

% Usage-profile sensitivity.
cfg.DCShare = 0.322;           % Table I: set to 0.50

n1=20;%number of cells

% Optional throughput-coupled temperature sensitivity. The reference case
% retains the original static per-cell temperature rise. For the thermal
% sensitivity case, set temperatureMode to "ah_proportional". The event-wise
% mean rise is then held at 10 deg C, while each cell's rise is proportional
% to its absolute Ah throughput during that charge or discharge event.
cfg.temperatureMode = "static";       % "static" or "ah_proportional"
cfg.meanOperatingDeltaT = 10.0;        % deg C for ah_proportional mode

% Runtime and robustness diagnostics.
cfg.collectRuntimeDiagnostics = true;
cfg.runtimeOutputPrefix = "runtime_summary";

%% Load Data
DATA = csvread("data_scott.csv",1,0);
original_time=3600*24*cfg.initialAgingOffsetDays;
stepnum=length(DATA);
stepcSOC=[100;DATA(:,6)]/100;
stepdSOC=DATA(:,3)/100;
chargetime=DATA(:,5);
chargestarttime=DATA(:,1)+original_time;
chargeendtime=DATA(:,2)+original_time;
dischargestarttime=DATA(:,2)+original_time;
DCcurrent=4*ones(stepnum,1);
dischargecurrent=2*ones(stepnum,1);
stepT=25*ones(stepnum,1);
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
% chargeendtime and dischargestarttime represent the same physical instant.
dischargestarttime=chargeendtime;
for i=1:stepnum
    if(stepcSOC(i+1)>0.995)
        stepcSOC(i+1)=0.995;
    end
end

% Common-mode SOC bias seen by the controller in each charging event.
% The target estimate is clipped at cfg.socEstimateMax, and the same
% applied bias is used for the initial estimate so that Delta SOC and the
% requested charging Ah remain unchanged.
socEstimateHeadroom=max(cfg.socEstimateMax-stepcSOC(2:end),0);
cfg.socBiasApplied=min(cfg.socBias*ones(stepnum,1),socEstimateHeadroom);
clear socEstimateHeadroom

%% Select DC fast-charging events
averageCrate=zeros(stepnum,1);
for i=1:stepnum
    averageCrate(i)=(stepcSOC(i+1)-stepdSOC(i))/chargetime(i)*3600;
end
[~,idx]=sort(averageCrate,"descend");
numDC=round(cfg.DCShare*stepnum);
numDC=min(max(numDC,0),stepnum);
whether_DC=zeros(stepnum,1);
whether_DC(idx(1:numDC))=1;

%% Optional comparison configuration
% Controls the optional second long-term comparison and the figure:
%   "none"     - skip Algorithm 2 and skip figure generation/export.
%   "exact"    - use the original proposed controller with exact SOH.
%   "capacity" - use ideal remaining-capacity balancing with exact capacity.
cfg.secondPanelMode = "capacity";
cfg.capacityBalanceMaxSOC = 1;
% Capacity balancing is used as an ideal comparison. By default, its
% analytical targets are not required to satisfy the LSPWM majorization
% condition. The charging-time constraint and SOC/C-rate limits remain.
cfg.capacityBalanceEnforceAchievability = false;
cfg.secondPanelMode = lower(string(cfg.secondPanelMode));
assert(any(cfg.secondPanelMode == ["none","exact","capacity"]), ...
    'cfg.secondPanelMode must be "none", "exact", or "capacity".');

assert(cfg.packSOHEOL >= cfg.cellSOHRetire, ...
    'Pack EOL must not be lower than the cell-retirement threshold.');
assert(cfg.calendarAgingScale>0 && cfg.cyclicAgingScale>0, ...
    'Aging scale factors must be positive.');
assert(cfg.DCShare>=0 && cfg.DCShare<=1, 'cfg.DCShare must be in [0,1].');
assert(cfg.socBias>=0, 'cfg.socBias must be nonnegative.');
assert(cfg.socEstimateMax>0 && cfg.socEstimateMax<=1, ...
    'cfg.socEstimateMax must lie in (0,1].');
assert(numel(cfg.socBiasApplied)==stepnum, ...
    'cfg.socBiasApplied must contain one value per charging event.');
assert(cfg.capacityErrorStd>=0, 'cfg.capacityErrorStd must be nonnegative.');
assert(cfg.epsSOH>0 && cfg.MbigFactor>1, ...
    'epsSOH must be positive and MbigFactor must exceed one.');
cfg.lpFormulation = lower(string(cfg.lpFormulation));
assert(isscalar(cfg.lpFormulation) && ...
    any(cfg.lpFormulation == ["full","ordered"]), ...
    'cfg.lpFormulation must be "full" or "ordered".');
assert(cfg.R0Physical>=0 && cfg.R0Design>=0, ...
    'Representative resistances must be nonnegative.');
assert(cfg.capacitySafetySOC>0 && cfg.capacitySafetySOC<=1, ...
    'cfg.capacitySafetySOC must lie in (0,1].');
cfg.temperatureMode = lower(string(cfg.temperatureMode));
assert(any(cfg.temperatureMode == ["static","ah_proportional"]), ...
    'cfg.temperatureMode must be "static" or "ah_proportional".');
assert(cfg.meanOperatingDeltaT>=0, ...
    'cfg.meanOperatingDeltaT must be nonnegative.');
if cfg.R0Design < cfg.R0Physical
    warning(['R0Design is below R0Physical. This is allowed for testing, ', ...
        'but it is not a conservative voltage-drop proxy.']);
end
if cfg.kneeAlpha > 0 && cfg.kneeSOH <= cfg.cellSOHRetire
    warning('The knee begins at or below the cell-retirement threshold.');
end

multiplier=24;
if(celltype==2)
    multiplier=8;
end
cellnum_in_group=1;
m=6;%SOC division number during CV stage
Uphase=50;

% Aliases retained to minimize changes to the original simulation code.
SOH_threshold1=cfg.kneeSOH;
Pack_SOH_end=cfg.packSOHEOL;
Cell_SOH_end=cfg.cellSOHRetire;

%% battery initial states generation
rng(cfg.randomSeed)%set random seed
options=optimoptions('linprog','Display','none', 'Algorithm', 'dual-simplex');
gamma=1+cfg.gammaStd*randn(n1,1);
initialCapacityStd=0.001;%0
Qmax=2.3;
Qmax1=Qmax.*(1+initialCapacityStd*randn(n1,1));
%Qmax22=Qmax.*(0.8+0.01*randn(n1,1));%linear capacity fade
%Qmax33=Qmax.*(0.7+0.01*randn(n1,1));%exponential capacity fade
% Original static per-cell operating-temperature rise used in the reference
% case. It is intentionally still sampled in ah_proportional mode so that the
% random-number sequence remains unchanged across sensitivity cases.
deltaT=10+2*randn(n1,1);
%deltaT=sort(deltaT,'descend');
%deltaT=zeros(n1,1);

Ahrecord1=zeros(1,stepnum);
Ahrecord2=zeros(1,stepnum);

fprintf('\nReference-case setup:\n');
fprintf('  Pack EOL / cell retirement = %.0f%% / %.0f%%\n', ...
    100*Pack_SOH_end, 100*Cell_SOH_end);
fprintf('  Knee SOH / alpha = %.0f%% / %.1f\n', ...
    100*SOH_threshold1, cfg.kneeAlpha);
fprintf('  Aging-rate variation sigma_gamma = %.3f\n', cfg.gammaStd);
fprintf('  SOH noise sigma = %.3f\n', cfg.sohNoiseStd);
fprintf('  Objective p / kappa = %.1f / %.2f\n', ...
    cfg.weightExponent, cfg.kappa);
fprintf('  Proposed-controller LP formulation = %s\n', cfg.lpFormulation);
fprintf('  Calendar / cyclic aging scales = %.2f / %.2f\n', ...
    cfg.calendarAgingScale, cfg.cyclicAgingScale);
fprintf('  SOC bias / capacity-error sigma / OCV bias = %.1f pp / %.3f / %.1f mV\n', ...
    100*cfg.socBias, cfg.capacityErrorStd, 1000*cfg.ocvBiasV);
if cfg.socBias>0
    clippedSOCBiasEvents=sum(cfg.socBiasApplied<cfg.socBias-1e-12);
    fprintf('  Applied SOC bias: median %.2f pp; %d of %d events clipped at %.1f%%.\n', ...
        100*median(cfg.socBiasApplied),clippedSOCBiasEvents,stepnum, ...
        100*cfg.socEstimateMax);
end
fprintf('  R0 physical / design = %.4f / %.4f Ohm\n', ...
    cfg.R0Physical, cfg.R0Design);
fprintf('  DC fast-charging events = %d of %d (%.1f%%)\n', ...
    numDC, stepnum, 100*numDC/stepnum);
if cfg.temperatureMode == "ah_proportional"
    fprintf('  Temperature mode / event-wise mean rise = %s / %.1f deg C\n', ...
        cfg.temperatureMode,cfg.meanOperatingDeltaT);
else
    fprintf('  Temperature mode / sampled mean rise = %s / %.2f deg C\n', ...
        cfg.temperatureMode,mean(deltaT));
end
fprintf('  Initial age offset / LMO DoD cutoff = %.1f days / %.3f\n\n', ...
    cfg.initialAgingOffsetDays, cfg.lmoDODCutoff);
%% algorithm 1: SOC-balancing baseline
group=n1;
n2=1;%retained for compatibility with the original code
Iremember=zeros(stepnum,1);
cfgRun=cfg;
cfgRun.strategyMode="soc";
[SOH_algorithm1, SOH_algorithm1_order, SOH_loss_cyc1, SOH_loss_cal1, Lifetime1, Iremember, iplot, runtime1]=simulation_main(celltype,n1,n2,m,group,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,DCcurrent,stepcSOC,stepdSOC,Qmax1,Qmax,SOH_threshold1,Pack_SOH_end,Cell_SOH_end,0,Iremember,0,gamma,multiplier,0,options,cfgRun);

%% optional algorithm 2 / second panel
runSecondPanel = cfg.secondPanelMode ~= "none";
SOH_algorithm2=[];
SOH_algorithm2_order=[];
SOH_loss_cyc2=[];
SOH_loss_cal2=[];
Lifetime2=NaN;
secondPanelName="";
runtime2=[];

if runSecondPanel
    group=cellnum_in_group;
    n2=n1/group;%remaining-capacity balancing is used during discharge
    cfgRun=cfg;

    switch cfg.secondPanelMode
        case "exact"
            cfgRun.strategyMode="proposed";
            secondPanelName="Ours w/o noise";
            secondPanelNoise=0;
        case "capacity"
            cfgRun.strategyMode="capacity";
            secondPanelName="Capacity balancing";
            secondPanelNoise=0; % ideal capacity/SOH information by design
    end

    [SOH_algorithm2, SOH_algorithm2_order, SOH_loss_cyc2, SOH_loss_cal2, Lifetime2, ~, ~, runtime2]=simulation_main(celltype,n1,n2,m,group,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,DCcurrent,stepcSOC,stepdSOC,Qmax1,Qmax,SOH_threshold1,Pack_SOH_end,Cell_SOH_end,1,Iremember,secondPanelNoise,gamma,multiplier,0,options,cfgRun);

    fprintf('%s lifetime change relative to SOC balancing: %.3f%%\n', ...
        secondPanelName,100*(Lifetime2-Lifetime1)/Lifetime1);
end

%% algorithm 3: proposed controller with the reference SOH noise
SOH_noise=cfg.sohNoiseStd;
group=cellnum_in_group;
n2=n1/group;%remaining-capacity balancing is used during discharge
cfgRun=cfg;
cfgRun.strategyMode="proposed";
[SOH_algorithm3, SOH_algorithm3_order, SOH_loss_cyc3, SOH_loss_cal3, Lifetime3, ~, ~, runtime3]=simulation_main(celltype,n1,n2,m,group,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,DCcurrent,stepcSOC,stepdSOC,Qmax1,Qmax,SOH_threshold1,Pack_SOH_end,Cell_SOH_end,1,Iremember,SOH_noise,gamma,multiplier,0,options,cfgRun);

fprintf('Proposed controller with %.1f-pp SOH noise: %.3f%% lifetime improvement.\n', ...
    100*SOH_noise,100*(Lifetime3-Lifetime1)/Lifetime1);


%% Runtime and robustness summary
if cfg.collectRuntimeDiagnostics
    proposedRuntimeName="Proposed controller ("+cfg.lpFormulation+")";
    runtimeNames=["SOC balancing"; proposedRuntimeName];
    runtimeList={runtime1; runtime3};
    if runSecondPanel
        runtimeNames=[runtimeNames(1); secondPanelName; runtimeNames(2)];
        runtimeList={runtime1; runtime2; runtime3};
    end
    runtimeSummary=build_runtime_summary(runtimeNames,runtimeList);
    disp(runtimeSummary)
    runtimeFile=append(cfg.runtimeOutputPrefix,"_celltype", ...
        int2str(celltype),".csv");
    writetable(runtimeSummary,runtimeFile);
    fprintf('Runtime summary written to %s.\n',runtimeFile);
    fprintf(['Online level assignment sorts %d cells per update: ', ...
        'O(n log n) time and O(n) memory.\n'],n1);
    normalizedMajorizationTol=1e-4;
    chargingTimeTolHours=1e-6;
    if runtime3.maxTrueStageMajorizationViolation>normalizedMajorizationTol || ...
            runtime3.maxTrueDischargeMajorizationViolation>normalizedMajorizationTol || ...
            runtime3.maxTrueChargingTimeOverrunHours>chargingTimeTolHours
        warning('simulation_main:NonnegligibleMismatchViolation', ...
            ['The proposed-controller schedule has a non-negligible ', ...
             'true-model achievability or charging-time violation. ', ...
             'Review the runtime CSV before reporting this sensitivity case.']);
    end
end

%% Figures
% In "none" mode, Algorithm 2 and all figure generation are skipped to save
% runtime. The SOC baseline and noisy proposed result are still computed.
if runSecondPanel
    % highest SOH, top 20% SOH, ..., worst SOH
    lifetimeForLimit=max([Lifetime1,Lifetime2,Lifetime3]);
    limit=ceil((multiplier*lifetimeForLimit*1.01)/365/3600/24/5)*5;

    % Legend labels with correct ordinal suffixes (2nd, 4th, ..., 20th).
    rankLabels=strings(10,1);
    for j=1:10
        k=n1/10*j;
        if mod(k,100)>=11 && mod(k,100)<=13
            suf="th";
        else
            switch mod(k,10)
                case 1, suf="st";
                case 2, suf="nd";
                case 3, suf="rd";
                otherwise, suf="th";
            end
        end
        rankLabels(j)=append(int2str(k),suf," highest SOH");
    end

    % Two-line titles: strategy name on line 1, lifetime change on line 2.
    % This avoids horizontal overflow in the 1x3 layout and reports the
    % metric as a lifetime change rather than "RUL".
    switch cfg.secondPanelMode
        case "exact"
            panel2Line1="(ii) Proposed controller (exact SOH)";
        case "capacity"
            panel2Line1="(ii) Ideal capacity balancing";
    end
    if SOH_noise>0
        panel3Line1=append("(iii) Proposed controller (\sigma_{SOH} = ", ...
            num2str(SOH_noise,'%.2f'),")");
    else
        panel3Line1="(iii) Proposed controller (exact SOH)";
    end
    panel2Line2=append("(",num2str(100*(Lifetime2-Lifetime1)/Lifetime1, ...
        '%+.1f'),"% lifetime)");
    panel3Line2=append("(",num2str(100*(Lifetime3-Lifetime1)/Lifetime1, ...
        '%+.1f'),"% lifetime)");

    f1=figure;
    tiledlayout(1,3,"TileSpacing","tight","Padding","compact");
    f1.Position = [100 100 980 520];

    x = (chargeendtime(1:stepnum)-chargeendtime(1))*multiplier;

    ax = nexttile;
    hold(ax,'on')
    nLines=15;
    C=parula(nLines);
    colororder(ax,C);
    ax.ColorOrderIndex=1;
    for j=1:10
        y=SOH_algorithm1_order(n1/10*j,:);
        plot(x(1:iplot)/365/24/3600,y(1:iplot)*100, ...
            DisplayName=rankLabels(j))
    end
    xlabel('Time (years)','FontSize',14)
    ylabel('SOH (%)','FontSize',14)
    set(gca,'Fontsize',14)
    legend(Location="southwest",FontSize=14)
    ylim([0 1]*100)
    xlim([0 limit])
    grid on
    title(["(i) Baseline method";"(SOC balancing)"], ...
        'FontSize',13.5,'FontWeight','normal');

    ax = nexttile;
    hold(ax,'on')
    colororder(ax,C);
    ax.ColorOrderIndex=1;
    for j=1:10
        plot([x;0]/365/24/3600,SOH_algorithm2_order(n1/10*j,:)*100, ...
            DisplayName=rankLabels(j))
    end
    xlabel('Time (years)','FontSize',14)
    set(gca,'Fontsize',14)
    legend(Location="southwest",FontSize=14)
    ylim([0 1]*100)
    xlim([0 limit])
    grid on
    title([panel2Line1;panel2Line2], ...
        'FontSize',13.5,'FontWeight','normal');

    ax = nexttile;
    hold(ax,'on')
    colororder(ax,C);
    ax.ColorOrderIndex=1;
    for j=1:10
        plot([x;0]/365/24/3600,SOH_algorithm3_order(n1/10*j,:)*100, ...
            DisplayName=rankLabels(j))
    end
    xlabel('Time (years)','FontSize',14)
    set(gca,'Fontsize',14)
    legend(Location="southwest",FontSize=14)
    ylim([0 1]*100)
    xlim([0 limit])
    grid on
    title([panel3Line1;panel3Line2], ...
        'FontSize',13.5,'FontWeight','normal');

    exportgraphics(f1,append('comparison',int2str(celltype),'.png'), ...
        'Resolution',600)
else
    fprintf('Second panel disabled: Algorithm 2 and figure export were skipped.\n');
end
%% functions
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

function factor = knee_multiplier(SOHnow, SOHknee, kneeAlpha)
%KNEE_MULTIPLIER Late-life acceleration used in the manuscript:
% factor = 1 + kneeAlpha*(SOHknee - SOHnow) below the knee, and 1 above it.
factor = 1;
if SOHnow < SOHknee
    factor = 1 + kneeAlpha*(SOHknee - SOHnow);
end
end

function Fade=lmo_capacity_loss_from_fd(fd)
%LMO_CAPACITY_LOSS_FROM_FD Xu et al. nonlinear SEI mapping.
alpha=0.0575;
beta=121;
Fade=1-alpha*exp(-beta*fd)-(1-alpha)*exp(-fd);
end

function ST=lmo_temperature_stress(T)
% Temperature T is in deg C, matching Xu et al.'s parameter table.
kT=0.0693;
Tref=25;
ST=exp(kT*(T-Tref)*(Tref+273.15)/(T+273.15));
end

function Ssig=lmo_soc_stress(SOC)
ks=1.04;
sigmaRef=0.5;
Ssig=exp(ks*(SOC-sigmaRef));
end

function fd0=lmo_initial_calendar_damage(ageSeconds,T,SOC)
% Latent linearized damage corresponding to the tunable initial-age offset.
% The associated capacity loss is intentionally not subtracted from the
% initialized SOH, preserving the original simulation convention.
kt=4.14e-10;
fd0=kt*max(ageSeconds,0)*lmo_soc_stress(SOC)*lmo_temperature_stress(T);
end

function [fdNew,FadePct]=agingmodelLMOcycle(T,SOC,DOD,fdOld,SOHnow,SOHknee,kneeAlpha,multiplier,dodCutoff)
% One call is presently treated as one rainflow HALF cycle (n=0.5).
% There is deliberately no second factor of 1/2 in agingcyc.
DOD=abs(DOD);
fdNew=fdOld;
FadePct=0;
if DOD<=dodCutoff
    return
end
if DOD>1+1e-9
    error('agingmodelLMOcycle:InvalidDOD', ...
        'LMO DoD must not exceed 1. Received %.6f.',DOD);
end
DOD=min(DOD,1);

kd1=1.4e5;
kd2=-0.501;
kd3=-1.23e5;
denominator=kd1*DOD^kd2+kd3;
if denominator<=0
    error('agingmodelLMOcycle:InvalidStress', ...
        'The Xu et al. DoD stress denominator is nonpositive.');
end
Sdelta=1/denominator;
Ssig=lmo_soc_stress(SOC);
ST=lmo_temperature_stress(T);
nCycle=0.5;
deltaFd=multiplier*nCycle*Sdelta*Ssig*ST;
fdNew=fdOld+deltaFd;

FadePct=100*(lmo_capacity_loss_from_fd(fdNew)- ...
    lmo_capacity_loss_from_fd(fdOld));
FadePct=FadePct*knee_multiplier(SOHnow,SOHknee,kneeAlpha);
end

function [fdNew,FadePct]=agingmodelLMOcalendar(deltat,T,SOC,fdOld,SOHnow,SOHknee,kneeAlpha,multiplier)
% Stateful calendar accumulation: past time is never reweighted using the
% current SOC or temperature.
fdNew=fdOld;
FadePct=0;
if deltat<=0
    return
end
kt=4.14e-10;
deltaFd=multiplier*kt*deltat* ...
    lmo_soc_stress(SOC)*lmo_temperature_stress(T);
fdNew=fdOld+deltaFd;

FadePct=100*(lmo_capacity_loss_from_fd(fdNew)- ...
    lmo_capacity_loss_from_fd(fdOld));
FadePct=FadePct*knee_multiplier(SOHnow,SOHknee,kneeAlpha);
end

function [agingStateNew,Fadecyc]=agingcyc(t,SOC,DOD,Qmax,deltaAh,T,crate,celltype,SOHnow,SOHknee,kneeAlpha,multiplier,agingState,lmoDODCutoff) %#ok<INUSD>
if(celltype==1)
    % Najera et al.: integrate absolute Ah throughput in both directions.
    Tem=T+273.15;
    a = 2.0916e-8;
    b = -1.2179e-5;
    c = 0.0018;
    d = -1.7082e-6;
    e = 0.0556;
    crate=abs(crate);
    k_crate=(a*Tem^2+b*Tem+c)*exp((d*Tem+e)*crate);
    Fadecyc=k_crate*abs(deltaAh)*multiplier/Qmax*100;
    Fadecyc=Fadecyc*knee_multiplier(SOHnow,SOHknee,kneeAlpha);
    agingStateNew=agingState;
else
    [agingStateNew,Fadecyc]=agingmodelLMOcycle( ...
        T,SOC,DOD,agingState,SOHnow,SOHknee,kneeAlpha, ...
        multiplier,lmoDODCutoff);
end
end

function [agingStateNew,Fadecal]=agingcal(t,deltat,T,SOC,celltype,Qmax,SOHnow,SOHknee,kneeAlpha,multiplier,agingState,initialAgingOffsetSeconds)
if(celltype==1)
    % Exact piecewise-constant integration of A(T,SOC)*t^(1/2).
    % The initial-age offset is not accelerated; only elapsed profile time is.
    agingStateNew=agingState;
    Fadecal=0;
    if deltat<=0
        return
    end
    Tem=T+273.15;
    f = 5.9808e6;
    g = 0.6898;
    h = -6.4647e3;
    A=f*exp(g*SOC)*exp(h/Tem);

    elapsedStart=max(t-initialAgingOffsetSeconds,0);
    effectiveT0=initialAgingOffsetSeconds+multiplier*elapsedStart;
    effectiveT1=effectiveT0+multiplier*deltat;
    day0=effectiveT0/86400;
    day1=effectiveT1/86400;

    Fadecal=A*(sqrt(day1)-sqrt(day0));
    Fadecal=Fadecal*100/Qmax;
    Fadecal=Fadecal*knee_multiplier(SOHnow,SOHknee,kneeAlpha);
else
    [agingStateNew,Fadecal]=agingmodelLMOcalendar( ...
        deltat,T,SOC,agingState,SOHnow,SOHknee,kneeAlpha,multiplier);
end
end


function [SOH_algorithm3, SOH_algorithm3_order, SOH_loss_cyc, SOH_loss_cal, Lifetime3, Irememberout, iplot, runtimeStats]=simulation_main(celltype,n1,n2,m,group,stepnum,stepT,deltaT,dischargecurrent,chargestarttime,chargeendtime,dischargestarttime,chargetime,Uphase,whether_DC,DCcurrent,stepcSOC,stepdSOC,Qmax1,Qmax,SOHknee,Pack_SOH_end,Cell_SOH_end,known_current,Iremember,SOH_noise,gamma,multiplier,AI,options,cfg)
rng(cfg.randomSeed)%set random seed

strategyMode=lower(string(cfg.strategyMode));
assert(any(strategyMode == ["soc","proposed","capacity"]), ...
    'cfg.strategyMode must be "soc", "proposed", or "capacity".');
runtimeStats=init_runtime_stats(strategyMode);

% The capacity-estimation error is persistent. A separate random stream keeps
% the existing SOH-noise sequence unchanged when capacityErrorStd is zero.
if strategyMode == "capacity"
    capacityEstimateFactor=ones(n1,1); % ideal benchmark uses exact capacity
else
    capacityEstimateFactor=make_capacity_estimate_factors(n1,cfg);
end

R0Physical=cfg.R0Physical;
R0Design=cfg.R0Design;
fOCVController=@(soc) OCV(soc,celltype)+cfg.ocvBiasV;

ttt=[];
agingState=zeros(n1,1);
initialAgingOffsetSeconds=cfg.initialAgingOffsetDays*86400;
if celltype==2 && initialAgingOffsetSeconds>0
    % The offset is evaluated at the initial storage condition: initial SOC
    % and the 25 deg C calendar-temperature assumption.
    initialFd=lmo_initial_calendar_damage( ...
        initialAgingOffsetSeconds,stepT(1),stepcSOC(1));
    agingState(:)=initialFd;
end
kneeAlpha=cfg.kneeAlpha;
Ahtotal3=100*ones(n1,1);

% Controller parameters. The surrogate SOH weight is referenced to the
% individual-cell retirement threshold, not to the pack-level EOL threshold.
params.SOH_EOL=Cell_SOH_end;
params.kappa=cfg.kappa;
params.weightExponent=cfg.weightExponent;
params.epsSOH=cfg.epsSOH;
params.Mbig=cfg.Mbig;
params.MbigFactor=cfg.MbigFactor;
params.socBaselineKappa=cfg.socBaselineKappa;
params.lpFormulation=cfg.lpFormulation;
params.linprogOptions=options;

Lifetime3=0;
Qnow=Qmax1;
Ahrecord3=zeros(1,stepnum);
memorysize=0;
if(AI==1)
    forgetcoe=0.9;
    SOH_noisy=[];
    Ah_memory=[];
end
SOH_loss_cyc=zeros(n1,stepnum);
SOH_loss_cal=zeros(n1,stepnum);
SOH_algorithm3=NaN(n1,stepnum+1);
SOH_algorithm3_order=NaN(n1,stepnum+1);
SOH_algorithm3(:,1)=Qmax1/Qmax;
SOH_algorithm3_order(:,1)=sort(SOH_algorithm3(:,1),'descend');
Irememberout=Iremember;
iplot=stepnum;
timesremember=1;
retirementWarningIssued=false;
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
    QcapTrueDis=Qmax.*SOH_algorithm3(:,i+1);
    capacityFactorDis=capacityEstimateFactor*sum(QcapTrueDis)/ ...
        sum(QcapTrueDis.*capacityEstimateFactor);
    QcapEstimatedDis=QcapTrueDis.*capacityFactorDis;
    if strategyMode == "soc"
        % Equalize the controller's estimated SOC. When capacityErrorStd=0,
        % this reduces exactly to the original SOC-balancing allocation.
        deltaAh=totalAh*QcapEstimatedDis/sum(QcapEstimatedDis);
    else
        % Remaining-capacity balancing uses the available capacity estimate.
        % Coulomb-counted discharge Ah is subtracted from this initial estimate.
        QremainingEstimated=Qnow.*capacityFactorDis;
        for ii=1:num
            [~,idx]=sort(QremainingEstimated-deltaAh,'descend');
            for jj=1:n1
                deltaAh(idx(jj))=deltaAh(idx(jj))+totalAh/num*duty2(jj)/sum(duty2);
            end
        end
    end
    deltaTDischarge=event_temperature_rise(deltaAh,deltaT,cfg);
    Ahtotal3=Ahtotal3+deltaAh;%*multiplier;
    Ahrecord3(i)=totalAh;
    % Aging during the discharge/idle interval preceding charge event i.
    % This interval is counted exactly once.
    if i==1
        calendarIntervalStart=initialAgingOffsetSeconds;
    else
        calendarIntervalStart=chargeendtime(i-1);
    end
    calendarIntervalDuration=max(chargestarttime(i)-calendarIntervalStart,0);

    for ii=1:n1
        DOD=deltaAh(ii)/SOH_algorithm3(ii,i+1)/Qmax;
        SOCtemp=(Qnow(ii)*2-deltaAh(ii))/SOH_algorithm3(ii,i+1)/Qmax/2;

        [agingState(ii),cyc_loss]=agingcyc( ...
            dischargestarttime(i),SOCtemp,DOD,Qmax,deltaAh(ii), ...
            stepT(i)+deltaTDischarge(ii),dischargecurrent(i)/Qmax,celltype, ...
            SOH_algorithm3(ii,i+1),SOHknee,kneeAlpha,multiplier, ...
            agingState(ii),cfg.lmoDODCutoff);
        cyc_loss=cfg.cyclicAgingScale*cyc_loss/100*gamma(ii);
        SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
        SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;

        if calendarIntervalDuration>0
            % Calendar temperature is assumed to be stepT (25 deg C here).
            [agingState(ii),cal_loss]=agingcal( ...
                calendarIntervalStart,calendarIntervalDuration,stepT(i), ...
                SOCtemp,celltype,Qmax,SOH_algorithm3(ii,i+1), ...
                SOHknee,kneeAlpha,multiplier,agingState(ii), ...
                initialAgingOffsetSeconds);
            cal_loss=cfg.calendarAgingScale*cal_loss/100*gamma(ii);
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
        deltaTDCCharge=event_temperature_rise(deltaAh,deltaT,cfg);
        % Instantaneous line-current C-rate while a cell is active.
        % Per-cell duty/throughput is already represented by deltaAh.
        DCrate=abs(DCcurrent(i))/Qmax;
        for ii=1:n1
           SOCtemp=(stepcSOC(i+1)+Qnow(ii)/SOH_algorithm3(ii,i+1)/Qmax)/2;
           DOD=stepcSOC(i+1)-Qnow(ii)/SOH_algorithm3(ii,i+1)/Qmax;
           [agingState(ii),cyc_loss]=agingcyc( ...
               chargestarttime(i),SOCtemp,DOD,Qmax,deltaAh(ii), ...
               stepT(i)+deltaTDCCharge(ii),DCrate,celltype, ...
               SOH_algorithm3(ii,i+1),SOHknee,kneeAlpha,multiplier, ...
               agingState(ii),cfg.lmoDODCutoff);
           cyc_loss=cfg.cyclicAgingScale*cyc_loss/100*gamma(ii);
           SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
           SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;

           % Calendar aging during the current charging interval only.
           [agingState(ii),cal_loss]=agingcal( ...
               chargestarttime(i),chargeendtime(i)-chargestarttime(i), ...
               stepT(i)+deltaTDCCharge(ii),SOCtemp,celltype,Qmax, ...
               SOH_algorithm3(ii,i+1),SOHknee,kneeAlpha,multiplier, ...
               agingState(ii),initialAgingOffsetSeconds);
           cal_loss=cfg.calendarAgingScale*cal_loss/100*gamma(ii);
           SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
           SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
        end
        Qnow=stepcSOC(i+1)*SOH_algorithm3(:,i+1)*Qmax;

    else
        %AC
        %calculate Qnow by optimization
        %record also the c-rates
        %calculate the duty cycles
        SOH_algorithm3(:,i+1)=max(0.1,SOH_algorithm3(:,i+1));
        acOptimizationTic=tic;
        lpCountBeforeSession=numel(runtimeStats.lpSolveTimes);

        % The controller sees the same positive SOC bias at the beginning
        % and target of this event. The physical simulation continues to use
        % the unshifted SOC trajectory from the recorded driving profile.
        if strategyMode == "capacity"
            socBiasApplied=0; % ideal benchmark uses exact SOC information
        else
            socBiasApplied=cfg.socBiasApplied(i);
        end
        SOCstartController=stepdSOC(i)+socBiasApplied;
        SOCtargetController=stepcSOC(i+1)+socBiasApplied;

        % Separate physical capacities from the values available to the
        % controller. The requested session Ah is kept fixed so that a
        % capacity-error case changes allocation, not the charging record.
        QcapTrueCharge=Qmax.*SOH_algorithm3(:,i+1);
        capacityFactorCharge=capacityEstimateFactor*sum(QcapTrueCharge)/ ...
            sum(QcapTrueCharge.*capacityEstimateFactor);
        QcapController=QcapTrueCharge.*capacityFactorCharge;
        % Separate the true stored charge from the controller estimate.
        % Adding socBiasApplied*QcapController shifts the estimated initial
        % SOC without changing the physical state. The target estimate is
        % shifted by the same amount because requestedAddedAh is still based
        % on the true recorded SOC rise.
        QinitTrueForLP=Qnow;
        QinitController=Qnow.*capacityFactorCharge + ...
            socBiasApplied*QcapController;
        requestedAddedAh=stepcSOC(i+1)*sum(QcapTrueCharge)-sum(QinitTrueForLP);
        requestedAddedAh=max(requestedAddedAh,0);
        QfinalSumController=sum(QinitController)+requestedAddedAh;
        SOCfinalController=QfinalSumController/sum(QcapController);

        % Capacity balancing is evaluated under ideal capacity/SOH
        % information. The SOC and proposed controllers retain their
        % original estimation treatment.
        if strategyMode == "capacity"
            SOH_now_noisy=SOH_algorithm3(:,i+1);
        else
            SOH_now_noisy=SOH_algorithm3(:,i+1)+cfg.sohBias+SOH_noise*randn(n1,1);
            if(AI==1)
                SOH_noisy=[SOH_noisy SOH_now_noisy];
                Ah_memory=[Ah_memory Ahtotal3];
                memorysize=memorysize+1;
            end
            deltaSOH=NaN(n1,1);
            if(AI==0)
                for iii=1:n1
                    deltaSOH(iii)=max(SOH_now_noisy(iii)-Pack_SOH_end,0.001);
                end
            else
                if(memorysize<3)
                    for iii=1:n1
                        deltaSOH(iii)=max(SOH_algorithm3(iii,i+1)+SOH_noise*randn()-Pack_SOH_end,0.001);
                    end
                else
                    w=ones(memorysize,1);
                    for iii=1:memorysize
                        w(iii)=forgetcoe^(memorysize-iii);
                    end
                    for iii=1:n1
                        y=SOH_noisy(iii,:).';
                        X=Ah_memory(iii,:).';
                        mdl=fitlm(X,y,'linear','Weights',w);
                        SOH_est=predict(mdl,Ahtotal3(iii));
                        deltaSOH(iii)=max(SOH_est-Pack_SOH_end,0.001);
                        SOH_now_noisy(iii)=SOH_est;
                        clear X y mdl SOH_est
                    end
                end
            end
        end
        if strategyMode == "capacity"
            % Ideal remaining-capacity balancing:
            %   1) determine final absolute charge targets by water filling,
            %   2) divide each cell's SOC rise across the CCCV stages, and
            %   3) reduce phase voltage only if needed for LSPWM achievability.
            % No LP is solved for this strategy.
            Tlimit_modified=chargetime(i)/3600;
            Q_final_sum_capacity=stepcSOC(i+1)*sum(SOH_algorithm3(:,i+1))*Qmax;

            if known_current && Iremember(i)>0
                IavgReference=Iremember(i);
            else
                IavgReference=NaN;
            end

            [Ahbest,Iavgbest,SOCgridbest,stageCrateForAging,Unow,capacityInfo] = ...
                build_capacity_balancing_schedule( ...
                    Qnow,Qmax.*SOH_algorithm3(:,i+1),Q_final_sum_capacity, ...
                    m,IavgReference,celltype,Qmax,R0Physical,Tlimit_modified, ...
                    cfg.capacityBalanceMaxSOC, ...
                    cfg.capacityBalanceEnforceAchievability);

            if known_current==0
                Irememberout(i)=Iavgbest;
            end
            if isfinite(IavgReference) && IavgReference>0
                timesremember=Iavgbest/IavgReference;
            else
                timesremember=1;
            end

            if capacityInfo.currentIncreased
                if cfg.capacityBalanceEnforceAchievability
                    reasonText = 'charging-time and achievability constraints';
                else
                    reasonText = 'charging-time constraint';
                end
                fprintf(['Capacity balancing increased the CC C-rate from ', ...
                    '%.4f to %.4f in event %d to satisfy the %s.\n'], ...
                    IavgReference,Iavgbest,i,reasonText);
            end
            % if ~cfg.capacityBalanceEnforceAchievability && ...
            %         capacityInfo.maxMajorizationViolation > 1e-6
            %     fprintf(['Ideal capacity balancing ignores an LSPWM ', ...
            %         'majorization violation of %.3e in event %d.\n'], ...
            %         capacityInfo.maxMajorizationViolation,i);
            % end
        else
        % Use maximum current, find the maximum voltage that makes charging
        % possible
        Tlimit_modified=chargetime(i)/3600;
        Unow = n1*fOCVController((SOCstartController+SOCtargetController)/2);
        Q_opt=[];
        Iavg=min(I_CV(0.1), I_CV(SOCstartController));
        while(known_current==1 && isempty(Q_opt) && Unow>4)
            SOC_CV=I_CV_inverse(Iavg);
            if(SOC_CV>=1)
                SOCgrid=ones(1,m+1);
            else
                SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
            end
            SOCgrid=max(SOCgrid,SOCstartController);
            SOCgrid2=[SOCstartController SOCgrid];
            SOCs=zeros(1,m+1);
            for ii=1:m+1
                SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
            end
            Is=I_CV(SOCs);
            Is(1)=Iavg;
            [Q_opt, solTry, probTry] = solve_charge_allocation_lp( ...
                SOH_now_noisy, QinitController, QcapController, ...
                QfinalSumController, Is, SOCgrid, ...
                Unow, 'sinusoidal', fOCVController, R0Design, Qmax, ...
                Tlimit_modified, Uphase, dischargecurrent(i)/Qmax, params);
            runtimeStats=record_lp_runtime(runtimeStats,solTry,probTry);
            Unow = Unow*0.95;
        end
        % Fix the voltage, find best current current

        %Find minimum current (C rate)
        if known_current && Iremember(i)>0
            Iavgmin=Iremember(i);
        else
            d=real(acos(fOCVController(0.5*SOCstartController+0.5*SOCtargetController)*(0:1:n1)/Unow))/pi*2;
            %optimization
            duty=zeros(1,n1);
            for ii=1:n1
                duty(ii)=0.5*(d(ii)+d(ii+1));
            end
            Iavgmin=sum(SOH_algorithm3(:,i+1))*(SOCtargetController-SOCstartController)/sum(duty)/Tlimit_modified;
        end
        % search for optimal current
        best=inf;
        Ahbest=[];
        stageCrateBest=[];
        probBest=[];
        while Iavg>Iavgmin*0.7
            SOC_CV=I_CV_inverse(Iavg);
            if(SOC_CV>=1)
                SOCgrid=ones(1,m+1);
            else
                SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
            end
            SOCgrid=max(SOCgrid,SOCstartController);
            SOCgrid2=[SOCstartController SOCgrid];
            SOCs=zeros(1,m+1);
            for ii=1:m+1
                SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
            end
            Is=I_CV(SOCs);
            Is(1)=Iavg;
            if strategyMode == "soc" % SOC balancing
                [Q_opt2, sol, prob] = solve_charge_allocation_equalSOC_lp( ...
                    SOH_now_noisy, QinitController, QcapController, ...
                    SOCfinalController, Is, SOCgrid, Unow, ...
                    'sinusoidal', fOCVController, R0Design, Qmax, ...
                    Tlimit_modified, params);
            else % SOH balancing
                [Q_opt2, sol, prob] = solve_charge_allocation_lp( ...
                    SOH_now_noisy, QinitController, QcapController, ...
                    QfinalSumController, Is, SOCgrid, Unow, ...
                    'sinusoidal', fOCVController, R0Design, Qmax, ...
                    Tlimit_modified, Uphase, dischargecurrent(i)/Qmax, params);
            end
            runtimeStats=record_lp_runtime(runtimeStats,sol,prob);
            if ~isempty(Q_opt2) && sol.exitflag>0 && ...
                    isfinite(sol.fval) && sol.fval<best
                best=sol.fval;
                Ahbest=Q_opt2;
                Iavgbest=Iavg;
                SOCgridbest=SOCgrid;
                stageCrateBest=Is;
                probBest=prob;
            end
            Iavg = Iavg*0.95;
        end
        if(isempty(Ahbest))
            disp(append('Failed: ',int2str(i)))
        else
            if(known_current==0)
                Irememberout(i)=Iavgbest;
            end
            timesremember=Iavgbest/Iavgmin;
        end
        if isempty(stageCrateBest)
            error('simulation_main:MissingStageCrate', ...
                'No AC charging-stage C-rate vector was stored for event %d.',i);
        end
        stageCrateForAging=stageCrateBest;
        end

        % Map the controller's biased stage boundaries back to the true
        % SOC scale before applying the existing online safety correction.
        SOCgridTrue=max(SOCgridbest-socBiasApplied,0);
        [Ahbest,safetyInfo]=enforce_true_stage_soc_limits( ...
            Ahbest,Qnow,QcapTrueCharge,SOCgridTrue,cfg.capacitySafetySOC);
        runtimeStats.capacitySafetyCorrections= ...
            runtimeStats.capacitySafetyCorrections+safetyInfo.corrected;
        runtimeStats.maxCommandedSOCViolation=max( ...
            runtimeStats.maxCommandedSOCViolation, ...
            safetyInfo.maxCommandedSOCViolation);

        % Use the executed total Ah of each cell to assign one effective
        % temperature rise over the AC charging event. This keeps the mean
        % thermal exposure fixed while penalizing cells that receive more Ah.
        deltaTACCharge=event_temperature_rise(sum(Ahbest,2),deltaT,cfg);

        % For diagnostics, evaluate the design OCV at the biased SOC
        % estimate while evaluating the physical OCV at the true SOC.
        fOCVControllerFromTrueSOC=@(socTrue) fOCVController( ...
            min(max(socTrue+socBiasApplied,0),cfg.socEstimateMax));
        trueDiag=evaluate_schedule_true_model( ...
            Ahbest,Qnow,QcapTrueCharge,stageCrateForAging,Unow, ...
            celltype,Qmax,R0Physical,Uphase,dischargecurrent(i)/Qmax, ...
            fOCVControllerFromTrueSOC,R0Design);
        runtimeStats.maxTrueStageMajorizationViolation=max( ...
            runtimeStats.maxTrueStageMajorizationViolation, ...
            trueDiag.maxStageMajorizationViolation);
        runtimeStats.maxTrueStageMajorizationDeficitAh=max( ...
            runtimeStats.maxTrueStageMajorizationDeficitAh, ...
            trueDiag.maxStageMajorizationDeficitAh);
        runtimeStats.maxTrueDischargeMajorizationViolation=max( ...
            runtimeStats.maxTrueDischargeMajorizationViolation, ...
            trueDiag.dischargeMajorizationViolation);
        runtimeStats.maxTrueChargingTimeOverrunHours=max( ...
            runtimeStats.maxTrueChargingTimeOverrunHours, ...
            max(trueDiag.timeUsed-Tlimit_modified,0));
        runtimeStats.maxVoltageProxyDifferenceV=max( ...
            runtimeStats.maxVoltageProxyDifferenceV, ...
            trueDiag.maxVoltageProxyDifferenceV);
        runtimeStats.lpSolvesPerACSession(end+1,1)= ...
            numel(runtimeStats.lpSolveTimes)-lpCountBeforeSession;
        runtimeStats.acOptimizationTimes(end+1,1)=toc(acOptimizationTic);

        SOC_copy=SOCnow;
        %AC charging cyclic aging
        SOH_memory=SOH_algorithm3(:,i+1);
        for ii=1:n1
            for iii=1:m+1
                %cyclic
                if(SOH_memory(ii)>0.11)
                    deltaAh=Ahbest(ii,iii);
                    SOCtemp=(2*SOCnow(ii)+deltaAh/Qmax/SOH_memory(ii))/2;
                    SOCnow(ii)=SOCnow(ii)+deltaAh/Qmax/SOH_memory(ii);
                    DOD=deltaAh/Qmax/SOH_memory(ii);
                    [agingState(ii),cyc_loss]=agingcyc( ...
                        chargestarttime(i),SOCtemp,DOD,Qmax,deltaAh, ...
                        stepT(i)+deltaTACCharge(ii),stageCrateForAging(iii), ...
                        celltype,SOH_algorithm3(ii,i+1),SOHknee,kneeAlpha, ...
                        multiplier,agingState(ii),cfg.lmoDODCutoff);
                    cyc_loss=cfg.cyclicAgingScale*cyc_loss/100*gamma(ii);
                else
                    cyc_loss=0;
                end
                SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cyc_loss;
                SOH_loss_cyc(ii,i)=SOH_loss_cyc(ii,i)+cyc_loss;
            end
        end
        %AC charging calendar aging
        for ii=1:n1
            deltat=chargeendtime(i)-chargestarttime(i);
            avgSOC=SOC_copy(ii)+sum(Ahbest(ii,:))/(SOH_memory(ii)*Qmax)/2;
            [agingState(ii),cal_loss]=agingcal( ...
                chargestarttime(i),deltat,stepT(i)+deltaTACCharge(ii),avgSOC, ...
                celltype,Qmax,SOH_algorithm3(ii,i+1),SOHknee,kneeAlpha, ...
                multiplier,agingState(ii),initialAgingOffsetSeconds);
            cal_loss=cfg.calendarAgingScale*cal_loss/100*gamma(ii);
            SOH_algorithm3(ii,i+1)=SOH_algorithm3(ii,i+1)-cal_loss;
            SOH_loss_cal(ii,i)=SOH_loss_cal(ii,i)+cal_loss;
        end
        Qnow = SOCnow.*SOH_memory.*Qmax;
        clear Ahbest SOC_copy SOH_memory
    end
    ttt=[ttt std(SOCnow) / mean(SOCnow)];

    SOH_algorithm3_order(:,i+1)=sort(SOH_algorithm3(:,i+1),'descend');
    [packSOH, retiredMask]=compute_pack_soh(SOH_algorithm3(:,i+1),Cell_SOH_end);

    % The original operating model is unchanged. This warning identifies the
    % uncommon case in which a cell retires before the pack reaches EOL, since
    % explicit post-retirement topology reconfiguration is not modeled here.
    if any(retiredMask) && packSOH>=Pack_SOH_end && ~retirementWarningIssued
        disp("A cell crossed the retirement threshold before pack EOL.");
        retirementWarningIssued=true;
    end

    disp(append(num2str(timesremember)," (",int2str(i)," cycles completed)"))
    if (packSOH<Pack_SOH_end && Lifetime3==0)
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
runtimeStats=finalize_runtime_stats(runtimeStats);
end


function [AhByStage,IccUsed,SOCstageUpper,stageCrate,UphaseUsed,info] = ...
    build_capacity_balancing_schedule( ...
    Q_init,Q_cap,Q_final_sum,m,IccReference,celltype,Q_no,R0,t_total,maxSOC, ...
    enforceAchievability)
%BUILD_CAPACITY_BALANCING_SCHEDULE
% Construct an ideal remaining-capacity-balancing charging schedule without
% solving an LP. Final absolute stored charge is equalized by water filling,
% subject to the current charge and a per-cell SOC ceiling. The reference
% CC C-rate from SOC balancing is retained whenever feasible. If necessary,
% it is increased only enough to meet the charging-time constraint. When
% enforceAchievability is true, the LSPWM majorization condition is also
% imposed by lowering the phase voltage. When false, this routine defines an
% ideal capacity-balancing benchmark and only reports the majorization error.

Q_init=Q_init(:);
Q_cap=Q_cap(:);
nCells=numel(Q_init);
tol=1e-9;

assert(numel(Q_cap)==nCells,'Q_init and Q_cap must have equal length.');
assert(all(Q_cap>0),'All active cells must have positive usable capacity.');
assert(maxSOC>0 && maxSOC<=1,'maxSOC must lie in (0,1].');
if nargin < 11 || isempty(enforceAchievability)
    enforceAchievability = true;
end
enforceAchievability = logical(enforceAchievability);

Q_target=equal_capacity_target(Q_init,Q_cap,Q_final_sum,maxSOC);
SOC_init=Q_init./Q_cap;
SOC_target=Q_target./Q_cap;

IccMax=min(I_CV(0.1),I_CV(min(SOC_init)));
IccMax=max(IccMax,1e-5);

if ~isfinite(IccReference) || IccReference<=0
    % Fallback only; in the intended comparison IccReference is inherited
    % from the SOC-balancing baseline.
    roughC=(Q_final_sum-sum(Q_init))/(Q_no*t_total*nCells);
    IccReference=min(IccMax,max(roughC,1e-3));
end

IccStart=min(max(IccReference,1e-5),IccMax);
IccCandidates=IccStart;
while IccCandidates(end)<IccMax*(1-1e-10)
    nextIcc=min(IccCandidates(end)*1.05,IccMax);
    if nextIcc<=IccCandidates(end)*(1+1e-12)
        break
    end
    IccCandidates(end+1)=nextIcc; %#ok<AGROW>
end

bestViolation=inf;
bestTime=inf;
bestUphase=NaN;

for ic=1:numel(IccCandidates)
    Icc=IccCandidates(ic);
    [AhCandidate,SOCupperCandidate,crateCandidate,SOCmidCandidate] = ...
        capacity_stage_allocation( ...
            SOC_init,SOC_target,Q_cap,Icc,m);

    SOCavg=(sum(Q_init)+sum(Q_target))/(2*sum(Q_cap));
    SOCavg=min(max(SOCavg,0),1);
    UphaseInitial=nCells*OCV(SOCavg,celltype);

    if enforceAchievability
        UphaseCandidates=UphaseInitial;
        while UphaseCandidates(end)>4
            nextU=0.95*UphaseCandidates(end);
            if nextU<=4
                break
            end
            UphaseCandidates(end+1)=nextU; %#ok<AGROW>
        end
    else
        % Ideal capacity-balancing benchmark: retain the maximum phase
        % voltage and do not reject the target based on majorization.
        UphaseCandidates=UphaseInitial;
    end

    for iu=1:numel(UphaseCandidates)
        UphaseCandidate=UphaseCandidates(iu);
        [timeUsed,maxViolation] = capacity_schedule_diagnostics( ...
            AhCandidate,SOCmidCandidate,crateCandidate,UphaseCandidate, ...
            celltype,Q_no,R0);

        if timeUsed<bestTime || ...
                (abs(timeUsed-bestTime)<1e-12 && maxViolation<bestViolation)
            bestViolation=maxViolation;
            bestTime=timeUsed;
            bestUphase=UphaseCandidate;
        end

        timeOK=timeUsed<=t_total*(1+1e-8);
        achievabilityOK=(~enforceAchievability) || maxViolation<=1e-8;

        if timeOK && achievabilityOK
            AhByStage=AhCandidate;
            IccUsed=Icc;
            SOCstageUpper=SOCupperCandidate;
            stageCrate=crateCandidate;
            UphaseUsed=UphaseCandidate;

            info=struct();
            info.Q_target=Q_target;
            info.SOC_target=SOC_target;
            info.timeUsed=timeUsed;
            info.maxMajorizationViolation=maxViolation;
            info.achievabilityEnforced=enforceAchievability;
            info.currentIncreased=IccUsed>IccReference*(1+1e-8);
            return
        end
    end
end

if enforceAchievability
    error('capacityBalancing:NoFeasibleSchedule', ...
        ['No analytical capacity-balancing schedule met both the charging ', ...
         'time and LSPWM-achievability requirements. Best maximum ', ...
         'majorization violation = %.3e, best time = %.3f h, phase voltage ', ...
         '= %.3f V.'],bestViolation,bestTime,bestUphase);
else
    error('capacityBalancing:NoFeasibleTimeSchedule', ...
        ['Ignoring LSPWM achievability was insufficient: no analytical ', ...
         'capacity-balancing schedule met the charging-time limit. Best ', ...
         'time = %.3f h, majorization violation = %.3e, phase voltage ', ...
         '= %.3f V.'],bestTime,bestViolation,bestUphase);
end
end


function Q_target=equal_capacity_target(Q_init,Q_cap,Q_final_sum,maxSOC)
%EQUAL_CAPACITY_TARGET Equalize final absolute stored charge by water filling.
% Cells may remain above the common level if charging cannot remove charge,
% and low-capacity cells may saturate at maxSOC.

Q_init=Q_init(:);
Q_cap=Q_cap(:);
Q_upper=maxSOC*Q_cap;
tol=1e-10;

if Q_final_sum<sum(Q_init)-tol
    error('capacityBalancing:TargetBelowInitial', ...
        'The requested final pack charge is below the initial charge.');
end
if Q_final_sum>sum(Q_upper)+tol
    error('capacityBalancing:TargetAboveLimit', ...
        'The requested final pack charge exceeds the SOC ceilings.');
end

qLow=min(Q_init);
qHigh=max(Q_upper);
for iter=1:100
    q=0.5*(qLow+qHigh);
    candidate=min(Q_upper,max(Q_init,q));
    if sum(candidate)<Q_final_sum
        qLow=q;
    else
        qHigh=q;
    end
end
Q_target=min(Q_upper,max(Q_init,0.5*(qLow+qHigh)));

% Remove the tiny bisection residual without changing the water-filling
% structure in any meaningful way.
residual=Q_final_sum-sum(Q_target);
if residual>tol
    idx=find(Q_target<Q_upper-tol,1,'first');
    if isempty(idx)
        error('capacityBalancing:Residual','No headroom remains for residual.');
    end
    Q_target(idx)=Q_target(idx)+residual;
elseif residual<-tol
    idx=find(Q_target>Q_init+tol,1,'first');
    if isempty(idx)
        error('capacityBalancing:Residual','No removable charge for residual.');
    end
    Q_target(idx)=Q_target(idx)+residual;
end

if max(abs(sum(Q_target)-Q_final_sum))>1e-7
    error('capacityBalancing:TargetMismatch', ...
        'The final capacity targets do not satisfy the pack target.');
end
end


function [AhByStage,SOCstageUpper,stageCrate,SOCstageMid] = ...
    capacity_stage_allocation(SOC_init,SOC_target,Q_cap,Icc,m)
%CAPACITY_STAGE_ALLOCATION Divide each cell's SOC rise across CCCV stages.

SOC_init=SOC_init(:);
SOC_target=SOC_target(:);
Q_cap=Q_cap(:);
nCells=numel(SOC_init);
nStages=m+1;
tol=1e-12;

SOCstart=min(SOC_init);
SOChigh=max(SOC_target);

if SOChigh<SOCstart-tol
    error('capacityBalancing:NegativeCharge', ...
        'At least one target would require net discharge during charging.');
end

SOCcc=I_CV_inverse(Icc);
SOCcc=min(max(SOCcc,SOCstart),SOChigh);

SOCstageUpper=zeros(1,nStages);
SOCstageUpper(1)=SOCcc;
if SOChigh>SOCcc+tol
    SOCstageUpper(2:end)=SOCcc+(1:m)/m*(SOChigh-SOCcc);
else
    SOCstageUpper(:)=SOChigh;
end

SOCstageLower=[SOCstart SOCstageUpper(1:end-1)];
SOCstageMid=0.5*(SOCstageLower+SOCstageUpper);

stageCrate=I_CV(SOCstageMid);
stageCrate(1)=Icc;
stageCrate=max(stageCrate,1e-6);

AhByStage=zeros(nCells,nStages);
for i=1:nCells
    for j=1:nStages
        deltaSOC=max(0, ...
            min(SOC_target(i),SOCstageUpper(j))- ...
            max(SOC_init(i),SOCstageLower(j)));
        AhByStage(i,j)=Q_cap(i)*deltaSOC;
    end
end

targetDelta=Q_cap.*(SOC_target-SOC_init);
if max(abs(sum(AhByStage,2)-targetDelta))>1e-8
    error('capacityBalancing:StageAllocationMismatch', ...
        'Stage allocations do not reproduce the final capacity targets.');
end
end


function [timeUsed,maxViolation] = capacity_schedule_diagnostics( ...
    AhByStage,SOCstageMid,stageCrate,U_phase,celltype,Q_no,R0)
%CAPACITY_SCHEDULE_DIAGNOSTICS Check charging time and Theorem-1
% majorization for the analytical stage allocations.

[nCells,nStages]=size(AhByStage);
timeUsed=0;
maxViolation=0;
tol=1e-12;

for j=1:nStages
    stageAh=sum(AhByStage(:,j));
    if stageAh<=tol
        continue
    end

    Uter=OCV(SOCstageMid(j),celltype)+Q_no*stageCrate(j)*R0;
    duty=zeros(nCells,1);
    for level=1:nCells
        ratio=((2*level-1)*Uter)/(2*U_phase);
        ratio=min(max(ratio,0),1);
        duty(level)=(2/pi)*acos(ratio);
    end

    sumDuty=sum(duty);
    if sumDuty<=tol
        timeUsed=inf;
        maxViolation=inf;
        return
    end

    timeUsed=timeUsed+stageAh/(Q_no*stageCrate(j)*sumDuty);

    lhs=cumsum(sort(AhByStage(:,j),'descend'))/stageAh;
    rhs=cumsum(sort(duty,'descend'))/sumDuty;
    maxViolation=max(maxViolation,max(lhs-rhs));
end
end


function deltaTEvent=event_temperature_rise(deltaAh,deltaTStatic,cfg)
%EVENT_TEMPERATURE_RISE Effective per-cell temperature rise for one event.
% In the reference "static" mode, this exactly returns the original sampled
% per-cell rises. In "ah_proportional" mode, the arithmetic mean rise is
% cfg.meanOperatingDeltaT and each cell's rise is proportional to its
% absolute Ah throughput during the event:
%   deltaT_i = meanDeltaT * |deltaAh_i| / mean(|deltaAh|).
% This is a controlled thermal-feedback sensitivity, not an electrothermal
% dynamics model.

deltaAh=abs(deltaAh(:));
deltaTStatic=deltaTStatic(:);
assert(numel(deltaAh)==numel(deltaTStatic), ...
    'deltaAh and deltaTStatic must have the same number of cells.');

switch lower(string(cfg.temperatureMode))
    case "static"
        deltaTEvent=deltaTStatic;
    case "ah_proportional"
        meanAh=mean(deltaAh);
        if meanAh<=1e-12
            deltaTEvent=cfg.meanOperatingDeltaT*ones(size(deltaAh));
        else
            deltaTEvent=cfg.meanOperatingDeltaT*deltaAh/meanAh;
        end
    otherwise
        error('Unknown cfg.temperatureMode: %s.',char(string(cfg.temperatureMode)));
end
end

function [packSOH, retiredMask] = compute_pack_soh(cellSOH, cellRetirementSOH)
%COMPUTE_PACK_SOH Usable pack SOH for the cell-level inverter topology.
% Retired cells contribute zero usable capacity while the nominal pack
% capacity in the denominator remains unchanged.
retiredMask = cellSOH <= cellRetirementSOH;
usableCellSOH = cellSOH;
usableCellSOH(retiredMask) = 0;
packSOH = mean(usableCellSOH);
end

function [Q_opt, sol, prob] = solve_charge_allocation_lp( ...
    SOH, Q_init, Q_max, Q_final_sum, I_no_avg, SOC_max, ...
    U_phase, modulationType, fOCV, R0, Q_no, t_total, ...
    U_dis_phase, I_no_dis_avg, params)
%SOLVE_CHARGE_ALLOCATION_LP  Solve the proposed-controller LP.
%
% Two formulations are available through params.lpFormulation:
%   "full"    - original epigraph representation of the stage-wise
%               majorization constraints. Decision variables are Q, t, and s.
%   "ordered" - reduced formulation. Each stage allocation follows the same
%               slowly varying capacity-based SOH order used by the hard final
%               charge and discharge constraints. The top-k sums are then
%               known, so majorization is imposed directly using Q only.
%
% In both cases, Q(i,j) is the stage-wise added charge (Ah) for cell i in
% stage j, with paper index j=0..m mapped to MATLAB j=1..(m+1).
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
%   .SOH_EOL  (default 0.70) cell-retirement SOH used in the surrogate weight
%   .kappa    (default 0.1)
%   .weightExponent (default 2)
%   .epsSOH   (default 1e-3)
%   .Mbig     (default 1e6, absolute floor)
%   .MbigFactor (default 10, enforces M >= MbigFactor/epsSOH^p)
%   .lpFormulation (default "full"; alternatives: "full", "ordered")
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
if nargin < 15 || isempty(params), params = struct(); end
if ~isfield(params,'SOH_EOL'),  params.SOH_EOL = 0.70; end
if ~isfield(params,'kappa'),    params.kappa   = 0.1; end
if ~isfield(params,'weightExponent'), params.weightExponent = 2; end
if ~isfield(params,'epsSOH'),   params.epsSOH  = 1e-3; end
if ~isfield(params,'Mbig'),     params.Mbig    = 1e6; end
if ~isfield(params,'MbigFactor'), params.MbigFactor = 10; end
if ~isfield(params,'lpFormulation'), params.lpFormulation = "full"; end
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
assert(params.weightExponent>0, 'weightExponent must be positive.');
assert(params.epsSOH>0 && params.MbigFactor>1, ...
    'epsSOH must be positive and MbigFactor must exceed one.');
lpFormulation=lower(string(params.lpFormulation));
assert(isscalar(lpFormulation) && ...
    any(lpFormulation == ["full","ordered"]), ...
    'params.lpFormulation must be "full" or "ordered".');

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
pWeight = params.weightExponent;
epsSOH = params.epsSOH;
% Keep the EOL branch strictly above the largest finite-branch weight near
% the cutoff for every tested exponent. This is essential for p=3.
Mbase = max(params.Mbig, params.MbigFactor/epsSOH^pWeight);

w = zeros(n1, nStages);
for j = 1:nStages
    factor = (1 + kappa*I_no_avg(j));
    for i = 1:n1
        if SOH(i) > SOH_EOL + epsSOH
            w(i,j) = factor / (SOH(i) - SOH_EOL)^pWeight;
        else
            w(i,j) = factor * Mbase;
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

% The hard ordering constraints use the slowly varying capacity-based SOH
% ranking. The instantaneous noisy SOH sample remains in the objective only,
% which avoids noise-induced hard rank reversals in small charging events.
SOH_order = Q_max./Q_no;
[~, idxAsc] = sortrows([SOH_order, (1:n1).'],[1 2]);
idxDesc = flipud(idxAsc);

if lpFormulation == "ordered"
    [Q_opt,sol,prob]=solve_charge_allocation_ordered_core( ...
        Q_init,Q_max,Q_final_sum,I_no_avg,SOC_max,Q_no,t_total, ...
        w,d,gamma,dprime,U_ter,idxAsc,idxDesc,params,Mbase);
    return
end

% -------------------- full epigraph formulation --------------------
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

% Objective. Scaling by a positive constant improves conditioning and
% does not change the optimizer. sol.fval is converted back to the original
% unscaled objective so candidate currents remain comparable.
objectiveScale=max(w(:));
fSolver = zeros(N,1);
fSolver(1:NQ) = w(:)/objectiveScale;

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

% (e) SOH-consistent ordering of final stored charge.
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
solveTic=tic;
try
    [x,fvalScaled,exitflag,output,lambda] = linprog( ...
        fSolver, Aineq, bineq, Aeq, beq, lb, ub, params.linprogOptions);
catch ME
    solveTime=toc(solveTic);
    % Save the failing case for debugging/support.
    save("linprog_crash_case.mat","fSolver","Aineq","bineq", ...
        "Aeq","beq","lb","ub","ME");

    % Graceful fallback so a batch run can continue.
    x = [];
    fvalScaled = NaN;
    exitflag = -999;  % custom "solver error" flag
    output = struct( ...
        "message", ME.message, ...
        "identifier", ME.identifier, ...
        "solver", "linprog", ...
        "note", "Caught exception; returned empty solution." );
    lambda = [];
end
if ~exist('solveTime','var')
    solveTime=toc(solveTic);
end

if isempty(x)
    Q_opt = [];
    Q_final = Q_init;
    fval = NaN;
else
    Q_opt = reshape(x(1:NQ), [n1, nStages]);
    Q_final = Q_init + sum(Q_opt, 2);
    fval = dot(w(:),x(1:NQ));
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
sol.objectiveScale = objectiveScale;
sol.fvalScaled = fvalScaled;
sol.solveTime = solveTime;
sol.Mbase = Mbase;
sol.formulation = "full";

prob = struct();
prob.formulation = "full";
prob.d = d;
prob.dprime = dprime;
prob.U_ter = U_ter;
prob.Aineq = Aineq; prob.bineq = bineq;
prob.Aeq = Aeq;     prob.beq = beq;
prob.lb = lb;       prob.ub = ub;
prob.f = fSolver;
prob.fUnscaled = [w(:); zeros(N-NQ,1)];
prob.index = struct('NQ',NQ,'Nt',Nt,'Ns',Ns,'N',N);
prob.nIneq = nIneq;
prob.nEq = size(Aeq,1);
prob.idxAscSOH = idxAsc;
prob.idxDescSOH = idxDesc;
modelInfo=whos('Aineq','bineq','Aeq','beq','lb','ub','fSolver');
prob.modelBytes=sum([modelInfo.bytes]);
end

function [Q_opt, sol, prob] = solve_charge_allocation_ordered_core( ...
    Q_init,Q_max,Q_final_sum,I_no_avg,SOC_max,Q_no,t_total, ...
    w,d,gamma,dprime,U_ter,idxAsc,idxDesc,params,Mbase)
%SOLVE_CHARGE_ALLOCATION_ORDERED_CORE Reduced proposed-controller LP.
%
% The only decision variables are Q(i,j). For each stage j, the added charge
% is constrained to follow the same capacity-based SOH order:
%   Q(highest-SOH,j) >= ... >= Q(lowest-SOH,j).
% With this fixed order, the k largest allocations are known in advance and
% the stage-wise majorization inequalities can be imposed directly. This
% removes all t(k,j) and s(i,k,j) epigraph variables. The resulting feasible
% set is more restrictive than the full formulation, but the SOH-balancing
% direction and all other controller constraints are retained.

Q_init=Q_init(:);
Q_max=Q_max(:);
I_no_avg=I_no_avg(:).';
SOC_max=SOC_max(:).';
dprime=dprime(:);

n1=numel(Q_init);
nStages=numel(I_no_avg);
NQ=n1*nStages;
N=NQ;
idxQ=@(i,j) (j-1)*n1+i;

sumd=sum(d,1);
sumdprime=sum(dprime);
Q_init_sum=sum(Q_init);

% Objective and nonnegative Q bounds.
objectiveScale=max(w(:));
fSolver=w(:)/objectiveScale;
lb=zeros(NQ,1);
ub=inf(NQ,1);

% Total added charge.
deltaQ_total=Q_final_sum-Q_init_sum;
Aeq=sparse(1,NQ);
Aeq(1,1:NQ)=1;
beq=deltaQ_total;

% Inequality blocks:
%   (a) charging time                                      1
%   (b) cumulative per-stage SOC limits                   n1*nStages
%   (c) per-stage SOH ordering                            (n1-1)*nStages
%   (d) direct stage-wise majorization                    (n1-1)*nStages
%   (e) final stored-charge ordering                      n1-1
%   (f) subsequent-discharge utilization                  n1-1
nAdj=max(n1-1,0);
nTime=1;
nSOC=n1*nStages;
nStageOrder=nAdj*nStages;
nStrong=nAdj*nStages;
nFinalOrder=nAdj;
nDis=nAdj;
nIneq=nTime+nSOC+nStageOrder+nStrong+nFinalOrder+nDis;

% Sparse triplet preallocation.
nnz_time=NQ;
nnz_soc=n1*(nStages*(nStages+1)/2);
nnz_stage_order=2*nStageOrder;
nnz_strong=n1*nStrong;
nnz_final_order=2*nStages*nFinalOrder;
nnz_dis=NQ*nDis;
nnz_est=nnz_time+nnz_soc+nnz_stage_order+nnz_strong+ ...
    nnz_final_order+nnz_dis;

rows=zeros(nnz_est,1);
cols=zeros(nnz_est,1);
vals=zeros(nnz_est,1);
bineq=zeros(nIneq,1);
p=0;
r=0;

% (a) Charging-time limit.
r=r+1;
for j=1:nStages
    cj=1/(Q_no*I_no_avg(j)*sumd(j));
    idxs=(j-1)*n1+(1:n1);
    rows(p+(1:n1))=r;
    cols(p+(1:n1))=idxs;
    vals(p+(1:n1))=cj;
    p=p+n1;
end
bineq(r)=t_total;

% (b) Per-stage SOC upper bounds.
for i=1:n1
    for k=1:nStages
        r=r+1;
        stageIdxs=(0:(k-1))*n1+i;
        nn=numel(stageIdxs);
        rows(p+(1:nn))=r;
        cols(p+(1:nn))=stageIdxs;
        vals(p+(1:nn))=1;
        p=p+nn;
        bineq(r)=Q_max(i)*SOC_max(k)-Q_init(i);
    end
end

% (c) Per-stage capacity allocation follows the capacity-based SOH order.
% idxAsc is lowest SOH to highest SOH, with deterministic index tie-breaking.
for j=1:nStages
    for rr=1:nAdj
        low=idxAsc(rr);
        high=idxAsc(rr+1);
        r=r+1;
        % Q(low,j) <= Q(high,j)
        rows(p+1)=r; cols(p+1)=idxQ(low,j);  vals(p+1)=1;  p=p+1;
        rows(p+1)=r; cols(p+1)=idxQ(high,j); vals(p+1)=-1; p=p+1;
        bineq(r)=0;
    end
end

% (d) Direct stage-wise majorization. Because of (c), idxDesc(1:k)
% identifies the k largest Q(i,j) values in stage j.
for j=1:nStages
    for k=1:nAdj
        gkj=gamma(k,j);
        topSet=idxDesc(1:k);
        coeffCell=-gkj*ones(n1,1);
        coeffCell(topSet)=1-gkj;

        r=r+1;
        idxs=(j-1)*n1+(1:n1);
        rows(p+(1:n1))=r;
        cols(p+(1:n1))=idxs;
        vals(p+(1:n1))=coeffCell;
        p=p+n1;
        bineq(r)=0;
    end
end

% (e) SOH-consistent ordering of final stored charge.
for rr=1:nAdj
    low=idxAsc(rr);
    high=idxAsc(rr+1);
    r=r+1;
    for j=1:nStages
        rows(p+1)=r; cols(p+1)=idxQ(low,j);  vals(p+1)=1;  p=p+1;
        rows(p+1)=r; cols(p+1)=idxQ(high,j); vals(p+1)=-1; p=p+1;
    end
    bineq(r)=Q_init(high)-Q_init(low);
end

% (f) Sufficient condition for subsequent discharge utilization.
for k=1:nAdj
    beta=sum(dprime(1:k))/sumdprime;
    topSet=idxDesc(1:k);
    coeffCell=-beta*ones(n1,1);
    coeffCell(topSet)=1-beta;

    r=r+1;
    for j=1:nStages
        idxs=(j-1)*n1+(1:n1);
        rows(p+(1:n1))=r;
        cols(p+(1:n1))=idxs;
        vals(p+(1:n1))=coeffCell;
        p=p+n1;
    end
    bineq(r)=beta*sum(Q_init)-sum(Q_init(topSet));
end

if r~=nIneq
    error('orderedLP:RowCountMismatch', ...
        'Built %d inequality rows, but expected %d.',r,nIneq);
end

rows=rows(1:p);
cols=cols(1:p);
vals=vals(1:p);
Aineq=sparse(rows,cols,vals,nIneq,NQ);

% Solve.
solveTic=tic;
try
    [x,fvalScaled,exitflag,output,lambda]=linprog( ...
        fSolver,Aineq,bineq,Aeq,beq,lb,ub,params.linprogOptions);
catch ME
    solveTime=toc(solveTic);
    save("linprog_crash_case_ordered.mat","fSolver","Aineq","bineq", ...
        "Aeq","beq","lb","ub","ME");
    x=[];
    fvalScaled=NaN;
    exitflag=-999;
    output=struct( ...
        "message",ME.message, ...
        "identifier",ME.identifier, ...
        "solver","linprog", ...
        "note","Caught exception in ordered formulation; returned empty solution.");
    lambda=[];
end
if ~exist('solveTime','var')
    solveTime=toc(solveTic);
end

if isempty(x)
    Q_opt=[];
    Q_final=Q_init;
    fval=NaN;
else
    Q_opt=reshape(x,[n1,nStages]);
    Q_final=Q_init+sum(Q_opt,2);
    fval=dot(w(:),x);
end

sol=struct();
sol.x=x;
sol.fval=fval;
sol.exitflag=exitflag;
sol.output=output;
sol.lambda=lambda;
sol.Q_final=Q_final;
sol.weights=w;
sol.objectiveScale=objectiveScale;
sol.fvalScaled=fvalScaled;
sol.solveTime=solveTime;
sol.Mbase=Mbase;
sol.formulation="ordered";

prob=struct();
prob.formulation="ordered";
prob.d=d;
prob.dprime=dprime;
prob.U_ter=U_ter;
prob.Aineq=Aineq; prob.bineq=bineq;
prob.Aeq=Aeq;     prob.beq=beq;
prob.lb=lb;       prob.ub=ub;
prob.f=fSolver;
prob.fUnscaled=w(:);
prob.index=struct('NQ',NQ,'Nt',0,'Ns',0,'N',N);
prob.nIneq=nIneq;
prob.nEq=size(Aeq,1);
prob.idxAscSOH=idxAsc;
prob.idxDescSOH=idxDesc;
modelInfo=whos('Aineq','bineq','Aeq','beq','lb','ub','fSolver');
prob.modelBytes=sum([modelInfo.bytes]);
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
% SOH           [n1x1] retained for interface compatibility; not used in the objective
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
%   .socBaselineKappa (default 0.1, fixed across proposed-objective cases)
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
if ~isfield(params,'socBaselineKappa'), params.socBaselineKappa = 0.1; end
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

% -------------------- cell-independent baseline weights --------------------
% Equal final estimated SOC fixes each cell's total added Ah. These weights
% only choose the charging stages. They are independent of SOH, p, and the
% proposed controller's kappa so the baseline remains unchanged in those
% sensitivity rows.
baselineKappa=params.socBaselineKappa;
w=repmat(1+baselineKappa*I_no_avg(:).',n1,1);

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

objectiveScale=max(w(:));
fSolver=w(:)/objectiveScale;
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
solveTic=tic;
try
    [x, fvalScaled, exitflag, output, lambda] = linprog( ...
        fSolver, Aineq, bineq, Aeq, beq, lb, ub, params.linprogOptions);
catch ME
    solveTime=toc(solveTic);
    x=[];
    fvalScaled=NaN;
    exitflag=-999;
    output=struct("message",ME.message,"identifier",ME.identifier, ...
        "solver","linprog","note","Caught exception; returned empty solution.");
    lambda=[];
end
if ~exist('solveTime','var')
    solveTime=toc(solveTic);
end
if isempty(x)
    Q_opt = [];
    Q_final = Q_init;
    fval=NaN;
else
    Q_opt = reshape(x, [n1, nStages]);
    Q_final = Q_init + sum(Q_opt, 2);
    fval=dot(w(:),x);
end
SOC_final_vec = Q_final ./ Q_max;

% -------------------- pack outputs --------------------
sol = struct();
sol.x = x;
sol.fval = fval;
sol.exitflag = exitflag;
sol.output = output;
sol.lambda = lambda;
sol.Q_final = Q_final;
sol.SOC_final = SOC_final_vec;
sol.weights = w;
sol.objectiveScale = objectiveScale;
sol.fvalScaled = fvalScaled;
sol.solveTime = solveTime;

prob = struct();
prob.d = d;
prob.U_ter = U_ter;
prob.Aineq = Aineq; prob.bineq = bineq;
prob.Aeq = Aeq;     prob.beq = beq;
prob.lb = lb;       prob.ub = ub;
prob.f = fSolver;
prob.fUnscaled = w(:);
prob.index = struct('NQ',NQ,'N',NQ);
prob.nIneq = nIneq;
prob.nEq = size(Aeq,1);
modelInfo=whos('Aineq','bineq','Aeq','beq','lb','ub','fSolver');
prob.modelBytes=sum([modelInfo.bytes]);
end

function factor=make_capacity_estimate_factors(nCells,cfg)
%MAKE_CAPACITY_ESTIMATE_FACTORS Persistent multiplicative capacity error.
% The sample is centered and normalized across the simulated pack so that
% cfg.capacityErrorStd controls cell-to-cell mismatch without changing the
% requested pack-level Ah. cfg.capacityErrorBias can add a common bias.

if cfg.capacityErrorStd==0
    relativeError=cfg.capacityErrorBias*ones(nCells,1);
else
    errorStream=RandStream('mt19937ar','Seed',cfg.randomSeed+104729);
    z=randn(errorStream,nCells,1);
    z=z-mean(z);
    zStd=std(z,0);
    if zStd<=eps
        z=zeros(nCells,1);
    else
        z=z/zStd;
    end
    relativeError=cfg.capacityErrorBias+cfg.capacityErrorStd*z;
end
factor=1+relativeError;
if any(~isfinite(factor)) || any(factor<=0)
    error('capacityEstimate:InvalidFactor', ...
        'Capacity-estimate factors must be finite and positive.');
end
end

function [AhSafe,info]=enforce_true_stage_soc_limits( ...
    AhCommand,QinitTrue,QcapTrue,SOCstageUpper,hardMaxSOC)
%ENFORCE_TRUE_STAGE_SOC_LIMITS Saturate and redistribute using true limits.
% This represents the fast online correction layer. It preserves the total
% requested Ah whenever the true pack has sufficient stage-wise headroom.
%
% A tiny SOC excess can arise from LP feasibility tolerances and floating-
% point roundoff. Correcting such an excess can move a negligible residual
% Ah into an otherwise empty stage and make a normalized majorization check
% look artificially large. Therefore, commands whose maximum excess is no
% larger than socNumericalTol are left unchanged.

AhCommand=max(AhCommand,0);
QinitTrue=QinitTrue(:);
QcapTrue=QcapTrue(:);
SOCstageUpper=SOCstageUpper(:).';
[nCells,nStages]=size(AhCommand);
assert(numel(QinitTrue)==nCells && numel(QcapTrue)==nCells, ...
    'AhCommand, QinitTrue, and QcapTrue must have compatible sizes.');
assert(numel(SOCstageUpper)==nStages, ...
    'SOCstageUpper must have one entry per charging stage.');

tol=1e-10;
socNumericalTol=1e-7; % fraction of SOC (0.00001 percentage points)
Qcommanded=QinitTrue;
maxViolation=0;
for j=1:nStages
    Qcommanded=Qcommanded+AhCommand(:,j);
    stageSOC=min(SOCstageUpper(j),hardMaxSOC);
    maxViolation=max(maxViolation, ...
        max(max(Qcommanded./QcapTrue-stageSOC,0)));
end

% Do not redistribute a schedule merely because of numerical roundoff.
if maxViolation<=socNumericalTol
    AhSafe=AhCommand;
    Qstate=QinitTrue+sum(AhSafe,2);
    info=struct();
    info.corrected=0;
    info.maxCommandedSOCViolation=maxViolation;
    info.totalRedistributedAh=0;
    info.maxFinalSOC=max(Qstate./QcapTrue);
    info.socNumericalTolerance=socNumericalTol;
    return
end

AhSafe=zeros(size(AhCommand));
Qstate=QinitTrue;
carry=0;
corrected=false;
for j=1:nStages
    desired=AhCommand(:,j);
    stageTotal=sum(desired)+carry;
    stageSOC=min(SOCstageUpper(j),hardMaxSOC);
    headroom=max(stageSOC*QcapTrue-Qstate,0);

    allocation=min(desired,headroom);
    residual=stageTotal-sum(allocation);
    available=max(headroom-allocation,0);
    if residual>tol && sum(available)>tol
        add=available*min(residual/sum(available),1);
        allocation=allocation+add;
        residual=residual-sum(add);
    end

    AhSafe(:,j)=allocation;
    Qstate=Qstate+allocation;
    carry=max(residual,0);
    if carry>tol || max(abs(allocation-desired))>tol
        corrected=true;
    end
end

if carry>1e-8
    error('capacityEstimate:TrueSOCLimitInfeasible', ...
        ['The estimated schedule exceeds the true pack headroom by %.6g Ah ', ...
         'after the final stage.'],carry);
end
if abs(sum(AhSafe(:))-sum(AhCommand(:)))>1e-7
    error('capacityEstimate:AhMismatch', ...
        'The true-SOC correction did not preserve total requested Ah.');
end

info=struct();
info.corrected=double(corrected);
info.maxCommandedSOCViolation=maxViolation;
info.totalRedistributedAh=0.5*sum(reshape(abs(AhSafe-AhCommand),[],1));
info.maxFinalSOC=max(Qstate./QcapTrue);
info.socNumericalTolerance=socNumericalTol;
end

function diag=evaluate_schedule_true_model( ...
    AhByStage,QinitTrue,QcapTrue,stageCrate,UphaseCharge,celltype,Qno, ...
    R0True,UdisPhase,IdisNo,fOCVDesign,R0Design)
%EVALUATE_SCHEDULE_TRUE_MODEL Diagnose OCV/R0 mismatch after optimization.
%
% The stage-wise majorization condition is homogeneous in stage Ah. A stage
% containing only a numerical residual can therefore show a large normalized
% violation even when the absolute infeasible Ah is negligible. We report the
% normalized violation only for active stages and separately retain the
% maximum absolute majorization deficit in Ah over all stages.

QinitTrue=QinitTrue(:);
QcapTrue=QcapTrue(:);
stageCrate=stageCrate(:).';
[nCells,nStages]=size(AhByStage);
assert(numel(stageCrate)==nStages,'stageCrate length mismatch.');

tol=1e-12;
sessionAh=sum(AhByStage(:));
% Ignore only stages smaller than one part per million of the session Ah,
% with an absolute floor of 1e-8 Ah. Such stages have no material effect on
% charging time or cell aging but can be dominated by solver roundoff.
activeStageAhTol=max(1e-8,1e-6*max(sessionAh,0));

QfinalTrue=QinitTrue+sum(AhByStage,2);
SOCavgTrue=(sum(QinitTrue)+sum(QfinalTrue))/(2*sum(QcapTrue));
SOCavgTrue=min(max(SOCavgTrue,0),1);
UocvTrue=OCV(SOCavgTrue,celltype);
UocvDesign=fOCVDesign(SOCavgTrue);

timeUsed=0;
maxStageViolation=0;
maxStageDeficitAh=0;
maxVoltageDifference=0;
worstActiveStageAh=NaN;
worstActiveStage=NaN;
worstActiveK=NaN;
for j=1:nStages
    Utrue=UocvTrue+Qno*stageCrate(j)*R0True;
    Udesign=UocvDesign+Qno*stageCrate(j)*R0Design;
    maxVoltageDifference=max(maxVoltageDifference,abs(Utrue-Udesign));

    duty=zeros(nCells,1);
    for level=1:nCells
        ratio=((2*level-1)*Utrue)/(2*UphaseCharge);
        ratio=min(max(ratio,0),1);
        duty(level)=(2/pi)*acos(ratio);
    end
    sumDuty=sum(duty);
    stageAh=sum(AhByStage(:,j));
    if stageAh<=tol
        continue
    end
    if sumDuty<=tol
        timeUsed=inf;
        maxStageViolation=inf;
        maxStageDeficitAh=inf;
        break
    end

    timeUsed=timeUsed+stageAh/(Qno*stageCrate(j)*sumDuty);

    sortedAh=sort(AhByStage(:,j),'descend');
    cumulativeRequested=cumsum(sortedAh);
    cumulativeAllowed=(cumsum(duty)/sumDuty)*stageAh;
    if nCells>1
        deficitVector=cumulativeRequested(1:end-1)- ...
            cumulativeAllowed(1:end-1);
        [stageDeficitAh,kWorst]=max(deficitVector);
        stageDeficitAh=max(stageDeficitAh,0);
    else
        stageDeficitAh=0;
        kWorst=1;
    end
    maxStageDeficitAh=max(maxStageDeficitAh,stageDeficitAh);

    % A normalized diagnostic is meaningful only for a stage carrying a
    % non-negligible fraction of the session throughput.
    if stageAh>activeStageAhTol
        stageViolation=stageDeficitAh/stageAh;
        if stageViolation>maxStageViolation
            maxStageViolation=stageViolation;
            worstActiveStageAh=stageAh;
            worstActiveStage=j;
            worstActiveK=kWorst;
        end
    end
end

SOCdisTrue=sum(QfinalTrue)/(2*sum(QcapTrue));
SOCdisTrue=min(max(SOCdisTrue,0),1);
UdisTrue=OCV(SOCdisTrue,celltype)-Qno*IdisNo*R0True;
UdisDesign=fOCVDesign(SOCdisTrue)-Qno*IdisNo*R0Design;
maxVoltageDifference=max(maxVoltageDifference,abs(UdisTrue-UdisDesign));

dprime=zeros(nCells,1);
for level=1:nCells
    ratio=((2*level-1)*UdisTrue)/(2*UdisPhase);
    ratio=min(max(ratio,0),1);
    dprime(level)=(2/pi)*acos(ratio);
end
if sum(dprime)<=tol || sum(QfinalTrue)<=tol
    dischargeViolation=inf;
else
    lhs=cumsum(sort(QfinalTrue,'descend'))/sum(QfinalTrue);
    rhs=cumsum(dprime)/sum(dprime);
    dischargeViolation=max(lhs(1:end-1)-rhs(1:end-1));
end

diag=struct();
diag.timeUsed=timeUsed;
diag.maxStageMajorizationViolation=max(maxStageViolation,0);
diag.maxStageMajorizationDeficitAh=max(maxStageDeficitAh,0);
diag.activeStageAhTolerance=activeStageAhTol;
diag.worstActiveStageAh=worstActiveStageAh;
diag.worstActiveStage=worstActiveStage;
diag.worstActiveK=worstActiveK;
diag.dischargeMajorizationViolation=max(dischargeViolation,0);
diag.maxVoltageProxyDifferenceV=maxVoltageDifference;
end

function stats=init_runtime_stats(strategyMode)
%INIT_RUNTIME_STATS Storage for solver, session, and mismatch diagnostics.
stats=struct();
stats.strategy=string(strategyMode);
stats.lpSolveTimes=zeros(0,1);
stats.acOptimizationTimes=zeros(0,1);
stats.lpSolvesPerACSession=zeros(0,1);
stats.maxVariables=0;
stats.maxInequalities=0;
stats.maxEqualities=0;
stats.maxModelBytes=0;
stats.capacitySafetyCorrections=0;
stats.maxCommandedSOCViolation=0;
stats.maxTrueStageMajorizationViolation=0;
stats.maxTrueStageMajorizationDeficitAh=0;
stats.maxTrueDischargeMajorizationViolation=0;
stats.maxTrueChargingTimeOverrunHours=0;
stats.maxVoltageProxyDifferenceV=0;
end

function stats=record_lp_runtime(stats,sol,prob)
%RECORD_LP_RUNTIME Append one LP solve and its model dimensions.
if isfield(sol,'solveTime') && isfinite(sol.solveTime)
    stats.lpSolveTimes(end+1,1)=sol.solveTime;
end
if isfield(prob,'index') && isfield(prob.index,'N')
    stats.maxVariables=max(stats.maxVariables,prob.index.N);
end
if isfield(prob,'nIneq')
    stats.maxInequalities=max(stats.maxInequalities,prob.nIneq);
end
if isfield(prob,'nEq')
    stats.maxEqualities=max(stats.maxEqualities,prob.nEq);
end
if isfield(prob,'modelBytes')
    stats.maxModelBytes=max(stats.maxModelBytes,prob.modelBytes);
end
end

function stats=finalize_runtime_stats(stats)
%FINALIZE_RUNTIME_STATS Derived summary statistics.
stats.lpCount=numel(stats.lpSolveTimes);
stats.acSessionCount=numel(stats.acOptimizationTimes);
if isempty(stats.lpSolveTimes)
    stats.medianLPSolveTime=NaN;
    stats.maxLPSolveTime=NaN;
else
    stats.medianLPSolveTime=median(stats.lpSolveTimes);
    stats.maxLPSolveTime=max(stats.lpSolveTimes);
end
if isempty(stats.acOptimizationTimes)
    stats.medianACOptimizationTime=NaN;
    stats.maxACOptimizationTime=NaN;
else
    stats.medianACOptimizationTime=median(stats.acOptimizationTimes);
    stats.maxACOptimizationTime=max(stats.acOptimizationTimes);
end
if isempty(stats.lpSolvesPerACSession)
    stats.medianLPSolvesPerACSession=NaN;
    stats.maxLPSolvesPerACSession=NaN;
else
    stats.medianLPSolvesPerACSession=median(stats.lpSolvesPerACSession);
    stats.maxLPSolvesPerACSession=max(stats.lpSolvesPerACSession);
end
end

function T=build_runtime_summary(strategyNames,statsList)
%BUILD_RUNTIME_SUMMARY Convert runtime structures to a CSV-ready table.
n=numel(statsList);
Strategy=strings(n,1);
LPCount=zeros(n,1);
ACSessionCount=zeros(n,1);
MedianLPSolve_ms=NaN(n,1);
MaxLPSolve_ms=NaN(n,1);
MedianLPSolvesPerACSession=NaN(n,1);
MaxLPSolvesPerACSession=NaN(n,1);
MedianACOptimization_s=NaN(n,1);
MaxACOptimization_s=NaN(n,1);
MaxLPVariables=zeros(n,1);
MaxLPInequalities=zeros(n,1);
MaxLPEqualities=zeros(n,1);
MaxStoredModel_MB=zeros(n,1);
SafetyCorrections=zeros(n,1);
MaxCommandedSOCViolation_pp=zeros(n,1);
MaxActiveStageMajorizationViolation=zeros(n,1);
MaxStageMajorizationDeficit_uAh=zeros(n,1);
MaxDischargeMajorizationViolation=zeros(n,1);
MaxChargingTimeOverrun_s=zeros(n,1);
MaxVoltageProxyDifference_mV=zeros(n,1);

for k=1:n
    st=statsList{k};
    Strategy(k)=string(strategyNames(k));
    LPCount(k)=st.lpCount;
    ACSessionCount(k)=st.acSessionCount;
    MedianLPSolve_ms(k)=1000*st.medianLPSolveTime;
    MaxLPSolve_ms(k)=1000*st.maxLPSolveTime;
    MedianLPSolvesPerACSession(k)=st.medianLPSolvesPerACSession;
    MaxLPSolvesPerACSession(k)=st.maxLPSolvesPerACSession;
    MedianACOptimization_s(k)=st.medianACOptimizationTime;
    MaxACOptimization_s(k)=st.maxACOptimizationTime;
    MaxLPVariables(k)=st.maxVariables;
    MaxLPInequalities(k)=st.maxInequalities;
    MaxLPEqualities(k)=st.maxEqualities;
    MaxStoredModel_MB(k)=st.maxModelBytes/1024^2;
    SafetyCorrections(k)=st.capacitySafetyCorrections;
    MaxCommandedSOCViolation_pp(k)=100*st.maxCommandedSOCViolation;
    MaxActiveStageMajorizationViolation(k)= ...
        st.maxTrueStageMajorizationViolation;
    MaxStageMajorizationDeficit_uAh(k)= ...
        1e6*st.maxTrueStageMajorizationDeficitAh;
    MaxDischargeMajorizationViolation(k)= ...
        st.maxTrueDischargeMajorizationViolation;
    MaxChargingTimeOverrun_s(k)=3600*st.maxTrueChargingTimeOverrunHours;
    MaxVoltageProxyDifference_mV(k)=1000*st.maxVoltageProxyDifferenceV;
end

T=table(Strategy,LPCount,ACSessionCount,MedianLPSolve_ms, ...
    MaxLPSolve_ms,MedianLPSolvesPerACSession,MaxLPSolvesPerACSession, ...
    MedianACOptimization_s,MaxACOptimization_s, ...
    MaxLPVariables,MaxLPInequalities,MaxLPEqualities,MaxStoredModel_MB, ...
    SafetyCorrections,MaxCommandedSOCViolation_pp, ...
    MaxActiveStageMajorizationViolation,MaxStageMajorizationDeficit_uAh, ...
    MaxDischargeMajorizationViolation,MaxChargingTimeOverrun_s, ...
    MaxVoltageProxyDifference_mV);
end
