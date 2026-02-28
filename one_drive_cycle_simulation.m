clear;
close all;
clc;
%% Parameters
n_cell=10;%number of cells
SOC_0=0;%initial cell SOC
SOC_end=0.7;%Final pack SOC after charging
Uphase=2.5*n_cell;%phase voltage
Tlimit=3600;%time constraint in seconds
deltat=3600*6;%relaxation time in seconds
m=6;%SOC division number during CV stage
R0=0.01;%Internal resistance
temperature=25;
Pack_SOH_end=0.7;
%% battery initial states generation
rng(100)%set random seed 
Qmax=2;
Qmax1=Qmax*(0.9:0.1/(n_cell-1):1);
Qmax1=sort(Qmax1,'descend');
SOH_now=Qmax1./Qmax;
Qnow=Qmax1*SOC_0;
Tlimit_h=Tlimit/3600;
ideal=1;
%% Charging optimization
% Use maximum current, find the maximum voltage that makes charging possible
Unow = n_cell*OCV((SOC_0+SOC_end)/2,1);
Q_opt=[];
Iavg=min(I_CV(0.1), I_CV(SOC_0));
params.SOH_EOL=Pack_SOH_end;
while(isempty(Q_opt))
    SOC_CV=I_CV_inverse(Iavg);
    if(SOC_CV>=1)
        SOCgrid=ones(1,m+1);
    else
        SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
    end
    SOCgrid=max(SOCgrid,SOC_0);
    SOCgrid2=[SOC_0 SOCgrid];
    SOCs=zeros(1,m+1);
    for ii=1:m+1
        SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
    end
    Is=I_CV(SOCs);
    Is(1)=Iavg;
    [Q_opt, ~, ~] = solve_charge_allocation_lp( ...
        SOH_now, Qnow, Qmax.*SOH_now, ...
        SOC_end*sum(SOH_now)*Qmax, Is, SOCgrid, ...
        Unow, 'sinusoidal', @OCV1, R0, Qmax, Tlimit_h, ...
        Uphase, 1, params);
    Unow = Unow*0.95;
end

% Fix the voltage, find best current current
%Find minimum current (C rate)
d=real(acos(OCV(0.5*SOC_0+0.5*SOC_end,1)*(0:1:n_cell)/Unow))/pi*2;
%optimization
duty=zeros(1,n_cell);
for ii=1:n_cell
    duty(ii)=0.5*(d(ii)+d(ii+1));
end
Iavgmin=sum(SOH_now)*(SOC_end-SOC_0)/sum(duty)/Tlimit_h;
% search for optimal current
best=inf;
Ahbest=[];
while Iavg>Iavgmin*0.7
    SOC_CV=I_CV_inverse(Iavg);
    if(SOC_CV>=1)
        SOCgrid=ones(1,m+1);
    else
        SOCgrid=SOC_CV:(1-SOC_CV)/m:1;
    end
    SOCgrid=max(SOCgrid,SOC_0);
    SOCgrid2=[SOC_0 SOCgrid];
    SOCs=zeros(1,m+1);
    for ii=1:m+1
        SOCs(ii)=0.5*(SOCgrid2(ii)+SOCgrid2(ii+1));
    end
    Is=I_CV(SOCs);
    Is(1)=Iavg;
    [Q_opt2, sol, ~] = solve_charge_allocation_lp( ...
        SOH_now, Qnow, Qmax.*SOH_now, ...
        SOC_end*sum(SOH_now)*Qmax, Is, SOCgrid, ...
        Unow, 'sinusoidal', @OCV1, R0, Qmax, Tlimit_h, ...
        Uphase, 1, params);
    if sol.fval<best
        best=sol.fval;
        Ahbest=Q_opt2;
        Iavgbest=Iavg;
        SOCgridbest=SOCgrid;
        hpAh=sol.hour_per_Ah;
    end
    Iavg = Iavg*0.95;
end
tplot=zeros(1,m+2);
for ii=1:m+1
    tplot(ii+1)=sum(Ahbest(:,ii))*hpAh(ii);
end
Qplot=Ahbest;

t=tplot;
Qplotnew=Qplot;
for i=1:m+2
    t(i)=sum(tplot(1:i));
end
for i=1:m+1
    for j=1:n_cell
        Qplotnew(j,i)=sum(Qplot(j,1:i));
    end
end
Qplotnew=Qplotnew*3600;
t=t*3600;
%% discharge data
S = load('drive_cycle_data.mat');
Iline = S.Iline;
Uline = S.Uline;
clear S
t2=1:length(Iline);
%% SOC balancing
l1=length(t);
t=[t t(l1)+deltat t(l1)+deltat+t2];
Qplotnew=[zeros(n_cell,1),Qplotnew Qplotnew(:,end)];
for i=l1+1+1:length(t)
    Qplotnew=[Qplotnew zeros(n_cell,1)];
    if(ideal==1 && Iline(i-l1-1)<0)
        [~,idx]=sort(Qplotnew(:,i-1),'ascend');
    else
        [~,idx]=sort(Qplotnew(:,i-1),'descend');
    end
    vtemp=zeros(1,n_cell+1);
    duty=NaN(1,n_cell);
    for j=1:n_cell
        vtemp(j+1)=vtemp(j)+OCV1(Qplotnew(idx(j),i-1)/3600/Qmax1(idx(j)))/Uline(i-l1-1);
        if(vtemp(j+1)>1)
            vtemp(j+1)=1;
        end
        duty(j)=(acos(vtemp(j))*0.5+acos(vtemp(j+1))*0.5)/pi*2;
        Qplotnew(idx(j),i)=Qplotnew(idx(j),i-1)-Iline(i-l1-1)*duty(j);
    end
    if(sum(Qplotnew(:,i))<0)
        break
    end
end

%% Figures
f1 = figure(1);
f1.Position = [100 100 500 380];

ax = gca;          % grab axes handle once (creates axes if needed)
hold(ax,'on');
ax.FontSize = 14;

nLines = 15;
C = parula(nLines);

colororder(ax, C);
ax.ColorOrderIndex = 1;

labels=strings(1,n_cell);
for j=1:n_cell
    plot(t(1:i)/3600,100*Qplotnew(j,1:i)/3600/Qmax1(j))
    labels(j)=append('Cell ',int2str(j));
end
x1 = Tlimit_h;  % your x location
xline(x1, 'r--', 'LineWidth', 2);
x2 = (deltat+Tlimit)/3600;  % your x location
xline(x2, 'r--', 'LineWidth', 2);

xlabel('Time (h)','FontSize',14)
ylabel('SOC (%)','FontSize',14)

legend(labels,Location="northeast",FontSize=12)
xlim([0 8+deltat/3600])
ylim([0 100])
grid on
%exportgraphics(f1,'one_cycle.png','Resolution',900)
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
sol.hour_per_Ah = 1 ./ (Q_no .* I_no_avg .* sumd);

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