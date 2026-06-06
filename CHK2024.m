%% Replication file
% Negative income tax and universal basic income in the eyes of Aiyagari
% Yongsung Chang, Jong-Suk Han, and Sun-Bin Kim
% Macroeconomic Dynamics, 2024, 28:813-825.
% Ye Chen, Xiamen University, Feb 2026

clear; clc; close all;
mypath = 'C:\Program Files\MATLAB\VFIToolkit';
addpath(genpath(mypath))

%% Scenarios
Params.Reform = 0; % 0: Benchmark
                   % 1: UBI reform
                   % 2: NIT reform

%% Grid size
n_d = 51;
n_a = 601;
n_z = 11;

%% Parameters
Params.beta = 0.94; % time discount factor
Params.alpha = 0.64; % labor share in Cobb-Douglas prod fun
Params.delta = 0.10; % depreciation rate
Params.sigma = 1.00; % relative risk aversion in utility, so actually is log(c)
Params.B = 4.5;   % disutility from working (I failed to replicate the results by B=22.5 in CHK paper,)
Params.gamma = 0.50; % labor supply elasticity
Params.rho = 0.91; % persistence of idiosyncratic productivity
Params.sigma_xi = 0.21; % std dev of innovation to prodctivity shock
% tax and transfer
if Params.Reform == 0
    Params.tau = 0;  % income tax rate
    Params.Tr = 0;   % initial guess for transfer (is determined in GE)
elseif Params.Reform == 1
    Params.tau = 0.20;
    Params.Tr = 0.080;
end
% prices
Params.r = 0.05; % initial guess for interest rate (is determined in GE)
Params.w = 1;    % initial guess for wage rate (is determined in GE)

%% Grids
% Deterministic steady state (no idiosyncratic risk)
%r_ss = 1/Params.beta - 1;
%K_ss = 

d_grid = linspace(0,1,n_d)';
a_grid = 30*linspace(0,1,n_a)'.^3;
Tauchen_q = 3;
[z_grid, pi_z] = discretizeAR1_Tauchen(0,Params.rho,Params.sigma_xi,n_z,Tauchen_q);
z_grid = exp(z_grid);
[mean_z,~,~,dist_z] = MarkovChainMoments(z_grid,pi_z);
z_grid = z_grid./mean_z;

%% ReturnFn
DiscountFactorParamNames = {'beta'};

ReturnFn = @(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr,Reform) CHK2024_ReturnFn(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr,Reform)

%% Solve for value function and policy function for checking that things are working
disp('Test ValueFnIter')
vfoptions = struct();
simoptions = struct();
%vfoptions.divideandconquer = 0;
%vfoptions.gridinterplayer = 1;
%vfoptions.ngridinterp = 15;
%simoptions.gridinterplayer = vfoptions.gridinterplayer;
%simoptions.ngridinterp = vfoptions.ngridinterp;
tic;
[V, Policy] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
toc

%% Genera eqm variables
if Params.Reform == 0
    GEPriceParamNames = {'r','w'};
elseif Params.Reform == 1
    GEPriceParamNames = {'r','w','Tr'};
end

%% Set up FnsToEvaluate

if Params.Reform == 0
    FnsToEvaluate.K = @(h,aprime,a,z) a; % Aggregate Capital
    FnsToEvaluate.L = @(h,aprime,a,z) z*h; % Aggregate Effective Labor
elseif Params.Reform == 1
    FnsToEvaluate.K = @(h,aprime,a,z) a; % Aggregate Capital
    FnsToEvaluate.L = @(h,aprime,a,z) z*h; % Aggregate Effective Labor
    FnsToEvaluate.taxrevenue = @(h,aprime,a,z,w,tau) tau*w*z*h; % Tax revenue
    FnsToEvaluate.transfer = @(h,aprime,a,z,Tr) Tr; % % Transfers
end

FnsToEvaluate2 = FnsToEvaluate;
FnsToEvaluate2.hours = @(h,aprime,a,z) h;  % working hours
FnsToEvaluate2.earnings = @(h,aprime,a,z,w) w*z*h; % labor earnings
FnsToEvaluate2.income = @(h,aprime,a,z,w,r) w*z*h+r*a; % before-tax income
FnsToEvaluate2.netincome = @(h,aprime,a,z,w,r,tau) (1-tau)*(w*z*h+r*a); % after-tax income
FnsToEvaluate2.wealth = @(h,aprime,a,z) a;  % assets
FnsToEvaluate2.consumption = @(h,aprime,a,z,r,w,tau,Tr) (1-tau)*(w*z*h+r*a)+a+Tr-aprime; % consumption
FnsToEvaluate2.taxes = @(h,aprime,a,z,r,w,tau) tau*(w*z*h+r*a); % taxes

%% General eqm equations
if Params.Reform == 0
    GeneralEqmEqns.capitalmarket = @(r,alpha,delta,K,L) r-((1-alpha)*(L^alpha)*(K^(-alpha))-delta);
    GeneralEqmEqns.labormarket = @(w,alpha,K,L) w-alpha*(L^(alpha-1)*K^(1-alpha));
elseif Params.Reform == 1
    GeneralEqmEqns.capitalmarket = @(r,alpha,delta,K,L) r-((1-alpha)*(L^alpha)*(K^(-alpha))-delta);
    GeneralEqmEqns.labormarket = @(w,alpha,K,L) w-alpha*(L^(alpha-1)*K^(1-alpha));
    GeneralEqmEqns.govbudgetbalance = @(taxrevenue,transfer) taxrevenue-transfer;
end

%% Solve for stationary general equilibrium
heteroagentoptions.verbose = 1;
heteroagentoptions.toleranceGEprices = 10^(-5);
heteroagentoptions.toleranceGEcondns = 10^(-5);
tic;
[p_eqm, GEcondns] = HeteroAgentStationaryEqm_Case1(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
toc

% Update the equilibrium values
% Params.r = p_eqm.r;
% Params.w = p_eqm.w;
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

%% Evaluate the model results
% Solve for value function and policy function
[V, Policy] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
% Values on the grid
ValuesOnGrid = EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% Convert policy from indices to values
PolicyValues = PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
% Stationary distribution
StationaryDist = StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
% Aggregate variables
% AggVars = EvalFnOnAgentDist_AggVars_Case1(StationaryDist,Policy,FnsToEvaluate2,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% Calculate various stats
AllStats = EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate2,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

%% Aggregate moments
agg.KK = AllStats.K.Mean;  % capital
agg.LL = AllStats.L.Mean;  % effective labor
agg.HH = AllStats.hours.Mean; % working hours
agg.ZZ = AllStats.income.Mean; % income
agg.CC = AllStats.consumption.Mean; % consumption
agg.YY = (agg.LL^Params.alpha)*(agg.KK^(1-Params.alpha)); % output
agg.II = Params.delta*agg.KK; % investment
agg.TT = AllStats.taxes.Mean; % taxes

% Calculate the Gini coefficient
gini.wealth = AllStats.wealth.Gini;
gini.income = AllStats.income.Gini;
gini.netincome = AllStats.netincome.Gini;
gini.earnings = AllStats.earnings.Gini;
gini.consumption = AllStats.consumption.Gini;

% Coefficient of variation for hours, earnings, income and wealth
cv.hours = AllStats.hours.StdDeviation/AllStats.hours.Mean;
cv.earnings = AllStats.earnings.StdDeviation/AllStats.earnings.Mean;
cv.income = AllStats.income.StdDeviation/AllStats.income.Mean;
cv.wealth = AllStats.wealth.StdDeviation/AllStats.wealth.Mean;
cv.consumption = AllStats.consumption.StdDeviation/AllStats.consumption.Mean;

%% Display results
if Params.Reform == 0
    disp('===========Benchmark=============')
elseif Params.Reform == 1
    disp('===========UBI reform============')
end 
% disp('==============Table 3. Steady states==================')
fprintf('Tax rate          : %.3f \n',Params.tau)
fprintf('Transfer at y=0   : %.3f \n',Params.Tr)
fprintf('relative to y_bar : %.3f \n',agg.TT/AllStats.income.Mean)
fprintf('Output            : %.3f \n',agg.YY)
fprintf('Capital           : %.3f \n',agg.KK)
fprintf('Hours worked      : %.3f \n',agg.HH)
fprintf('Effective labor   : %.3f \n',agg.LL)
fprintf('Wage rate         : %.3f \n',Params.w)
fprintf('Interest rate     : %.3f \n',Params.r*100)
% fprintf('Income            : %f \n',agg.ZZ)
fprintf('Mean(y_bar)       : %.3f \n',AllStats.income.Mean)
fprintf('Median            : %.3f \n',AllStats.income.Median)
fprintf('Mean-Median ratio : %.3f \n',AllStats.income.Mean/AllStats.income.Median)
fprintf('Before-tax Gini   : %.3f \n',gini.income)
fprintf('After-tax Gini    : %.3f \n',gini.netincome)
fprintf('Wealth Gini       : %.3f \n',gini.wealth)
fprintf('Tax/GDP           : %.3f \n',agg.TT/agg.YY)
disp('==================================')

%% Plots

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,PolicyValues(2,:,1),'LineWidth',2)
plot(a_grid,PolicyValues(2,:,n_z),'LineWidth',2)
hold off
legend('45 line','low z','high z')
grid on
title('Policy function a''(a,z)')
ylabel('Next-period assets')
xlabel('Current assets')

figure
plot(a_grid,PolicyValues(1,:,1),'LineWidth',2)
hold on
plot(a_grid,PolicyValues(1,:,n_z),'LineWidth',2)
hold off
legend('low z','high z')
grid on
title('Policy function h(a,z)')
ylabel('Working hours')
xlabel('Current assets')

% Asset distribution function
figure
plot(sum(StationaryDist, 2), 'LineWidth', 2)
grid on
title('Probability distribution of assets')
ylabel('Probability')
xlabel('Current assets')
ylim([0,1])

% Asset cumulative distribution function
figure
plot(cumsum(sum(StationaryDist, 2)), 'LineWidth', 2)
grid on
title('Cumulative distribution of assets')
ylabel('Probability')
xlabel('Current assets')
ylim([0,1])









