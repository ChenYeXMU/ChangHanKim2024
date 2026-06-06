%% Replication file
% Negative income tax and universal basic income in the eyes of Aiyagari
% Yongsung Chang, Jong-Suk Han, and Sun-Bin Kim
% Macroeconomic Dynamics, 2024, 28:813-825.
% Ye Chen, Xiamen University, Feb 2026

clear; clc; close all;
mypath = 'C:\Program Files\MATLAB\VFIToolkit';
addpath(genpath(mypath))

%% Grid size
n_d = 101;
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
% prices
Params.r = 0.05; % initial guess for interest rate (is determined in GE)
Params.w = 1;    % initial guess for wage rate (is determined in GE)
% benchmark 
Params.tau = 0;
Params.Tr = 0;

%% Grids
% Deterministic steady state (no idiosyncratic risk)
%r_ss = 1/Params.beta - 1;
%K_ss = 

d_grid = linspace(0,1,n_d)';
a_grid = 20*linspace(0,1,n_a)'.^3;
Tauchen_q = 3;
[z_grid, pi_z] = discretizeAR1_Tauchen(0,Params.rho,Params.sigma_xi,n_z,Tauchen_q);
z_grid = exp(z_grid);
[mean_z,~,~,dist_z] = MarkovChainMoments(z_grid,pi_z);
z_grid = z_grid./mean_z;

%% ReturnFn
DiscountFactorParamNames = {'beta'};

ReturnFn = @(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr) benchmark_ReturnFn(h,aprime,a,z,r,w,sigma,B,gamma,tau,Tr)

%% Solve for value function and policy function for checking that things are working
disp('Test ValueFnIter')
vfoptions = struct();
simoptions = struct();
% vfoptions.divideandconquer = 0;
% vfoptions.gridinterplayer = 1;
% vfoptions.ngridinterp = 15;
% simoptions.gridinterplayer = vfoptions.gridinterplayer;
% simoptions.ngridinterp = vfoptions.ngridinterp;
tic;
[V, Policy] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
toc

%% Genera eqm variables
GEPriceParamNames = {'r','w'};

%% Set up FnsToEvaluate

FnsToEvaluate.K = @(h,aprime,a,z) a; % Aggregate Capital
FnsToEvaluate.L = @(h,aprime,a,z) z*h; % Aggregate Effective Labor

FnsToEvaluate0 = FnsToEvaluate;
FnsToEvaluate0.hours = @(h,aprime,a,z) h;  % working hours
FnsToEvaluate0.earnings = @(h,aprime,a,z,w) w*z*h; % labor earnings
FnsToEvaluate0.income = @(h,aprime,a,z,w,r) w*z*h+r*a; % before-tax income
FnsToEvaluate0.netincome = @(h,aprime,a,z,w,r,tau,Tr) (1-tau)*(w*z*h+r*a)+Tr; % after-tax income
FnsToEvaluate0.wealth = @(h,aprime,a,z) a;  % assets
FnsToEvaluate0.consumption = @(h,aprime,a,z,r,w,tau,Tr) (1-tau)*(w*z*h+r*a)+a+Tr-aprime; % consumption
FnsToEvaluate0.taxes = @(h,aprime,a,z,r,w,tau) tau*(w*z*h+r*a); % taxes

%% General eqm equations
GeneralEqmEqns.capitalmarket = @(r,alpha,delta,K,L) r-((1-alpha)*(L^alpha)*(K^(-alpha))-delta);
GeneralEqmEqns.labormarket = @(w,alpha,K,L) w-alpha*(L^(alpha-1)*K^(1-alpha));

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
AllStats = EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate0,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

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

% Quintile in Wealth
quintileK = [sum(AllStats.K.QuantileMeans(1:4)),sum(AllStats.K.QuantileMeans(5:8)),sum(AllStats.K.QuantileMeans(9:12)),sum(AllStats.K.QuantileMeans(13:16)),sum(AllStats.K.QuantileMeans(17:20))];
quintileKshare = quintileK./sum(quintileK).*100;

% Welfare
Wel_0 = V(1,ceil(n_z/2)).*StationaryDist(1,ceil(n_z/2));

%% Display results
disp('=======Benchmark=========')
% disp('=======Table 3. Steady states=========')
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
disp('======================================')

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
plot(a_grid,sum(StationaryDist, 2), 'LineWidth', 2)
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

categories = {'1st Quintile', '2nd Quintile', '3rd Quintile','4th Quintile', '5th Quintile'};
figure
b = bar(categorical(categories),quintileKshare,0.5);
grid on
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = compose('%.2f%%', quintileKshare);      
text(xtips, ytips, labels, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, 'FontWeight', 'bold');
%ylim([0 max(quintileKshare) * 1.15]);
%b.FaceColor = [0.3 0.6 0.9];


%% Reform for the UBI
Params.tau = 0.2;
Params.Tr = 0.080; % initial guess for transfer (is determined in GE)

%% Genera eqm variables
GEPriceParamNames = {'r','w','Tr'};

%% Set up FnsToEvaluate
FnsToEvaluate1 = FnsToEvaluate;
FnsToEvaluate1.taxrevenue = @(h,aprime,a,z,r,w,tau) tau*(w*z*h+r*a); % Tax revenue
FnsToEvaluate1.transfer = @(h,aprime,a,z,Tr) Tr; % % Transfers

%% General eqm equations
% GeneralEqmEqns1.capitalmarket = @(r,alpha,delta,K,L) r-((1-alpha)*(L^alpha)*(K^(-alpha))-delta);
% GeneralEqmEqns1.labormarket = @(w,alpha,K,L) w-alpha*(L^(alpha-1)*K^(1-alpha));
GeneralEqmEqns1 = GeneralEqmEqns;
GeneralEqmEqns1.govbudgetbalance = @(taxrevenue,transfer) taxrevenue-transfer;

%% Solve for stationary general equilibrium
tic;
[p_eqm1, GEcondns1] = HeteroAgentStationaryEqm_Case1(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate1,GeneralEqmEqns1,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
toc

% Update the equilibrium values
% Params.r = p_eqm.r;
% Params.w = p_eqm.w;
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm1.(GEPriceParamNames{pp});
end

%% Evaluate the model results
% Solve for value function and policy function
[V1, Policy1] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
% Values on the grid
ValuesOnGrid1 = EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy1,FnsToEvaluate1,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% Convert policy from indices to values
PolicyValues1 = PolicyInd2Val_Case1(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
% Stationary distribution
StationaryDist1 = StationaryDist_Case1(Policy1,n_d,n_a,n_z,pi_z,simoptions,Params);
% Aggregate variables
% AggVars = EvalFnOnAgentDist_AggVars_Case1(StationaryDist,Policy,FnsToEvaluate0,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% Calculate various stats
AllStats1 = EvalFnOnAgentDist_AllStats_Case1(StationaryDist1,Policy1,FnsToEvaluate0,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

%% Aggregate moments
agg.KK1 = AllStats1.K.Mean;  % capital
agg.LL1 = AllStats1.L.Mean;  % effective labor
agg.HH1= AllStats1.hours.Mean; % working hours
agg.ZZ1 = AllStats1.income.Mean; % income
agg.CC1 = AllStats1.consumption.Mean; % consumption
agg.YY1 = (agg.LL1^Params.alpha)*(agg.KK1^(1-Params.alpha)); % output
agg.II1 = Params.delta*agg.KK1; % investment
agg.TT1 = AllStats1.taxes.Mean; % taxes

% Calculate the Gini coefficient
gini.wealth1 = AllStats1.wealth.Gini;
gini.income1 = AllStats1.income.Gini;
gini.netincome1 = AllStats1.netincome.Gini;
gini.earnings1 = AllStats1.earnings.Gini;
gini.consumption1 = AllStats1.consumption.Gini;

% Coefficient of variation for hours, earnings, income and wealth
cv.hours1 = AllStats1.hours.StdDeviation/AllStats.hours.Mean;
cv.earnings1 = AllStats1.earnings.StdDeviation/AllStats.earnings.Mean;
cv.income1 = AllStats1.income.StdDeviation/AllStats.income.Mean;
cv.wealth1 = AllStats1.wealth.StdDeviation/AllStats.wealth.Mean;
cv.consumption1 = AllStats1.consumption.StdDeviation/AllStats.consumption.Mean;

% Quintile in Wealth
quintileK1 = [sum(AllStats1.K.QuantileMeans(1:4)),sum(AllStats1.K.QuantileMeans(5:8)),sum(AllStats1.K.QuantileMeans(9:12)),sum(AllStats1.K.QuantileMeans(13:16)),sum(AllStats1.K.QuantileMeans(17:20))];
quintileKshare1 = quintileK1./sum(quintileK1).*100;

% Welfare
Wel_1 = V1(1,ceil(n_z/2)).*StationaryDist1(1,ceil(n_z/2));
del_Wel = (Wel_1-Wel_0)/(-Wel_0);

% CEV
CEV = Wel_1-Wel_0;
CEV = exp(CEV)-1;

%% Display results
disp('=======UBI case=========')
fprintf('Tax rate          : %.3f \n',Params.tau)
fprintf('Transfer at y=0   : %.3f \n',Params.Tr)
fprintf('relative to y_bar : %.3f \n',agg.TT1/AllStats1.income.Mean)
fprintf('Output            : %.3f \n',agg.YY1)
fprintf('Capital           : %.3f \n',agg.KK1)
fprintf('Hours worked      : %.3f \n',agg.HH1)
fprintf('Effective labor   : %.3f \n',agg.LL1)
fprintf('Wage rate         : %.3f \n',Params.w)
fprintf('Interest rate     : %.3f \n',Params.r*100)
% fprintf('Income            : %f \n',agg.ZZ)
fprintf('Mean(y_bar)       : %.3f \n',AllStats1.income.Mean)
fprintf('Median            : %.3f \n',AllStats1.income.Median)
fprintf('Mean-Median ratio : %.3f \n',AllStats1.income.Mean/AllStats1.income.Median)
fprintf('Before-tax Gini   : %.3f \n',gini.income1)
fprintf('After-tax Gini    : %.3f \n',gini.netincome1)
fprintf('Wealth Gini       : %.3f \n',gini.wealth1)
fprintf('Tax/GDP           : %.3f \n',agg.TT1/agg.YY1)
disp('======================================')

%% Plots

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,PolicyValues1(2,:,1),'LineWidth',2)
plot(a_grid,PolicyValues1(2,:,n_z),'LineWidth',2)
hold off
legend('45 line','low z','high z')
grid on
title('Policy function a''(a,z)')
ylabel('Next-period assets')
xlabel('Current assets')

figure
plot(a_grid,PolicyValues1(1,:,1),'LineWidth',2)
hold on
plot(a_grid,PolicyValues1(1,:,n_z),'LineWidth',2)
hold off
legend('low z','high z')
grid on
title('Policy function h(a,z)')
ylabel('Working hours')
xlabel('Current assets')

% Asset distribution function
figure
plot(sum(StationaryDist1, 2), 'LineWidth', 2)
grid on
title('Probability distribution of assets')
ylabel('Probability')
xlabel('Current assets')
ylim([0,1])

% Asset cumulative distribution function
figure
plot(cumsum(sum(StationaryDist1, 2)), 'LineWidth', 2)
grid on
title('Cumulative distribution of assets')
ylabel('Probability')
xlabel('Current assets')
ylim([0,1])

categories = {'1st Quintile', '2nd Quintile', '3rd Quintile','4th Quintile', '5th Quintile'};
figure
b = bar(categorical(categories),quintileKshare1,0.5);
grid on
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = compose('%.2f%%', quintileKshare1);      
text(xtips, ytips, labels, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, 'FontWeight', 'bold');
%ylim([0 max(quintileKshare) * 1.15]);
%b.FaceColor = [0.3 0.6 0.9];
%% Optimal tax rate for UBI
Params.tau_opt = [0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20];
Wel_opt = zeros(1,length(Params.tau_opt));
CEV_opt = zeros(1,length(Params.tau_opt));
%%
for tt = 1:length(Params.tau_opt)
    disp(['iteration number', num2str(tt), '/', num2str(length(Params.tau_opt))]);
    disp('======================================')
    Params.tau = Params.tau_opt(tt);
    [p_eqm_opt, GEcondns_opt] = HeteroAgentStationaryEqm_Case1(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate1,GeneralEqmEqns1,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    [V_opt,Policy_opt] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist_opt = StationaryDist_Case1(Policy_opt,n_d,n_a,n_z,pi_z,simoptions,Params);
    Wel_opt(tt) = V_opt(1,ceil(n_z/2)).*StationaryDist_opt(1,ceil(n_z/2));
    CEV_opt(tt) = Wel_opt(tt)-Wel_0;
    CEV_opt(tt) = exp(CEV_opt(tt))-1;
end
%%
Wel_opt
CEV_opt

%% plot


%% Reform for the NIT
