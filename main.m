clc
close all
clear

%% load data

crsp = load("crsp_vwret.mat");

csrpdates = crsp.crspdates;
[year,moth,day,date_num] = decode_date(csrpdates);

ret = crsp.ret;

figure 
plot(date_num, ret,'LineWidth',1.0)
datetick

%% GJR_GARCH(1,1)

% initial guess : mu, omega, alpha, gamma, beta
Theta0 = [mean(ret); 0.0002; 0.1; 0.3; 0.6];

% constriant 1
LB = [-1 0 0 0 0];

% constriant 2
A = [0 0 1 0.5 1];
b = 1;

opt = optimset('display','none');

Thetahat = fmincon(@(Theta) GJR_GARCH11lgl(Theta,ret),Theta0, ...
    A,b,[],[],LB,[],[],opt);

[L_GJR, ~, sigsqhat_GJR] = GJR_GARCH11lgl(Thetahat,ret);

% se
S = score_matrix(@GJR_GARCH11lgl,Thetahat,ret);
H = hessian_2sided(@GJR_GARCH11lgl,Thetahat,ret);
seQMLE = diag(sqrt(H\S/H));
zstat = Thetahat./seQMLE;

z_crit = norminv(1-0.05/2);

string = [cellstr(char(956)); cellstr(char(969));cellstr(char(945)) ... 
    ;cellstr(char(947));cellstr(char(946))];
t1 = table(string,Thetahat,seQMLE,zstat,'VariableNames', ...
    {'name' 'point estimate' 'standard error' 'zstat'});


%% GRACH(1,1)

% Use estimated parameters as starting values
Theta0 = [mean(ret); 0.0001; 0.1; 0.6];

% Constraints
% 1. Upper and lower bound of parameters
LB = [-1 0 0 0]; % -1 for lowest mean, 0 for omega, alpha, and beta

% 2. Sum constraints (alpha + beta) < 1
A = [0 0 1 1];
b = 1;

Theta = fmincon(@(Theta) restricted_GARCH11lgl(Theta,ret), Theta0, ...
    A,b,[],[],LB,[],[],opt);

[restricted_L, ~, sigsqhat_garch] = restricted_GARCH11lgl(Theta,ret);

% Likelihood ratio test
% H0: gamma = 0

LRT = 2*((-L_GJR) - (-restricted_L))

chi_crit = chi2inv(1-0.05,1)

% reject null

figure
subplot(2,1,1)
plot(date_num, sqrt(sigsqhat_GJR)*100,'LineWidth',1.0,'Color','red')
hold on
    plot(date_num, sqrt(sigsqhat_garch)*100,'LineWidth',1.0,'Color','blue')
hold off
legend('GJR-GARCH(1,1) model','GARCH(1,1) model')
ylabel('Monthly Volatility (%)')
xlabel('Time')
title('Volatility model comparsion')
datetick

subplot(2,1,2)
plot(date_num, ret*100,'LineWidth',1.0,'Color','blue')
legend("Return")
ylabel('Monthly return (%)')
xlabel('Time')
datetick



