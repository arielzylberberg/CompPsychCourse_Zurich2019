addpath(genpath('../functions_addtopath/'))

load ../test_data/test_data.mat

%% run one set of parameters
kappa = 15;
ndt_m = 0.2;
ndt_s = 0.01;
B0    = 1.2;
coh0  = 0;
y0    = 0;
ndt_m_delta = 0;
plot_flag = 1;
pars = [];
theta = [kappa,ndt_m,ndt_s,B0,coh0,y0,ndt_m_delta];
[err,P] = wrapper_dtb_rt_analytic(theta,rt,coh,choice,c,pars,plot_flag);


%% fit

% kappa, ndt_mu, ndt_sigma, B0, coh0, y0, ndt_m_delta
tl = [5,  0.1, 0 ,.5  ,0,0,0];
th = [40, 0.5, 0 ,3  ,0,0,0];
tg = [15, 0.2, 0 ,1  ,0,0,0];

plot_flag = true;
pars = struct('optim_method',5); % MSE
% pars = struct('optim_method',4); % fit mean reaction times, predict RT
% pars = struct('optim_method',3); % fit means choice & RT

c(coh==0) = 1; % to use all 0% coh trials in the fits
fn_fit = @(theta) (wrapper_dtb_rt_analytic(theta,rt,coh,choice,c,pars,plot_flag));

options = optimset('Display','final','TolFun',.01,'FunValCheck','on');
ptl = tl;
pth = th;
[theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

% save fitted params
save('fit_output','theta','fval','exitflag','output','tl','th','tg','pars');