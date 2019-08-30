function [err,P] = wrapper_dtb_rt_analytic(theta,rt,coh,choice,c,pars,plot_flag)
% changes: (1) check if pars has the bound type, and the integration method
% (2) changed -coh0 to coh0

%%
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
B0     = theta(4);
coh0   = theta(5);
y0a    = theta(6);
ndt_mu_delta = theta(7);

%%

if ~isempty(pars) && isfield(pars,'t')
    t = pars.t;
    dt = t(2)-t(1);
else
    dt = 0.005;
    t  = 0:dt:10;
end

t = t(:);

%%
Bup = B0;
drift = kappa * unique(coh + coh0);

yp = y0a/B0; % as a proportion of the bound height
P =  analytic_dtb(drift,t,Bup,yp);
% for legacy:
P.drift = drift;
P.Bup = B0;
P.Blo = -B0;
P.up.pdf_t = P.up.pdf_t';
P.lo.pdf_t = P.lo.pdf_t';


%% likelihood
if isfield(pars,'optim_method')
    optim_method = pars.optim_method;
else
    optim_method = 3;
end

switch optim_method
    case 1
        [err,pPred] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);
    case 2
        [err,pPred] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);
        % for incorrect trials, I just take the prob. or error (ignore RT)
        ignore_RT_incorrectTrials_flag = 1;
        if (ignore_RT_incorrectTrials_flag)
            ucoh = unique(coh);
            ncoh = length(ucoh);
            pError = nan(length(ucoh),1);
            pError(ucoh>0) = P.lo.p(ucoh>0);
            pError(ucoh<0) = P.lo.p(ucoh<0);
            pError(ucoh==0) = 0.5;
            for i=1:ncoh
                inds = c==0 & coh==ucoh(i);
                pPred(inds) = pError(i);
            end
            
            %clip
            pPred(pPred<eps) = eps;
            
            % logl = -sum(log(pPred));
            err = -nanmean(log(pPred));
        end
        
    case 3
        
        err = logl_choice_meanRT_1d(P,choice,rt,coh,c,ndt_m,ndt_m+ndt_mu_delta,ndt_s);
        
    case 4
        
        err = logl_choice_meanRT_1d(P,choice,rt,coh,c,ndt_m,ndt_m+ndt_mu_delta,ndt_s,'fit_choices',false);
        
    case 5 % MSE
        
        [~,mRT,eRT] = curva_media(rt,coh,c==1,0);
        udrift = P.drift;
        mean_dt = nan(size(udrift));
        mean_dt(udrift>0) = P.up.mean_t(udrift>0);
        mean_dt(udrift<0) = P.lo.mean_t(udrift<0);
        mean_dt(udrift==0) = 0.5 * (P.lo.mean_t(udrift==0) + P.up.mean_t(udrift==0)); % assume equal number of right and left

        ndtm = nan(size(udrift));
        ndtm(udrift>0) = ndt_m + ndt_mu_delta;
        ndtm(udrift<0) = ndt_m;
        ndtm(udrift==0) = ndt_m + ndt_mu_delta/2; % assume equal number of right and left

        
        err = mean((mean_dt + ndtm - mRT).^2);
        
        
end


%%
%% print
fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f B0=%.2f coh0=%.2f y0=%.2f ndt_delta=%.2f \n',...
    err,kappa,ndt_m,ndt_s,B0,coh0,y0a,ndt_mu_delta);

%%
if plot_flag
    
    figure(1);clf
    set(gcf,'Position',[263  338  377  563])
    subplot(2,1,1);
    curva_media(choice,coh,[],3);
    hold all
    ucoh = unique(coh);
    plot(ucoh,P.up.p,'k-');
    xlim([min(ucoh),max(ucoh)]);
    xlabel('Motion coherence');
    ylabel('P rightward choice')
    
    subplot(2,1,2);
    
    % only correct trials
    rt_model_c = P.up.mean_t;
    rt_model_c(ucoh<0) = P.lo.mean_t(ucoh<0);
    rt_model_c(ucoh==0) = (P.up.mean_t(ucoh==0)+P.lo.mean_t(ucoh==0))/2;
    rt_model_c = rt_model_c + ndt_m;
    
    curva_media(rt,coh,c==1,3);
    hold all
    plot(ucoh,rt_model_c,'k-');
    xlim([min(ucoh),max(ucoh)]);
    xlabel('Motion coherence');
    ylabel('RT (s)')
    
    %set(gcf,'Position',[270   793  1084   293])
    
    format_figure(gcf);
    
    drawnow
    
end