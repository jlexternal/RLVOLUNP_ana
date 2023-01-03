% RLVOLUNP Analyses
clear all
clc
% -------------- Input: --------------
%
nameofdataset = 'discovery'; % 'discovery' or 'replication'
%
% ------------------------------------
% default constants
condrgb = [41 52 29; 19 46 55; 58 22 22]/100;
condstr = {'Ref','V+','S+'};
condorder = [3 1 2]; % order of condition presentation in paper
parstr = {'learning rate','decay rate','learning noise','choice temperature'};
parstrgreek = {'alpha','delta','zeta','tau'};

% add directories
addpath(genpath('./data'));
addpath(genpath('./toolbox'));

% processed data:
% Note: condition order here is 1:Ref, 2:V+, 3:S+
load(sprintf('pcorrect_average_%s.mat',nameofdataset)); % loads pcor_raw  (subject,condition)
load(sprintf('pcorrect_curve_%s.mat',nameofdataset));   % loads pcor_cond (subject,trial,condition)
load(sprintf('prepeat_average_%s.mat',nameofdataset));  % loads prep_raw  (subject,condition)
load(sprintf('prepeat_curve_%s.mat',nameofdataset));    % loads prep_cond (subject,trial,condition)
load(sprintf('parameter_fits_%s.mat',nameofdataset));    % loads pars     (subject,parameter,condition)

% binomial test for exclusion criterion 
ntrialspercond = 160;
% Input: You may try different values of ncorrect (threshold = 89)
ncorrect = 90; % one lower will pass binomial test below
passBinomialTest = binocdf(ncorrect,ntrialspercond,.5,'upper') >= 0.05;
if passBinomialTest
    fprintf('%d/%d correct trials in any condition PASSES the binomial test FOR random choices.\n',...
        ncorrect,ntrialspercond);
else
    fprintf('%d/%d correct trials in any condition FAILS the binomial test FOR random choices.\n',...
        ncorrect,ntrialspercond);
end
pcorrect_threshold = ncorrect/ntrialspercond;
idx_incl = ~(any(pcor_raw < pcorrect_threshold,2) | any(isnan(pcor_raw),2));
nsubj = sum(idx_incl);

clearvars ncorrect pcorrect_threshold ntrialspercond passBinomialTest

%% Figure 1C,D (discovery) | Supplementary Figure 1A,B (replication)

% Figure 1C (left): Values of accuracy in each condition
fprintf('\nMedian values [& interquartile range] of participant accuracy:\n');
for icond = 1:3
    cond = condorder(icond);
    fprintf('%s: ',pad(condstr{cond},3));
    fprintf('%.2f [%.2f %.2f]\n',quantile(pcor_raw(idx_incl,cond),[.5 .25 .75]));
end
% Figure 1C (right): Accuracy curves in each condition
fprintf('\nMean values (& SEM) of participant accuracy curves\n');
for icond = 1:3
    cond = condorder(icond);
    fprintf('%s:\n',condstr{cond})
    fprintf(' %.2f  ',mean(pcor_cond(idx_incl,:,cond))');
    disp(' ');
    fprintf('(%.2f) ',std(pcor_cond(idx_incl,:,cond))/sqrt(nsubj));
    disp(' ')
end

% Figure 1D (left): Values of switch rate in each condition
fprintf('\nMedian values [& interquartile range] of participant switch rate:\n');
for icond = 1:3
    cond = condorder(icond);
    fprintf('%s: ',pad(condstr{cond},3));
    fprintf('%.2f [%.2f %.2f]\n',quantile(1-prep_raw(idx_incl,cond),[.5 .25 .75]));
end
% Figure 1D (right): Switch rate curves in each condition
fprintf('\nMean values (& SEM) of participant switch curves\n');
for icond = 1:3
    cond = condorder(icond);
    fprintf('%s:\n  ',condstr{cond})
    fprintf(' %.2f  ',mean(1-prep_cond(idx_incl,:,cond))');
    fprintf('\n  ');
    fprintf('(%.2f) ',std(1-prep_cond(idx_incl,:,cond))/sqrt(nsubj));
    disp(' ')
end

% Statistics for Figure 1C,D
cs = nchoosek(1:3,2);
for imeas = 1:2
    if imeas == 1
        fprintf('\nSigned-rank tests on p(correct)...\n');
    else
        fprintf('Signed-rank tests on p(switch)...\n');
    end
    for ic = 1:size(cs,1)
        icond = cs(ic,1);
        jcond = cs(ic,2);
        if jcond == 2
            jcond = icond;
            icond = 2;
        end
        switch imeas
        case 1
            x = pcor_raw(idx_incl,icond);
            y = pcor_raw(idx_incl,jcond);
        case 2
            x = 1-prep_raw(idx_incl,icond);
            y = 1-prep_raw(idx_incl,jcond);
        end
        [p,~,stats] = signrank(x,y);
        fprintf('%s vs %s: p=%.4f, z=%+.4f\n',condstr{icond},condstr{jcond},p,stats.zval);
    end
    disp(' ');
end

%% Figure 2

% See directory /models for model fitting code

%% Figure 3A (discovery) | Supplementary Figure 1C (replication)
clc
% Figure 3A
fprintf('\nMedian values [& interquartile range] of model fit parameters:\n');
for ipar = 1:4
    fprintf('%s\n',parstr{ipar});
    for icond = 1:3
        cond = condorder(icond);
        fprintf('  %s: ',pad(condstr{cond},3));
        fprintf('%.3f [%.3f %.3f]\n',quantile(pars(idx_incl,ipar,cond),[.5 .25 .75]));
    end
end

% Statistics for Figure 3A
fprintf('\nStatistics: condition-wise parameter differences (signed-rank tests)\n');
cs = nchoosek(1:3,2);
for ipar = 1:4
    fprintf('%s\n',parstr{ipar});
    parx = squeeze(pars(:,ipar,:));
    for ic = 1:size(cs,1)
        icond = cs(ic,1);
        jcond = cs(ic,2);
        if jcond == 2
            jcond = icond;
            icond = 2;
        end
        x = parx(idx_incl,icond);
        y = parx(idx_incl,jcond);
        
        [p,~,stats] = signrank(x,y);
        fprintf('  %s: p=%.4f, z=%+.4f\n',pad(sprintf('%s vs %s',condstr{icond},condstr{jcond}),9),p,stats.zval);
    end
    disp(' ');
end

% Figures 3B,C
%   See directory /models for model fitting code

clearvars cs ic icond imeas ipar jcond p parx stats x y
%% Figure 4

% Figure 4A: Correlation learning noise/choice temperature ~ behavior (discovery) | 
% Supplementary Figure 3 (replication) 
parscat = [pars(idx_incl,:,1);pars(idx_incl,:,2);pars(idx_incl,:,3)];
meascat = [pcor_raw(idx_incl,1);pcor_raw(idx_incl,2);pcor_raw(idx_incl,3)];
meascat = cat(2,meascat,1-[prep_raw(idx_incl,1);prep_raw(idx_incl,2);prep_raw(idx_incl,3)]);
measstr = {'accuracy','switch rate'};
for ipar = 3:4
    x = parscat(:,ipar);
    for imeas = 1:2
        y = meascat(:,imeas);
        [r,p] = corr(x,y,'Type','Spearman');
        fprintf('%s rho = %+.2f, p = %.4f\n',...
            pad(sprintf('Correlation (Spearman) %s - %s:',parstr{ipar},measstr{imeas}),56), ...
            r,p);
    end
end

splitstr = {'upper','lower'};
% Figure 4B,C
for ipar = 3:4
    % Figure 4B,C (bar plots): Values of accuracy in each condition split on
    % median value of parameter
    fprintf('\nMedian values [& interquartile range] of participant accuracy split on median value of %s:\n',parstr{ipar});
    for icond = 1:3
        cond = condorder(icond);
        idx_split = pars(:,ipar,cond) < median(pars(idx_incl,ipar,cond));
        fprintf('%s: Upper split          Lower split\n',pad(condstr{cond},3));
        for isplit = 1:2
            % isplit == 1/2 : upper split/lower split
            idx_split = ~idx_split;
            idx = idx_split & idx_incl;
            fprintf('     %.2f [%.2f %.2f]',quantile(pcor_raw(idx,cond),[.5 .25 .75]));
        end
        disp(' ');
    end
    
    % Figure 4B,C (curves): Accuracy curves in each condition split on median
    % value of parameter
    fprintf('\nMean values (& SEM) of participant accuracy curves split on median value of %s:\n',parstr{ipar});
    for icond = 1:3
        cond = condorder(icond);
        idx_split = pars(:,ipar,cond) < median(pars(idx_incl,ipar,cond));
        for isplit = 1:2
            idx_split = ~idx_split;
            idx = idx_split & idx_incl;
            fprintf('%s (%s split):\n',condstr{cond},splitstr{isplit});
            fprintf('  ');
            fprintf(' %.2f  ',mean(pcor_cond(idx,:,cond))');
            disp(' ');
            fprintf('  ');
            fprintf('(%.2f) ',std(pcor_cond(idx,:,cond))/sqrt(nsubj));
            disp(' ');
        end
    end
end

clearvars parscat meascat ipar imeas x y r p idx idx_split isplit splitstr measstr icond cond
%% Figure 5

% Figure 5A: Fit parameter correlation matrix (discovery) | 
% Supplementary Figure 4A (replication)
parscat = [pars(idx_incl,:,3) pars(idx_incl,:,1) pars(idx_incl,:,2)];
[r,p] = corr(parscat,parscat,'Type','Spearman');
r = tril(r);
r(logical(eye(size(r)))) = nan;
r(r==0) = nan;
% control for false discovery rate
[~,~,~,p_adj] = fdr_bh(p,0.05,'dep','yes');
figure(1);
clf
heatmap(r,'MissingDataColor','w','GridVisible','off','MissingDataLabel',' ','Colormap',parula)
clim([-.5 .5]);
title('Figure 5A');
fprintf('Adjusted p-values as shown on Figure 5A:\n');
p_adj = flip(flip(tril(p_adj),2)',2)

% Figure 5B: Null correlation matrix
xhat    = nan(size(pars));
rholim  = .05; % measure of decorrelatedness
for icond = 1:3
    x           = pars(idx_incl,:,icond); % original parameters
    toshuffle   = true;
    itry        = 0;
    % this can take some time depending on the value of rholim
    while toshuffle
        % shuffle parameter matrix
        for ipar = 1:size(pars,2)
            xhat(idx_incl,ipar,icond) = x(randperm(nsubj),ipar);
        end
        if isequal(~eye(4),abs(corr(xhat(idx_incl,:,icond),'Type','Spearman'))<rholim)
            toshuffle = false;
        end
        itry = itry + 1;
        if mod(itry,1e4) == 0
            fprintf('%d shuffles have been made...\n',itry);
        end
    end
end
fprintf('Parameter shuffling finished.\n');
parscat = [xhat(idx_incl,:,3) xhat(idx_incl,:,1) xhat(idx_incl,:,2)];
[r,p] = corr(parscat,parscat,'Type','Spearman');
r = tril(r);
r(logical(eye(size(r)))) = nan;
r(r==0) = nan;
% control for false discovery rate
[~,~,~,p_adj] = fdr_bh(p,0.05,'dep','yes');
figure(2);
clf
heatmap(r,'MissingDataColor','w','GridVisible','off','MissingDataLabel',' ','Colormap',parula);
clim([-.5 .5]);
title('Figure 5B');
fprintf('Adjusted p-values as shown on Figure 5B:\n');
p_adj = flip(flip(tril(p_adj),2)',2)

% Figure 5C: Covariations between conditions (discovery) | 
% Supplementary Figure 4B (replication)
for ipar = [1 3 4]
    x = pars(idx_incl,ipar,1); % Reference condition parameter
    fprintf('Covariation of %s between\n',parstr{ipar});
    for icond = 2:3
        fprintf('Reference & %s: ',condstr{icond});
        y = pars(idx_incl,ipar,icond);
        [r,p] = corr(x,y,'Type','Spearman');
        fprintf('rho = %.2f, p = %.4f\n',r,p);
    end
end

% Figure 5D: Covariations between parameters (discovery) | 
% Supplementary Figure 5 (replication)
parcombinations = [1 3; 1 4];
for icombi = 1:2
    fprintf('Covariation between %s and %s\n',parstr{parcombinations(icombi,1)},parstr{parcombinations(icombi,2)});
    for icond = 1:3
        cond = condorder(icond);
        if icond == 1 
            x = pars(idx_incl,parcombinations(icombi,1),:);
            y = pars(idx_incl,parcombinations(icombi,2),:);
            [r,p] = corr(x(:),y(:),'Type','Spearman');
            fprintf('%s(All): rho = %+.2f, p = %.4f\n',pad(condstr{cond},3),r,p)
        end
        x = pars(idx_incl,parcombinations(icombi,1),cond);
        y = pars(idx_incl,parcombinations(icombi,2),cond);
        [r,p] = corr(x,y,'Type','Spearman');
        fprintf('%s: rho = %+.2f, p = %.4f\n',pad(condstr{cond},8),r,p)
    end
end

clearvars cb cond icond ipar itry p p_adj parscat r rholim x xhat toshuffle y parcombinations icombi
%% Figure 6 (discovery) | Supplementary Figure 6 (replication)

parscat = [pars(idx_incl,:,3) pars(idx_incl,:,1) pars(idx_incl,:,2)];
% run PCA
fprintf('Running PCA...\n');
[coeffs,scores,~,~,expl] = pca(zscore(parscat,[],1));
coeffs_orig = coeffs;
scores_orig = scores;

% Figure 6A: Percent variance explained
fprintf('Percent variance of the parameters explained by \n');
for ipc = 1:3
    fprintf(' PC%d: %.2f\n',ipc,expl(ipc));
end

% sign coefficients by sign of 1st ingredient (alpha)
ind_negalpha = coeffs(1,:) < 0; 
ind_negalpha = repmat(ind_negalpha,[12 1]);
coeffs(ind_negalpha) = -coeffs(ind_negalpha);
idx_pca = scores(:,1) < 0;
scores(idx_pca,:) = -scores(idx_pca,:);

% r-squared values for each ingredient on PC1 and 2
r_pc_par = nan(2,size(parscat,2));
p_pc_par = nan(2,size(parscat,2));
for ipc = 1:2
    x = scores_orig(:,ipc);
    fprintf('PC %d\n',ipc);
    pstr = '';
    for ipar = 1:size(parscat,2)
        y = parscat(:,ipar);
        [r_pc_par(ipc,ipar),p_pc_par(ipc,ipar)] = corr(x,y,'type','Spearman');
        pstr = strcat(pstr,num2str(sprintf('%.04f ',p_pc_par(ipc,ipar))),',');
    end
    fprintf('%s\n\n',pstr);
end

% Figure 6B: Two dimensions of adaptation
% (left) Variance explained in each condition by each principal component
for ipc = 1:2
    fprintf('Percent variance of the parameters explained by PC%d in \n',ipc)
    fprintf('  %s    %s    %s\n',condstr{3},condstr{1},condstr{2})
    xpar = zscore(parscat,[],1);
    xhat = scores_orig(:,ipc)*coeffs_orig(:,ipc)';
    xres = xhat-xpar;
    vexp_cnd = [ ...
            1-sum(var(xres(:,1:4)))/sum(var(xpar(:,1:4))), ...
            1-sum(var(xres(:,5:8)))/sum(var(xpar(:,5:8))), ...
            1-sum(var(xres(:,9:12)))/sum(var(xpar(:,9:12)))];
    fprintf(' %.2f ',vexp_cnd*100);
    disp(' ');
end

% (right) Coefficients associated with each parameter and its r-squared value
for ipc = 1:2
    fprintf('Coefficients of PC%d (& r^2)\n',ipc);
    fprintf(['|           %s             |           %s             |            %s            |'...
        '\n'],condstr{3},condstr{1},condstr{2})
    for i = 1:3
        fprintf('%s    ',strjoin(parstrgreek,{'  ','  ','   '}));
    end
    disp(' ');
    fprintf('%+.2f  ',coeffs(:,ipc)');
    disp(' ');
    fprintf('(%.2f) ',r_pc_par(ipc,:).^2);
    disp(' ');
    % r-squared values
    fprintf('\n\n');
end

%% Figure 7

% see directory /models for model simulation code
