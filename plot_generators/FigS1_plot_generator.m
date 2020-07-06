% FigS1_plot_generator.m
% This script uses complete_dataset.xlsx to generate figure S1
% Figure S1 - a control analysis for figure 4: a subset of data with ellipse threshold

% Author: Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

% load data and info files
% if use complete_dataset.xlsx
[~,~,Exp2_info] = xlsread('complete_dataset.xlsx', 4); % sheet #4
[~,~,Exp2_data] = xlsread('complete_dataset.xlsx', 5); % sheet #5
% from info file
allsessions = Exp2_info(2:end,1);
% from data file
allsessions_data = Exp2_data(:,1);
      
%
doremoveoutliers = 1;  % remove outliers of value range difference and relative value difference
limout=3;              % outlier threshold

%
h(1) = figure;
set(gcf,'position',[510 265 1250 550], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['Figure S1. Ellipse subset control'],'fontsize',10);axis off   
                                                                                        
%                                                                                    
nsessions = length(allsessions);
%
% initialization
nsessions_ana    = 0;
rho_all          = [];
rho_OFF_all      = [];
rho_ON_all       = [];
rangeA_all       = [];
rangeB_all       = [];
%
choiceA_prop_all = [];
%
% for demining regression
%
error_rho_all = [];
error_deltarho_all = [];
%
error_choiceA_prop_all = [];

%
for isession = 1:nsessions
    session = ['ST',allsessions{isession}];
    ind_trials = ismember(allsessions_data,session);
    QtyA = cell2mat(Exp2_data(ind_trials,3)); % quantity of A
    QtyB = cell2mat(Exp2_data(ind_trials,4)); % quantity of B
    ChosenID = cell2mat(Exp2_data(ind_trials,5)); % chosen juice: 0 for juice A, 1 for juice B
    stim = cell2mat(Exp2_data(ind_trials,9)); % -1 for stimOFF, 1 for stimON

    % logstic regression
    warning off
    [psyphy, table01] = probitLink_regression(QtyA, QtyB, ChosenID, stim);
    tableOFF = table01.stimOFF;
    tableON  = table01.stimON;

    nsessions_ana = nsessions_ana + 1;
    

    % calculate relative values and ranges
    beta0 = psyphy.beta(1);
    beta1 = psyphy.beta(2);
    beta2 = psyphy.beta(3);
    relvalue = exp(-beta0/beta1);
    relvalue_OFF = exp(-(beta0-beta2)/beta1);
    relvalue_ON  = exp(-(beta0+beta2)/beta1);
    rangeA = max(QtyA)-min(QtyA);
    rangeB = max(QtyB)-min(QtyB);
    %
    rho_all(nsessions_ana,:)     = relvalue;
    rho_OFF_all(nsessions_ana,:) = relvalue_OFF;
    rho_ON_all(nsessions_ana,:)  = relvalue_ON;
    rangeA_all(nsessions_ana,:)  = rangeA;
    rangeB_all(nsessions_ana,:)  = rangeB;     
    
    % calculate percentage of choice A
    trialnum_off = sum(tableOFF(:,4));
    trialnum_on  = sum(tableON(:,4));
    choiceBnum_off = sum(tableOFF(:,4).*tableOFF(:,3));
    choiceBnum_on  = sum(tableON(:,4).*tableON(:,3));
    choiceB_prop = (choiceBnum_off+choiceBnum_on)./(trialnum_off+trialnum_on);
    choiceA_prop = 1-choiceB_prop;
    %
    choiceA_prop_all(nsessions_ana,:) = choiceA_prop;
    
    % calculate errors for deming regression
    se_beta0 = psyphy.se(1);
    se_beta1 = psyphy.se(2);
    se_beta2 = psyphy.se(3);
    error_rho = abs(relvalue * se_beta0) + abs(beta0 * relvalue * se_beta1 / (beta1^2)); % d? = |? * d?0| + |?0 * ? /?1^2 * d?1| 
    error_rho_all(nsessions_ana,:) = error_rho;
    % error_deltarho = ?_ON - ?_OFF ~= |2?2 / ?1 * d?| + |2 * ? / ?1 * d?2| + |2? * ?2 / (?1^2) * d?1|
    error_deltarho = abs(2*beta2*error_rho/beta1)+...
                     abs(2*relvalue*se_beta2/beta1)+...
                     abs(2*relvalue*beta2*se_beta1/(beta1^2));
    error_deltarho_all(nsessions_ana,:) = error_deltarho;
    error_choiceA_prop_all(nsessions_ana,:) = 1/sqrt(trialnum_off+trialnum_on); % 1/sqrt(N)
    
end

% %
% plot the results
% %
% Figure S1
% step1: plot all sessions
deltarange = rangeA_all .* rho_all - rangeB_all;
deltarho   = rho_ON_all - rho_OFF_all;
error_deltarange = rangeA_all .* error_rho_all;

% % %
if doremoveoutliers  % remove outlier sessions with too large value range difference and relative value difference
    K=limout;
    OUT_x = (deltarange>(mean(deltarange)+std(deltarange)*K)) | (deltarange<(mean(deltarange)-std(deltarange)*K));
    OUT_y = (deltarho > (mean(deltarho)+std(deltarho)*K))     | (deltarho < (mean(deltarho)-std(deltarho)*K));
    %
    ind_goodsession = ~OUT_x & ~OUT_y;
else
    ind_goodsession = logical(ones(size(deltarange)));
end
% % %

%
YYY_plots = {'100*choiceA_prop_all', 'rho_all', 'deltarho'};
errorYYY_plots = {'error_choiceA_prop_all', 'error_rho_all', 'error_deltarho_all'};
YYY_labels = {'% choice A', 'Relative value (\rho)', '\rhostim ON - \rhostim OFF'};
YY2_ranges = [0 100; 0.5 5.5; -0.45 0.45];
YY3_dashlines = [50 50; 1 1; 0 0];
YY4_text = [-31; -1.15; -0.75];
%
for iplot = 1:3
    subplot(1,3,iplot)
    XXX = deltarange(ind_goodsession);
    eval(['YYY = ', YYY_plots{iplot},'(ind_goodsession);']);
    errorXXX = error_deltarange(ind_goodsession);
    eval(['errorYYY = ', errorYYY_plots{iplot},'(ind_goodsession);']);
    hold on
    plot(XXX,YYY,'ko','MarkerSize',6);
    [RR_Spe_all,pp_Spe_all] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea_all,pp_Pea_all] = corr(XXX,YYY,'Type','Pearson');
    [~,~,~,~,stats] = deming(XXX, YYY, mean(errorYYY)/mean(errorXXX),0.05); % deming regression
    slope_deming = stats.b_ci(2,:);
    intcept_deming = stats.b_ci(1,:);
    aa(1) = mean(slope_deming);
    aa(2) = mean(intcept_deming);
    XX  = [-15,15];
    YY  = aa(1)*XX+aa(2);
    YY2 = YY2_ranges(iplot,:);
    YY3 = YY3_dashlines(iplot,:);
    plot(XX,YY,'Color',[.7 .7 .7],'LineWidth',2);
    plot(XX,YY3,'k--','LineWidth',0.5);
    plot([0 0],YY2,'k--','LineWidth',0.5);
    YY4 = YY4_text(iplot,:);
    if pp_Pea_all<10^-5 & pp_Spe_all<10^-5
        text(-14,YY4,{['Pearson:'];...
                      ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.0e')];...
                      ['Spearman:'];...
                      ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.0e')];...
                      ['session # = ',num2str(length(XXX))]}...
                      , 'fontsize', 10, 'color', 'k');  
    else
        text(-14,YY4,{['Pearson:'];...
                      ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.4f')];...
                      ['Spearman:'];...
                      ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.4f')];...
                      ['session # = ',num2str(length(XXX))]}...
                      , 'fontsize', 10, 'color', 'k');  
    end
    axis([XX,YY2]); axis square
    box off
    xlabel('\Delta VA - \Delta VB');
    ylabel(YYY_labels{iplot});
end
%
% step2: plot subset sessions
% parameters for ellipse 
a_ell = 7; % theshold for choiceB percentage 43-57%
b_ell = 4.5; % threshold for value range difference
if a_ell>b_ell, c_ell = sqrt(a_ell*a_ell-b_ell*b_ell); e_ell = c_ell/a_ell; 
else, c_ell = sqrt(b_ell*b_ell-a_ell*a_ell); e_ell = c_ell/b_ell; end
if a_ell>b_ell, x1=0; x2=0; y1=50+a_ell; y2=50-a_ell;
else, x1=-b_ell; x2=b_ell; y1=50; y2=50;end
t_ell = -pi:0.01:pi;
xellipse = 0+b_ell*cos(t_ell);
yellipse = 50+a_ell*sin(t_ell);  
%
XXdata_ell = deltarange(ind_goodsession);
YYdata_ell = 100*choiceA_prop_all(ind_goodsession);
data_ell = XXdata_ell.^2/b_ell^2 + (YYdata_ell-50).^2/a_ell^2;
ind_thre = data_ell<=1;    
%
hold on
for iplot = 1:3
    subplot(1,3,iplot)
    %
    % plot the ellipse
    if iplot == 1
        hold on
        fill(xellipse,yellipse, [1 0.9 1], 'EdgeColor',[1 0.6 1]); 
        XXX = deltarange(ind_goodsession);
        eval(['YYY = ', YYY_plots{iplot},'(ind_goodsession);']);
        plot(XXX,YYY,'ko','MarkerSize',6);
    end
    %
    XXX = deltarange(ind_goodsession);
    XXX = XXX(ind_thre,:);
    eval(['YYY = ', YYY_plots{iplot},'(ind_goodsession);']);
    YYY = YYY(ind_thre,:);
    errorXXX = error_deltarange(ind_goodsession,:);
    errorXXX = errorXXX(ind_thre,:);
    eval(['errorYYY = ', errorYYY_plots{iplot},'(ind_goodsession);']);
    errorYYY = errorYYY(ind_thre,:);
    hold on
    plot(XXX,YYY,'mo','MarkerSize',6);
    [RR_Spe_all,pp_Spe_all] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea_all,pp_Pea_all] = corr(XXX,YYY,'Type','Pearson');
    [~,~,~,~,stats] = deming(XXX, YYY, mean(errorYYY)/mean(errorXXX),0.05); % deming regression
    slope_deming = stats.b_ci(2,:);
    intcept_deming = stats.b_ci(1,:);
    aa(1) = mean(slope_deming);
    aa(2) = mean(intcept_deming);
    XX  = [-15,15];
    YY  = aa(1)*XX+aa(2);
    YY2 = YY2_ranges(iplot,:);
    YY3 = YY3_dashlines(iplot,:);
    plot(XX,YY,'Color',[1 .8 1],'LineWidth',2);
    plot(XX,YY3,'k--','LineWidth',0.5);
    plot([0 0],YY2,'k--','LineWidth',0.5);
    YY4 = YY4_text(iplot,:);
    if pp_Pea_all<10^-5 & pp_Spe_all<10^-5
        text(1,YY4,{['Pearson:'];...
                      ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.0e')];...
                      ['Spearman:'];...
                      ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.0e')];...
                      ['session # = ',num2str(length(XXX))]}...
                      , 'fontsize', 10, 'color', 'm');  
    else
        text(1,YY4,{['Pearson:'];...
                      ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.4f')];...
                      ['Spearman:'];...
                      ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.4f')];...
                      ['session # = ',num2str(length(XXX))]}...
                      , 'fontsize', 10, 'color', 'm');  
    end
    axis([XX,YY2]); axis square
    box off
end



% % %
% function 1
% % %
function [psyphy, table01] = probitLink_regression(QtyA, QtyB, ChosenID, stim)
%
% This function is to do logistic regression with "probit" link function
% This regression assumes stimON and stimOFF have the same steepness,
% therefore X = ?0 + ?1 * log(#B/#A) + ?2 *(?_stimON - ?_stimOFF)
%
% inputs of this function are:
%   QtyA:     quantity of juice A
%   QtyB:     quantity of juice B
%   ChosenID: chosen juice, 0 for juice A, 1 for juice B
%   stim:     stimulate or not, -1 for stimOFF, 1 for stimON

% outputs of this function are structures "psyphy" and "table01":
%   psyphy.vars:     names of varibles
%   psyphy.bera:     betas of regression
%   psyphy.se:       standard errors of betas
%   psyphy.pval:     p-values of betas
%   table01.stimOFF: decision table for stimulation OFF
%   table01.stimON:  decision table for stimulation ON

% remove forced choices
ii = ~prod([QtyA,QtyB],2)==0;
%
iOF = stim==-1;
iON = stim==1;
%
logB_A = log(abs(QtyB)./QtyA);
%
% logistic regression
jj = ii;
[~ , ~, stats] = glmfit([logB_A(jj,1),stim(jj,1)], ChosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
psyphy.vars	= {'constant','logit(#B/#A)','stimEffect'};
psyphy.beta = stats.beta;
psyphy.se   = stats.se;
psyphy.pval = stats.p;
%
% get decision tables
table01.stimOFF = get_table01(QtyA(iOF,:), QtyB(iOF,:), ChosenID(iOF,:));
table01.stimON  = get_table01(QtyA(iON,:), QtyB(iON,:), ChosenID(iON,:)); 

end



% % %
% function 2
% % %
function [table01] = get_table01(QtyA, QtyB, ChosenID)
%
pairtrials = [QtyA, QtyB, ChosenID];
[offer, ~, groups] = unique(abs(pairtrials(:,1:2)),'rows');
noffs = size(offer,1);
choiz = pairtrials(:,3);
%
perc_B = nan(noffs,1);
Ntrials = nan(noffs,1);
for i = 1:noffs
	ind = find(groups== i);
	perc_B(i,1) = length(find(choiz(ind)== 1))/length(ind);	%perc of ch b
	Ntrials(i,1) = length(ind);
end
table01 = [offer, perc_B, Ntrials];
%
%sort choices
eps = 0.001;
aux = abs(table01) + eps;
[~, jnd] = sort(aux(:,2)./aux(:,1));
table01 = table01(jnd,:);

end



