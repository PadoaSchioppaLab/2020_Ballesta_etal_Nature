% FigS3_plot_generator.m
% This script uses complete_dataset.xlsx to generate figure S3
% Figure S3 - Reaction time analysis of Exp2

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
allmonkeys  = cell2mat(Exp2_info(2:end,2)); % 1 for monkey D(acourt); 2 for monkey G(ervinho)
nmonkeys = 2; % length(unique(allmonkeys));
monkeynames = {'D','G'};
% from data file
allsessions_data = Exp2_data(:,1);

%
doremoveoutliers = 1;  % remove outliers of value range difference and relative value difference
limout=3;              % outlier threshold

% remove trials/trialtypes that has less than atleast_nntrials
atleast_nntrials = 3;

% example session
examsess_names = {'ST180414b';'ST180515d'}; % monkey D and monkey G 
                       
h(1) = figure;
set(gcf,'position',[510 265 1350 850], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['Figure S3. Reaction time analysis'],'fontsize',10);axis off  
h(2) = figure;
set(gcf,'position',[510 265 1050 650], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['reaction difference'],'fontsize',10);axis off  
                                                                                        
for imonkey = 1:nmonkeys
    ind_imonkey = allmonkeys == imonkey;
    sessions_imonkey = allsessions(ind_imonkey);
    nsessions = length(sessions_imonkey);
    %
    % initialization
    nsessions_ana    = 0;
    sessionnames_ana = {};
    %
    rho_OFF_all      = [];
    rho_ON_all       = [];
    chV_OFF_all      = {};
    chV_ON_all       = {};
    Rxt_OFF_all      = {};
    Rxt_ON_all       = {};
    %
    slope_OFF_all    = [];
    intcept_OFF_all  = [];
    slope_ON_all     = [];
    intcept_ON_all   = [];
    %
    meanRxt_OFF_all  = [];
    meanRxt_ON_all   = [];
    %
    rangeA_all       = [];
    rangeB_all       = [];
    
    %
    for isession = 1:nsessions
        session = ['ST',sessions_imonkey{isession}];
        ind_trials = ismember(allsessions_data,session);
        QtyA = cell2mat(Exp2_data(ind_trials,3)); % quantity of A
        QtyB = cell2mat(Exp2_data(ind_trials,4)); % quantity of B
        ChosenID = cell2mat(Exp2_data(ind_trials,5)); % chosen juice: 0 for juice A, 1 for juice B
        stim = cell2mat(Exp2_data(ind_trials,9)); % -1 for stimOFF, 1 for stimON
        Rxtime = cell2mat(Exp2_data(ind_trials,12)); % Reaction time
        
        % logstic regression
        warning off
        [psyphy, table01, ~] = probitLink_regression(QtyA, QtyB, ChosenID, stim); % table_trialtype is trialtype table
        tableOFF = table01.stimOFF;
        tableON  = table01.stimON;       
       
        nsessions_ana = nsessions_ana + 1;
        sessionnames_ana(nsessions_ana,:) = {session};
        
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
        
        % %
        % remove trials with outlier reaction time: for stimOFF and stimON separately
        % stimOFF
        meanRxt_OFF = nanmean(Rxtime(stim==-1));
        meanRxt_ON  = nanmean(Rxtime(stim== 1));
        upper_limit = mean(Rxtime(stim==-1))+4*std(Rxtime(stim==-1));
        lower_limit = mean(Rxtime(stim==-1))-4*std(Rxtime(stim==-1));
        outlier_ind = (Rxtime<=lower_limit | Rxtime>=upper_limit) & stim==-1;
        % stimON
        upper_limit = mean(Rxtime(stim== 1))+4*std(Rxtime(stim== 1));
        lower_limit = mean(Rxtime(stim== 1))-4*std(Rxtime(stim== 1));
        outlier_ind = outlier_ind | ((Rxtime<=lower_limit | Rxtime>=upper_limit) & stim == 1);
        QtyA(outlier_ind) = [];
        QtyB(outlier_ind) = [];
        ChosenID(outlier_ind) = [];
        stim(outlier_ind) = [];
        Rxtime(outlier_ind) = [];
                              
        % %
        % remove trialtypes that has less than atleast_nntrials
        [~, ~, table_trialtype] = probitLink_regression(QtyA, QtyB, ChosenID, stim); % table_trialtype is trialtype table
        tableOFF_trialtype = table_trialtype.stimOFF;
        tableON_trialtype  = table_trialtype.stimON;
        %
        ind_bad_OFF = tableOFF_trialtype(:,4)<atleast_nntrials;
        bad_trialtype_OFF = tableOFF_trialtype(ind_bad_OFF,1:3);
        ind_romove_OFF = ismember(tableOFF_trialtype(:,1:3),bad_trialtype_OFF,'rows');
        ind_romove_ON  = ismember(tableON_trialtype(:,1:3),bad_trialtype_OFF,'rows');
        tableOFF_trialtype(ind_romove_OFF,:) = [];
        tableON_trialtype(ind_romove_ON,:)   = [];
        %
        ind_bad_ON  = tableON_trialtype(:,4) <atleast_nntrials;
        bad_trialtype_ON = tableON_trialtype(ind_bad_ON,1:3);
        ind_romove_OFF = ismember(tableOFF_trialtype(:,1:3),bad_trialtype_ON,'rows');
        ind_romove_ON  = ismember(tableON_trialtype(:,1:3),bad_trialtype_ON,'rows');
        tableOFF_trialtype(ind_romove_OFF,:) = [];
        tableON_trialtype(ind_romove_ON,:)   = [];
        
        % %
        % calculate chosen value as in the trialtype table
        beta0 = psyphy.beta(1);
        beta1 = psyphy.beta(2);
        beta2 = psyphy.beta(3);
        relvalue = exp(-beta0/beta1);
        relvalue_OFF = exp(-(beta0-beta2)/beta1);
        relvalue_ON  = exp(-(beta0+beta2)/beta1);
                   
        % 
        % get chosen value as in trialtype table 
        chV_OFF_trialtype = tableOFF_trialtype(:,1).*relvalue;
        chV_OFF_trialtype(tableOFF_trialtype(:,3)==1) = tableOFF_trialtype(tableOFF_trialtype(:,3)==1,2);
        chV_ON_trialtype  = tableON_trialtype(:,1).*relvalue;
        chV_ON_trialtype(tableON_trialtype(:,3)==1)   = tableON_trialtype(tableON_trialtype(:,3)==1,2);
        %
        chV_OFF_all(nsessions_ana,:) = {chV_OFF_trialtype};
        chV_ON_all(nsessions_ana,:)  = {chV_ON_trialtype};
        
        %
        % calculate average reaction time as in trialtype table
        trialtype_all = [QtyA,QtyB,ChosenID,Rxtime];
        trialtype_OFF_all = trialtype_all(stim==-1,:);
        trialtype_ON_all  = trialtype_all(stim== 1,:);
        %
        ntrialtype_OFF = size(tableOFF_trialtype,1);
        Rxtime_OFF = [];
        for itrialtype = 1:ntrialtype_OFF
            trialtype = tableOFF_trialtype(itrialtype,1:3);
            ind_trials = ismember(trialtype_OFF_all(:,1:3),trialtype,'rows');
            Rxtime_OFF(itrialtype,:) = nanmean(trialtype_OFF_all(ind_trials,4));
        end
        %
        ntrialtype_ON = size(tableON_trialtype,1);
        Rxtime_ON = [];
        for itrialtype = 1:ntrialtype_ON
            trialtype = tableON_trialtype(itrialtype,1:3);
            ind_trials = ismember(trialtype_ON_all(:,1:3),trialtype,'rows');
            Rxtime_ON(itrialtype,:) = nanmean(trialtype_ON_all(ind_trials,4));
        end
        %
        Rxt_OFF_all(nsessions_ana,:) = {Rxtime_OFF};
        Rxt_ON_all(nsessions_ana,:)  = {Rxtime_ON};
        %
        meanRxt_OFF_all(nsessions_ana,:) = meanRxt_OFF;
        meanRxt_ON_all(nsessions_ana,:)  = meanRxt_ON;
        
        %
        % get intercept and slope of linear regression of chosen value v.s. reaction time
        mdl = fitlm(chV_OFF_trialtype, Rxtime_OFF);
        intcept_OFF = mdl.Coefficients.Estimate(1);
        slope_OFF   = mdl.Coefficients.Estimate(2);
        mdl = fitlm(chV_ON_trialtype, Rxtime_ON);
        intcept_ON  = mdl.Coefficients.Estimate(1);
        slope_ON    = mdl.Coefficients.Estimate(2);
        %
        slope_OFF_all(nsessions_ana,:)   = slope_OFF;
        intcept_OFF_all(nsessions_ana,:) = intcept_OFF;
        slope_ON_all(nsessions_ana,:)    = slope_ON;
        intcept_ON_all(nsessions_ana,:)  = intcept_ON;
        
                           
        % %
        % plot the results - PART1: example sessions
        % %
        if ismember(session,examsess_names)
            figure(1)
            % Figure S3AD
            subplot(nmonkeys,3,(imonkey-1)*3+1)
            hold on
            plot(chV_OFF_trialtype,Rxtime_OFF,'ko','markersize',8);
            plot(chV_ON_trialtype, Rxtime_ON, 'ro','markersize',8);
            XX = [min([chV_OFF_trialtype;chV_ON_trialtype]),max([chV_OFF_trialtype;chV_ON_trialtype])];
            YY_OFF = slope_OFF * XX + intcept_OFF;
            YY_ON  = slope_ON  * XX + intcept_ON ;
            plot(XX,YY_OFF,'k-','LineWidth',2);
            plot(XX,YY_ON, 'r-','LineWidth',2);
            xlabel('Chosen value (uB)');
            ylabel('Reaction time (ms)');
            title(['monkey ', monkeynames{imonkey}]);
            if imonkey == 1
                axis([0 10 140 180]);
            elseif imonkey == 2
                axis([0 14 150 185]);
            end
            legend({'stim OFF', 'stim ON'});
        end
        
    end % for isession
    
    % % %
    if doremoveoutliers  % remove outlier sessions with too large value range difference and relative value difference
        deltarange = rangeA_all .* rho_all - rangeB_all;
        deltarho   = rho_ON_all - rho_OFF_all;
        K=limout;
        OUT_x = (deltarange>(mean(deltarange)+std(deltarange)*K)) | (deltarange<(mean(deltarange)-std(deltarange)*K));
        OUT_y = (deltarho > (mean(deltarho)+std(deltarho)*K))     | (deltarho < (mean(deltarho)-std(deltarho)*K));
        %
        ind_goodsession = ~OUT_x & ~OUT_y;
    else
        ind_goodsession = logical(ones(size(deltarange)));
    end
    % % %
        
    
    % %
    % plot the results - PART2: polulation results
    % %
    % Figure S3BCEF
    subplot(nmonkeys,3,(imonkey-1)*3+2)
    XXX = intcept_OFF_all(ind_goodsession);
    YYY = intcept_ON_all(ind_goodsession);
    ind_exam = ismember(sessionnames_ana(ind_goodsession),examsess_names);
    %
    hold on
    plot(XXX,YYY,'ko','MarkerSize',8);
    plot(XXX(ind_exam),YYY(ind_exam),'ko','MarkerSize',8,'MarkerFaceColor','c');
    p_signrank = signrank(XXX,YYY);
    [~,p_ttest]= ttest(XXX,YYY);
    XX  = [140,200]; YY = XX;
    plot(XX,YY,'k--', 'LineWidth',1);
    text(145,195,{['mean difference = ',num2str(mean(YYY-XXX),'%.2f')];...
                  ['ttest p = ',num2str(p_ttest,'%.3f')];...
                  ['Wilcoxon p = ',num2str(p_signrank,'%.3f')]},...
                  'fontsize', 10, 'color', 'k');  
    axis([XX,YY]); axis square
    box off
    xlabel('Intercept, stimOFF (ms)');
    ylabel('Intercept, stimON (ms)');
    title(['monkey ', monkeynames{imonkey}]);
    %
    %
    subplot(nmonkeys,3,(imonkey-1)*3+3)
    XXX = slope_OFF_all(ind_goodsession);
    YYY = slope_ON_all(ind_goodsession);
    %
    hold on
    plot(XXX,YYY,'ko','MarkerSize',8);
    plot(XXX(ind_exam),YYY(ind_exam),'ko','MarkerSize',8,'MarkerFaceColor','c');
    p_signrank = signrank(XXX,YYY);
    [~,p_ttest]= ttest(XXX,YYY);
    XX  = [-4,3]; YY = XX;
    plot(XX,YY,'k--', 'LineWidth',1);
    plot([0 0],YY,'k--', 'LineWidth',1);
    plot(XX,[0 0],'k--', 'LineWidth',1);
    text(-3.8,2.5,{['mean difference = ',num2str(mean(YYY-XXX),'%.2f')];...
                  ['ttest p = ',num2str(p_ttest,'%.3f')];...
                  ['Wilcoxon p = ',num2str(p_signrank,'%.3f')]},...
                  'fontsize', 10, 'color', 'k');  
    axis([XX,YY]); axis square
    box off
    xlabel('Slope, stimOFF (ms)');
    ylabel('Slope, stimON (ms)');
    title(['monkey ', monkeynames{imonkey}]);
     
    %
    % figure 2
    % direct comparison of reaction time difference
    figure(2)
    subplot(1,nmonkeys,imonkey)
    XXX = meanRxt_OFF_all(ind_goodsession);
    YYY = meanRxt_ON_all(ind_goodsession);
    ind_exam = ismember(sessionnames_ana,examsess_names);
    %
    hold on
    plot(XXX,YYY,'ko','MarkerSize',8);
    plot(XXX(ind_exam),YYY(ind_exam),'ko','MarkerSize',8,'MarkerFaceColor','c');
    p_signrank = signrank(XXX,YYY);
    [~,p_ttest]= ttest(XXX,YYY);
    XX  = [140,200]; YY = XX;
    plot(XX,YY,'k--', 'LineWidth',1);
    text(145,195,{['mean difference = ',num2str(mean(YYY-XXX),'%.2f')];...
                  ['ttest p = ',num2str(p_ttest,'%.3f')];...
                  ['Wilcoxon p = ',num2str(p_signrank,'%.3f')]},...
                  'fontsize', 10, 'color', 'k');  
    axis([XX,YY]); axis square
    box off
    xlabel('meanRxTime, stimOFF (ms)');
    ylabel('meanRxTime, stimON (ms)');
    title(['monkey ', monkeynames{imonkey}]);
    
end % for imonkey



% % %
% function 1
% % %
function [psyphy, table01, table_trialtype] = probitLink_regression(QtyA, QtyB, ChosenID, stim)
%
% This function is to do logistic regression with "probit" link function
% This regression assumes stimON and stimOFF have the same steepness,
% therefore X = b0 + b1 * log(#B/#A) + b2 *(delta_stimON - delta_stimOFF)
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

% get trialtype tables
table_trialtype.stimOFF = get_tabletrialtype(QtyA(iOF,:), QtyB(iOF,:), ChosenID(iOF,:));
table_trialtype.stimON  = get_tabletrialtype(QtyA(iON,:), QtyB(iON,:), ChosenID(iON,:)); 

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



% % %
% function 3
% % %
function [table01] = get_tabletrialtype(QtyA, QtyB, ChosenID)
%
pairtrials = [QtyA, QtyB, ChosenID];
[offer, ~, groups] = unique(abs(pairtrials(:,1:3)),'rows');
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
table01 = [offer, Ntrials];
%
%sort choices
eps = 0.001;
aux = abs(table01) + eps;
[~, jnd] = sort(aux(:,2)./aux(:,1));
table01 = table01(jnd,:);

end
