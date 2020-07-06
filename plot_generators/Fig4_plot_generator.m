% Fig4_plot_generator.m
% This script uses complete_dataset,xlsx to generate figure 4
% Figure 4 - range dependent effects of Exp2

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
                                                                    
% example session
examsess_names = {'ST180403a';'ST180515d'}; % monkey D and monkey G 
           
%
h(1) = figure;
set(gcf,'position',[510 265 1250 1250], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['Figure 4. Range Dependent Effects'],'fontsize',10);axis off   
                                                                                        
for imonkey = 1:nmonkeys
    ind_imonkey = allmonkeys == imonkey;
    sessions_imonkey = allsessions(ind_imonkey);
    nsessions = length(sessions_imonkey);
    %
    % initialization
    nsessions_ana    = 0;
    sessionnames_ana = {};
    rho_all          = [];
    rho_OFF_all      = [];
    rho_ON_all       = [];
    rangeA_all       = [];
    rangeB_all       = [];
    
    %
    for isession = 1:nsessions
        
        nsessions_ana = nsessions_ana + 1;
        
        session = ['ST',sessions_imonkey{isession}];
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
        %
        sessionnames_ana(nsessions_ana,:) = {session};
          
        
        % %
        % plot the results - PART1: examples
        % %
        if ismember(session,examsess_names)
            figure(1)
            % Figure 4AB
            subplot(2,nmonkeys,imonkey)
            hold on
            set(gca,'fontsize',9, 'fontweight','normal')
            % stimOFF 
            p1 = plot_psyphy_spec(tableOFF, [beta0-beta2,beta1], 'o', 'k', 8);
            % stimON
            p2 = plot_psyphy_spec(tableON,  [beta0+beta2,beta1], '^', 'r', 8);
            %           
            lgd = legend([p1,p2],'stim OFF','stim ON','Location','NorthWest');
            lgd.FontSize = 14;           
            xtickangle(45);
            axis square
            xlabel('log(qB/qA)');
            ylabel('% B choices');
            text(2,0.5,{['\rho = ',num2str(relvalue,'%.1f')];...
                        ['\DeltaVA = ',num2str(rangeA*relvalue,'%.1f')];...
                        ['\DeltaVB = ',num2str(rangeB)];...
                        ['\delta\rho = ',num2str(relvalue_ON-relvalue_OFF,'%.2f')]},...
                        'fontsize', 10, 'color', 'k');  
        end % for ismember
        
    end % for isession
       
    
    % %
    % plot the results - PART2: population results
    % %
    % Figure 4CD
    subplot(2,nmonkeys,imonkey+nmonkeys)
    deltarange = rangeA_all .* rho_all - rangeB_all;
    deltarho   = rho_ON_all - rho_OFF_all;
    ind_exam = ismember(sessionnames_ana,examsess_names);
    %
    XXX = deltarange;
    YYY = deltarho;
    
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
    
    XXX = XXX(ind_goodsession);
    YYY = YYY(ind_goodsession);
    hold on
    plot(XXX,YYY,'ko','MarkerSize',8);
    plot(XXX(ind_exam),YYY(ind_exam),'ko','MarkerSize',8,'MarkerFaceColor','g');
    aa = polyfit(XXX,YYY,1);
    [RR_Spe_all,pp_Spe_all] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea_all,pp_Pea_all] = corr(XXX,YYY,'Type','Pearson');
    % XX  = [min(XXX),max(XXX)];
    XX = [-15,15];
    YY  = aa(1)*XX+aa(2);
    XX2 = [-15,15];
    YY2 = [-0.4,0.4];
    plot(XX,YY,'Color',[.7 .7 .7],'LineWidth',3);
    plot(XX2,[0 0],'k--','LineWidth',1);
    plot([0 0],YY2,'k--','LineWidth',1);
    text(-13,0.3,{['Pearson:'];...
                  ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.3f')];...
                  ['Spearman:'];...
                  ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.3f')]}...
                  , 'fontsize', 10, 'color', 'k');  
    axis([XX2,YY2]); axis square
    box off
    xlabel('\Delta VA - \Delta VB');
    ylabel('\rhostimON - \rhostimOFF');
    title({['monkey ',monkeynames{imonkey},', N=',num2str(length(XXX))]});
    
end % for imonkey



% % %
% function 1
% % %
function [psyphy, table01] = probitLink_regression(QtyA, QtyB, ChosenID, stim)
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
function [po] = plot_psyphy_spec(table01, beta, dots, clr, dotssize)

table01_all = table01;

table02 = table01(table01(:,1) & table01(:,2),:);
xx = log(table02(:,2)./table02(:,1));
yy = table02(:,3);

%domain for plot
x = min(xx)-log(2)/0.75 : .05 : max(xx)+log(2)/0.75;
%plot fitted sigmoid
y_choicefit = 1./(1+exp(-(beta(1)+beta(2)*x)));
plot(x,y_choicefit,'-','color',clr,'linewidth',3);
%plot choice pattern
po = plot(xx, yy, dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr); %,'markerfacecolor',clr)

%add forced choices
forcedAtab = table01_all( table01_all(:,1) & ~table01_all(:,2),:);
nfA = size(forcedAtab,1);
xx_forcedA = min(xx)-log(2)*(1:nfA);	xx_forcedA = sort(xx_forcedA)';
plot(xx_forcedA, forcedAtab(:,3)', dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr)

forcedBtab = table01_all(~table01_all(:,1) &  table01_all(:,2),:);
nfB = size(forcedBtab,1);
xx_forcedB = max(xx)+log(2)*(1:nfB);	xx_forcedB = sort(xx_forcedB)';
plot(xx_forcedB, forcedBtab(:,3)', dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr)

%xlims
set(gca,'xlim',[min(xx)-log(2)*(nfA+.5) max(xx)+log(2)*(nfB+.5)])

%cosmetics
[~,ind,~] = unique(xx);	%remove doubles in xx
xxx = xx(ind);

%xlabels
xlab = cell(1,nfA);
for ifA = 1:nfA
	xlab{ifA} = [num2str(forcedAtab(ifA,2)),':',num2str(forcedAtab(ifA,1))];
end
for i = 1:size(xxx,1)
	xlab{nfA+i} = [num2str(table02(ind(i),2)),':',num2str(table02(ind(i),1))];
end
for ifB = 1:nfB
	xlab{nfA+i+ifB} = [num2str(forcedBtab(ifB,2)),':',num2str(forcedBtab(ifB,1))];
end

%add forced choices
xxx = [xx_forcedA;xxx;xx_forcedB];

set(gca,'xtick',xxx,'xticklabel',xlab)
set(gca,'ylim',[0,1]);
set(gca,'ytick',0:.25:1,'yticklabel',0:25:100)

end