% FigS2_plot_generator.m
% This script uses complete_dataset.xlsx to generate figure S2
% Figure S2 - session pairs

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

% example session
examsess_names = {'ST180517a','ST180517b'};
nexamples = 2;

h(1) = figure;
set(gcf,'position',[510 265 1250 1250], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['Figure S2. Range dependent effects in session pairs'],'fontsize',10);axis off   

%
% remove sessions first
if doremoveoutliers
    nsessions = length(allsessions);
    allsessions_ana = {};
    nsessions_ana    = 0;
    rho_all          = [];
    rho_OFF_all      = [];
    rho_ON_all       = [];
    rangeA_all       = [];
    rangeB_all       = [];

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
    end
    
    deltarange = rangeA_all .* rho_all - rangeB_all;
    deltarho   = rho_ON_all - rho_OFF_all;
    K=limout;
    OUT_x = (deltarange>(mean(deltarange)+std(deltarange)*K)) | (deltarange<(mean(deltarange)-std(deltarange)*K));
    OUT_y = (deltarho > (mean(deltarho)+std(deltarho)*K))     | (deltarho < (mean(deltarho)-std(deltarho)*K));
    %
    ind_goodsession = ~OUT_x & ~OUT_y;
    allsessions = allsessions(ind_goodsession);   
end

%
%
% find paired sessions
% CAUTION!: these sessions pairs are identified manually based on lab notebook instead of automatically
badsessions = {'ST180405a','ST180510a','ST180515a','ST180518a','ST180524a','ST180414c'}; % these sessions do not have consecutive sessions
sessions_withpairs = {'inta','intb'}; % initialization: inta - initialized a session; intb - initialized b session
nsessions = length(allsessions);
for isession = 1:nsessions
    session_target = ['ST',allsessions{isession}];
    %
    if ~ismember(session_target,sessions_withpairs) % avoid count one session twice
        if isequal(session_target(9),'a') & ~ismember(session_target,badsessions)
            session_next=[session_target(1:8),'b'];
        elseif isequal(session_target(9),'c') & ~ismember(session_target,badsessions)
            session_next=[session_target(1:8),'d'];
        elseif isequal(session_target(9),'e') & ~ismember(session_target,badsessions)
            session_next=[session_target(1:8),'f'];
        elseif isequal(session_target,'ST180405b')
            session_next='ST180405c';
        elseif isequal(session_target,'ST180518d')
            session_next='ST180518e';
        elseif isequal(session_target,'ST180522d')
            session_next='ST180522e';
        else
            continue
        end
        %
        if ismember(session_next(3:9),allsessions)
            if ismember('inta',sessions_withpairs)
                sessions_withpairs{1,1} = session_target;
                sessions_withpairs{1,2} = session_next;
            else
                sessions_withpairs{end+1,1} = session_target;
                sessions_withpairs{end,2} = session_next;
            end
        end   
    end
end

%
%
% generate data and plot figures
%
% initialization
nsessionpairs_ana = 0;
sessionpairnames_ana = {};
rho_all           = [];
rho_OFF_all       = [];
rho_ON_all        = [];
rangeA_all        = [];
rangeB_all        = [];
% 
iexample          = 0;
%
nsessionpairs = size(sessions_withpairs,1);
for isessionpair = 1:nsessionpairs
    nsessionpairs_ana = nsessionpairs_ana + 1;
    %
    for ipair = 1:2  
        session = sessions_withpairs{isessionpair,ipair};
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
        rho_all(nsessionpairs_ana,ipair)     = relvalue;
        rho_OFF_all(nsessionpairs_ana,ipair) = relvalue_OFF;
        rho_ON_all(nsessionpairs_ana,ipair)  = relvalue_ON;
        rangeA_all(nsessionpairs_ana,ipair)  = rangeA;
        rangeB_all(nsessionpairs_ana,ipair)  = rangeB;     
        %
        sessionpairnames_ana(nsessionpairs_ana,ipair) = {session};
        
        % %
        % plot the results - PART1: two example
        % %
        if ismember(session,examsess_names)
            iexample = iexample + 1;
            figure(1)
            % Figure S2AB
            subplot(2,nexamples,iexample)
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
        end
        
    end % for ipair
end % for isessionpair

% %
% plot the results - PART2
% %
% Figure 4CD
subplot(2,nexamples,iexample+1)
deltarange = rangeA_all .* rho_all - rangeB_all;
deltarho   = rho_ON_all - rho_OFF_all;
ind_exam = ismember(sessionpairnames_ana, examsess_names);

%
% % only analyze sessions pairs that have opposite delta_ranges
% ind_tgtpairs = ismember(sign(deltarange(:,1)).*sign(deltarange(:,2)),-1,'rows');
% % analyze all session pairs
ind_tgtpairs = logical(ones(size(deltarange,1),1));
%
deltarange = deltarange(ind_tgtpairs,:);
deltarho   = deltarho(ind_tgtpairs,:);
ind_exam   = ind_exam(ind_tgtpairs,:);
%

% plot each pair
ind_exam = find(ind_exam(:,1)==1);
npairs = size(deltarange,1);
slopes_all = [];
for ipair = 1:npairs
    XXX = deltarange(ipair,:);
    YYY = deltarho(ipair,:);
    hold on
    plot(XXX,YYY,'ko-','MarkerSize',8);
    if ipair == ind_exam
        plot(XXX,YYY,'ko','MarkerSize',8,'MarkerFaceColor',[.4 .6 .4]);
        plot(XXX,YYY,'-','Color',[.4 .6 .4]);
    end
    slopes_all(ipair,1) = (YYY(2)-YYY(1))/(XXX(2)-XXX(1));
end
%
XXX = [ deltarange(:,2); deltarange(:,1) ];
YYY = [ deltarho(:,2); deltarho(:,1) ];
aa = polyfit(XXX,YYY,1);
[RR_Spe_all,pp_Spe_all] = corr(XXX,YYY,'Type','Spearman');
[RR_Pea_all,pp_Pea_all] = corr(XXX,YYY,'Type','Pearson');
XX  = [-8,11];
YY  = aa(1)*XX+aa(2);
YY2 = [-0.25,0.4];
plot(XX,YY,'Color',[.7 .7 .7],'LineWidth',3);
plot(XX,[0 0],'k--','LineWidth',1);
plot([0 0],YY2,'k--','LineWidth',1);
p_slope_signrank = signrank(slopes_all,0,'tail','both','method','exact');
[~,p_slope_ttest]= ttest(slopes_all,0,'tail','both');
text(6,-0.125,{['Pearson:'];...
              ['r = ',num2str(RR_Pea_all,'%.2f'),', p = ',num2str(pp_Pea_all,'%.3f')];...
              ['Spearman:'];...
              ['r = ',num2str(RR_Spe_all,'%.2f'),', p = ',num2str(pp_Spe_all,'%.3f')];...
              [''];['mean slope = ',num2str(mean(slopes_all),'%.3f')];...
              ['ttest: p = ',num2str(p_slope_ttest,'%.2f')];...
              ['Wilcoxon: p = ',num2str(p_slope_signrank,'%.3f')]},...
              'fontsize', 8, 'color', 'k');  
axis([XX,YY2]); axis square
box off
xlabel('\Delta VA - \Delta VB');
ylabel('\rhostimON - \rhostimOFF');
title({['both monkeys,N pairs = ',num2str(npairs)]});


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