% FigS5_plot_generator.m
% This script uses complete_dataset.xlsx to generate figure S5
% Figure S5 - no steepness change of Exp2

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
doremoveoutliers = 1;  % remove outliers
limout=3;              % outlier threshold

h(1) = figure;
set(gcf,'position',[510 265 550 550], 'PaperPositionMode','auto')
axes('position',[.1 .97 .2 .05]);text(0,0,['Figure S5. No steepness change from stimulation'],'fontsize',10);axis off   
                                                                                        
%                                                                                    
nsessions = length(allsessions);
%
% initialization
nsessions_ana     = 0;
steepness_ON_all  = [];
steepness_OFF_all = [];
rho_all          = [];
    
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

    % calculate steepness
    beta0 = psyphy.beta(1);
    beta1 = psyphy.beta(2);
    beta2 = psyphy.beta(3);
    beta3 = psyphy.beta(4);    
    %
    steepness_ON_all(nsessions_ana,1)  = beta1;
    steepness_OFF_all(nsessions_ana,1) = beta3;

end
    
% %
% plot the results
% %
% Figure S5
figure(1)
subplot(1,1,1)

% % % remove outliers
if doremoveoutliers
    K = limout;
    upper_off = mean(steepness_OFF_all)+K*std(steepness_OFF_all);
    lower_off = mean(steepness_OFF_all)-K*std(steepness_OFF_all);
    ind_off = [steepness_OFF_all<upper_off & steepness_OFF_all>lower_off ];
    upper_on = mean(steepness_ON_all)+K*std(steepness_ON_all);
    lower_on = mean(steepness_ON_all)-K*std(steepness_ON_all);
    ind_on = [steepness_ON_all<upper_on & steepness_ON_all>lower_on ];
    ind_steepness = ind_off & ind_on; 
else
    ind_steepness = logical(ones(size(steepness_OFF_all)));
end
XXX = steepness_OFF_all(ind_steepness);
YYY = steepness_ON_all(ind_steepness);
plot(XXX,YYY,'ro','MarkerSize',10);hold on
[~,p_ttest]=ttest(XXX,YYY);
[p_wil,~] = signrank(XXX,YYY);
XX = [0, 14]; 
YY = XX;
plot(XX,YY,'k--','LineWidth',1); 
Sigma_ell = cov(XXX,YYY);
mu_ell(1) = mean(XXX);
mu_ell(2) = mean(YYY);       
plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
text(1,12,{['ttest:'];...
              ['p = ',num2str(p_ttest,'%.3f')];...
              ['Wilcoxon:'];...
              ['p = ',num2str(p_wil,'%.3f')];...
              ['session # = ',num2str(length(XXX))]}...
              , 'fontsize', 10, 'color', 'k');  
axis([XX,YY]); axis square
box off
xlabel('\eta stim OFF');
ylabel('\eta stim ON');




% % %
% function 1
% % %
function [psyphy, table01] = probitLink_regression(QtyA, QtyB, ChosenID, stim)
%
% This function is to do logistic regression with "probit" link function
% This regression assumes stimON and stimOFF have different steepness,
% therefore X = (b0 + b1 * log(#B/#A))*delta_stimON + (b2 + b3 * log(#B/#A))*delta_stimOFF
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
[~ , ~, stats] = glmfit([double(iON(jj)), logB_A(jj,1).*double(iON(jj)), double(iOF(jj)), logB_A(jj,1).*double(iOF(jj))], ChosenID(jj,1), 'binomial', 'link','probit', 'constant','off');
psyphy.vars	= {'const_stimOn', 'logit(#B/#A)_stimOn', 'const_stimOff', 'logit(#B/#A)_stimOff'};
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
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2),'Color',[.6 .6 .6],'LineWidth',3);
end