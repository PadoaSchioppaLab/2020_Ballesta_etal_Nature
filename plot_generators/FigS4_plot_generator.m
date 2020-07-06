% FigS4_plot_generator.m
% This script uses complete_dataset.xlsx to generate figure S4
% Figure S4 - Exp1, analysis of range-dependent bias

% Authors: Sebastien Ballesta & Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

% load data and info files
[~,~,Exp1_info] = xlsread('complete_dataset.xlsx', 2); % sheet 2
[~,~,Exp1_data] = xlsread('complete_dataset.xlsx', 3); % sheet 3
info=cell2table(Exp1_info(2:end,:),'VariableNames',Exp1_info(1,:));
alldata=cell2table(Exp1_data(2:end,:),'VariableNames',Exp1_data(1,:));

% DO FITs
[info] = probitLink_regression(alldata, info);

% %
AR = info.relative_value_ON-info.relative_value_OFF; % relative value difference
mrho = info.relative_value;                          % mean relative value
RR = info.rangeA-info.rangeB;                        % value range difference
isplit = info.isplit;                                % index for good sessions fulfilling incomplete separation
isat = info.isat;                                    % index for good sessions fulfilling saturation
stimcurr = info.StimCurrent;                         % stimulation current level
stimlat = info.StimLateral;                          % uni or bilateral micro-stimulation

%
doremoveoutliers = 1;  % remove outliers based on RR and AR
limout=3;              % outlier threshold


% current levels
limlow=25;
limhigh=50;
limhigh2=100;

%
for k = [1 2 3] % different current levels
    
    mask=[];
    switch k    % define mask1 to index sessions with target current level
        case 1
            mask1= stimcurr>=limlow & stimcurr<limhigh ;
            tit='25\muA';
        case 2
            mask1= stimcurr>=limhigh & stimcurr<limhigh2;
            tit='50\muA';
        case 3
            mask1= stimcurr>=limhigh2 ;
            tit='\geq100\muA';
    end
    mask = logical(mask1);

    % % 
    % REMOVE OUTLIERS
    if doremoveoutliers
        K = limout;
        OUT_RR = (RR>(mean(RR)+std(RR)*K)) | (RR<(mean(RR)-std(RR)*K));  
        OUT_AR = (AR>(mean(AR)+std(AR)*K)) | (AR<(mean(AR)-std(AR)*K));  
        %
        OUT = OUT_RR | OUT_AR;
        %
        mask=logical(~OUT.*mask);
    else
        OUT=~mask;
    end
    disp([' Removed outliers = ' num2str(sum(mask1(OUT)))])
    mask=logical(mask);

    % %
    % plot
    figure(1)
    set(gcf,'position',[110 65 1550 550], 'PaperPositionMode','auto');
    subplot(1,3,k)
    hold on;
    XX=RR;
    YY=AR;
    plot(XX(mask),YY(mask),'ko','markersize',6); 
    plot([-8 16], [0 0],'k--','LineWidth',1);
    plot([0 0], [-0.5 1], 'k--','LineWidth',1);
    [R,P]=corrcoef(XX(mask),YY(mask));
    [Rsp,Psp]=corr(XX(mask),YY(mask),'Type','Spearman');
    xlabel('\DeltaVA - VB' );
    ylabel('\rhostim ON - \rhostim OFF' ) ;
    title(tit)
    axis([-8 16 -0.5 1])
    axis square
    box off
    %linear fit
    warning off
    xx=XX(mask); yy=YY(mask);
    g = fittype('m*x+q','coeff',{'m','q'});
    [linfit,goodness] = fit(xx,yy,g);
    intercept = linfit.q;
    slope= linfit.m;
    interval = confint(linfit,0.95);
    XXplot = [min(xx),max(xx)];
    YYplot = slope.*XXplot + intercept;
    plot(XXplot,YYplot,'LineWidth',3,'Color',[0.7 0.7 0.7]);

    text(5,-0.25,{ ['Pearson: r = ', num2str(round(R(2),3)), ', p = ', num2str(round(P(2),3))]...
                   ['Spearman: r = ', num2str(round(Rsp,3)), ', p = ', num2str(round(Psp,3))]...
                   ['n sessions= ' num2str(sum(mask))] ['Removed sess= ' num2str(sum(mask1(~mask)))]}...
                    ,'Color', 'k');

end



% % %
% function 1
% % %
function [info] = probitLink_regression(alldata, info)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A)
% regression is done on stimOFF and stimON trials seperately
%
nsessions = size(info,1);
%
for isession = 1:nsessions
    session = info.session{isession};
    % load raw data from alldata
    ind_trials = ismember(alldata.Session, ['ST',session]);
    QtyA = alldata.QtyA(ind_trials,:);
    QtyB = alldata.QtyB(ind_trials,:);
    stim = alldata.Stim(ind_trials,:);
    chosenID = alldata.ChosenID(ind_trials,:);
    order = -alldata.Order(ind_trials,:);
    %
    % remove forced choices
    ii = ~prod([QtyA,QtyB],2)==0;
    %
    iOF = stim==-1;
    iON = stim==1;
    
    % logistic regression - based on single trials
    % for stimOFF
    logB_A = log(abs(QtyB)./abs(QtyA));
    jj = ii & iOF;
    mdl = fitglm([logB_A(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
    betaOFF = mdl.Coefficients.Estimate;
    rhoOFF = exp(-betaOFF(1)/betaOFF(2));
    steepnessOFF = betaOFF(2);
    % for stimON
    jj = ii & iON;
    mdl = fitglm([logB_A(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
    betaON = mdl.Coefficients.Estimate;
    rhoON = exp(-betaON(1)/betaON(2));
    steepnessON = betaON(2); 
               
    % save data
    relative_value(isession,:) = (rhoOFF + rhoON)/2;
    relative_value_OFF(isession,:) = rhoOFF;
    relative_value_ON(isession,:) = rhoON;
    steepness_OFF(isession,:) = steepnessOFF;
    steepness_ON(isession,:)  = steepnessON;
    rangeA_all(isession,:) = max(QtyA)*(rhoOFF + rhoON)/2;
    rangeB_all(isession,:) = max(QtyB);
    
    
    %%%% FILTER SESSIONS WITH NO SPLIT OR NOT SATURATED OFFERTYPES
    %%%%% CREATE QUALITY FILTER % 2 SPLITS TRIALS + GOOD SATURATION
    %%%%% (>90%)
    nlow_split = .1; nhigh_split = .9;
    nlow_satur = .1; nhigh_satur = .9;
    ordermerge = zeros(size(order)); % AB BA trials are merged
    jj = ii & iOF;
    [tableOFF, binosizeOFF] = gettable(QtyA(jj,1),QtyB(jj,1),ordermerge(jj,1),chosenID(jj,1));
    jj = ii & iON;
    [tableON,  binosizeON ] = gettable(QtyA(jj,1),QtyB(jj,1),ordermerge(jj,1),chosenID(jj,1));    
    
    if 1 % FILTER ON ALL TRIALS
        % % stim OFF
        % % based on BAratio
        % uni_logBA = unique(tableOFF.offer);
        % nuni = size(uni_logBA,1);
        % BpercOFF = [];
        % for iuni = 1:nuni
        %     ind = ismember(tableOFF.offer,uni_logBA(iuni),'rows');
        %     BpercOFF(iuni,1) = sum(tableOFF.choice(ind))./sum(binosizeOFF(ind));
        % end
        % based on offer type
        BpercOFF = tableOFF.choice./binosizeOFF;
        t1OFF = sum(BpercOFF>nlow_split & BpercOFF<nhigh_split); % SPLIT TRIALS
        t2OFF = sum(BpercOFF<=nlow_satur)>=1 | sum(BpercOFF>=nhigh_satur)>=1; % SATURATED TRIALS
        % t2OFF = sum(BpercOFF<=nlow_satur)>=1; % SATURATED TRIALS
        % % stim ON
        % % based on BAratio
        % uni_logBA = unique(tableON.offer);
        % nuni = size(uni_logBA,1);
        % BpercON = [];
        % for iuni = 1:nuni
        %     ind = ismember(tableON.offer,uni_logBA(iuni),'rows');
        %     BpercON(iuni,1) = sum(tableON.choice(ind))./sum(binosizeON(ind));
        % end
        % based on offer type
        BpercON = tableON.choice./binosizeON;
        t1ON = sum(BpercON>nlow_split & BpercON<nhigh_split); % SPLIT TRIALS
        t2ON = sum(BpercON<=nlow_satur)>=1 | sum(BpercON>=nhigh_satur)>=1; % SATURATED TRIALS
        % t2ON = sum(BpercON<=nlow_satur)>=1; % SATURATED TRIALS
        %
        isplit(isession,1) = t1OFF>=1 & t1ON>=1 ;
        isat(isession,1) = t2OFF>=1 & t2ON>=1 ;        
    end
    
end % for isession

% save data to info
info.relative_value_OFF = relative_value_OFF;
info.relative_value_ON = relative_value_ON;
info.relative_value = relative_value;
info.steepnessOFF = steepness_OFF;
info.steepnessON = steepness_ON;
info.rangeA = rangeA_all;
info.rangeB = rangeB_all;
info.isat=isat;
info.isplit=isplit;
end


% % %
% function 2
% % %
function [table02,binosize] = gettable(QtyA,QtyB,order,chosenID)

alltrials = [QtyA,QtyB,order];
[trialtypes,~,~] = unique(alltrials,'rows');
ntrialtypes = size(trialtypes,1);

table01 = [];
for itrialtype = 1:ntrialtypes
    trialtype = trialtypes(itrialtype,:);
    table01(itrialtype,1:3) = trialtype;
    %
    ind_trialtypes = ismember(alltrials,trialtype,'rows');
    ntrials = sum(ind_trialtypes);
    nchBtrials = sum(chosenID(ind_trialtypes,1));
    %
    table01(itrialtype,4) = nchBtrials;
    table01(itrialtype,5) = ntrials;
end

logB_A = log(table01(:,2)./table01(:,1));
order2 = table01(:,3);
choice2 = table01(:,4);
binosize = table01(:,5);

table02=table(logB_A,order2,choice2,'VariableNames',{'offer','order','choice'});

end


