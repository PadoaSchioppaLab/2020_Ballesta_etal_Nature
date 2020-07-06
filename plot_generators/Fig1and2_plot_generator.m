% Fig1and2_plot_generator.m
% This script uses complete_dataset.xlsx to generate figures 1 and 2
% Figure 1 - Exp1, primary analysis, 100 uA
% Figure 2 - Exp1, primary analysis, all current levels

% Authors: Sebastien Ballesta & Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

% load data and info files
% if use complete_dataset.xlsx
[~,~,Exp1_info] = xlsread('complete_dataset.xlsx', 2); % sheet 2
[~,~,Exp1_data] = xlsread('complete_dataset.xlsx', 3); % sheet 3
info=cell2table(Exp1_info(2:end,:),'VariableNames',Exp1_info(1,:));
alldata=cell2table(Exp1_data(2:end,:),'VariableNames',Exp1_data(1,:));

examplesessions = {'170619c','170828d'};

% %
% % FIGURE 1
% %

% define analysis type and figure axis labels
ylab_f1={ '\geq100\muA offer1'; '\geq100\muA offer2';};
% current levels - figure 1 only analyzes >=100uA
limhigh2=100;

% %
f(1)=figure(1);
set(gcf,'position',[110 65 1250 1250], 'PaperPositionMode','auto')
hold on;
paroff=[  1 2  ]; % offer1 or offer2
parsup=[  limhigh2 limhigh2];  % use 100uA as lower boundary to choose target sessions
parinf=[  Inf Inf ];           % use inf as upper boundary to include all high current sessions
nplot=numel(paroff);
%
for k=1:3 % 3: three analysis target - relative value, steepness and order bias
	for n=1:numel(paroff)
		% mask contains index for target sessions
		mask=(info.StimOffer(:,1)==paroff(n) & info.StimCurrent>=parsup(n) & info.StimCurrent<parinf(n)) ;
		ind_exam = ismember(info.session(mask),examplesessions);
		% load data accordingly
		try
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 2
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 3
					tit={'Order bias'};
					data=[info.order_bias_OFF info.order_bias_ON ];
			end
		catch
			[info] = probitLink_regression(alldata, info);
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 2
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 3
					tit={'Order bias'};
					data=[info.order_bias_OFF info.order_bias_ON ];
			end
		end
		% find the limits for two axes for plotting
		fact=1;
		for i=1:size(data,1)
			limi(i,1)=min(min(data))*fact; limi(i,2)=max(max(data))*fact;
		end
		% only looks at target sessions
		data=data(mask,:);
		% plot
		ax{n,k}=subplot(nplot,3,k+(n-1)*3);
		hold on;
		% add ellipse
		eli=error_ellipse(cov(data),'mu',nanmean(data),'conf',.9);
		eli.LineWidth=2; eli.Color=[0.7 0.7 0.7];
		% add data points
		plot(data(:,1),data(:,2),'ko','MarkerSize',8);
		plot(data(ind_exam,1),data(ind_exam,2),'ko','MarkerSize',8,'MarkerFaceColor',[0.4 0.6 0.4]);
		plot(limi(k,:),limi(k,:),'k--','HandleVisibility','off','LineWidth',1);
		axis ([limi(k,:), limi(k,:)])
		axis square
		xlabel('Stim OFF');
		title(ylab_f1{n});
		if n ==1
			ylabel({tit{1} 'Stim ON'});
		else
			ylabel('Stim ON');
		end
	end
end
%
% statistical analysis
for k=1:3
	for n=1:numel(paroff)
		mask=(info.StimOffer(:,1)==paroff(n) & info.StimCurrent>=parsup(n) & info.StimCurrent<parinf(n));
		switch k
			case 1
				tit={'Rho'};
				data=[ info.relative_value_OFF info.relative_value_ON ];
			case 2
				tit={'Steepness'};
				data=[info.steepness_OFF info.steepness_ON ];
			case 3
				tit={'Order bias'};
				data=[info.order_bias_OFF info.order_bias_ON ];
		end
		% only looks at target sessions
		data=data(mask,:);
		%
		[~, p2(k,n)]=ttest(data(:,1),data(:,2));
		[ p(k,n) ] = signrank(data(:,1),data(:,2));
		pos=ax{n,k}.Position; AX{n,k}=axes('Parent',f,'Position',[pos(1) pos(2)-.1 pos(3) .1 ]);
		text(0,0.8,['  n sessions= ' num2str(sum(mask))],'Color','k');
		if p(k,n)>.05
			text(0,.6,['   Wilcoxon: p = ', num2str(p(k,n),'%.4f') ])
		else
			text(0,.6,['   Wilcoxon: p = ', num2str(p(k,n),'%.4f') ], 'FontWeight', 'bold')
		end
		if p2(k,n)>.05
			text(0,.4 ,['  ttest: p = ', num2str(p2(k,n),'%.4f') ])
		else
			text(0,.4, ['   ttest: p = ', num2str(p2(k,n),'%.4f') ], 'FontWeight', 'bold')
		end
		AX{n,k}.Visible ='off';
	end
end
suptitle('Figure 1')


% %
% % FIGURE 2
% %

% define analysis type and figure axis labels
ylab_f2={ 'NoStim'; '25\muA'; '25\muA Offer2'; '50\muA'; '50\muA Offer2'; '\geq100\muA'; '\geq100\muA offer2';};
% different current levels
limlow=25;
limhigh=50;
limhigh2=100;
%
paroff=[0 1 2 1 2 1 2]; % 0: non stimulation control; 1: offer1 stimulation; 2: offer2 stimulation
parsup=[0 limlow limlow limhigh limhigh limhigh2 limhigh2]; % lower boundary to choose target current level
parinf=[1 limhigh limhigh limhigh2 limhigh2 Inf Inf ]; % upper boundary to choose target current level
%
for k=1:3 % 3: three analysis target - relative value, steepness and order bias
	for n=1:numel(paroff)
		% mask contains index for target sessions
		mask=(info.StimOffer(:,1)==paroff(n) & info.StimCurrent>=parsup(n) & info.StimCurrent<parinf(n)) ;
		try
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 2
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 3
					tit={'Order bias'};
					data=[info.order_bias_OFF info.order_bias_ON ];
			end
		catch
			[info] = probitLink_regression(alldata, info);
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 2
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 3
					tit={'Order bias'};
					data=[info.order_bias_OFF info.order_bias_ON ];
			end
		end
		% only looks at target sessions
		data=data(mask,:);
		mdata(k,n) = mean(diff(data'));
		stdata(k,n) =std(diff(data'))./sqrt(sum(mask));
		[~,pval_ttest(k,n)] = ttest(diff(data'));
		pval_wil(k,n) = signrank(diff(data'));
		if isnan(mdata(k,n))
			keyboard
		end
		
	end
end
%
% plot
tit={'Relative value' , 'Steepness' , 'Order bias'};
f(2) = figure(2);
set(gcf,'position',[110 65 1450 850], 'PaperPositionMode','auto')
for k=1:3
	mdata_off1=mdata(k,paroff==1 | paroff==0);
	mdata_off2=mdata(k,paroff==2 | paroff==0);
	stddata_off1=stdata(k,paroff==1 | paroff==0);
	stddata_off2=stdata(k,paroff==2 | paroff==0);
	pttest_off1=pval_ttest(k,paroff==1 | paroff==0);
	pttest_off2=pval_ttest(k,paroff==2 | paroff==0);
	pwil_off1=pval_wil(k,paroff==1 | paroff==0);
	pwil_off2=pval_wil(k,paroff==2 | paroff==0);
	%
	subplot(1,3,k)
	plot([0 numel(parsup)/2+1], [0 0],'k--','HandleVisibility','off'); hold on
	h=errorbar([1:numel(parsup)/2+1],mdata_off1,stddata_off1,'o','Linewidth',1,'HandleVisibility','off');
	h.MarkerEdgeColor = [0.0 0.2 0.8];
	h.MarkerFaceColor = [0.0 0.2 0.8];
	h.MarkerSize = 10;
	h = errorbar([1:numel(parsup)/2+1],mdata_off2,stddata_off2,'o','Linewidth',1,'HandleVisibility','off');
	h.MarkerEdgeColor = [1.0 0.4 0.0];
	h.MarkerFaceColor = [1.0 0.4 0.0];
	h.MarkerSize = 10;
	pl1 = plot([1:numel(parsup)/2+1],mdata_off1,'-','Linewidth',1);
	pl1.Color = [0.0 0.2 0.8];
	pl2 = plot([1:numel(parsup)/2+1],mdata_off2,'-','Linewidth',1);
	pl2.Color = [1.0 0.4 0.0];
	title(tit{k});
	legend([ pl1 pl2],{'Offer1', 'Offer2'},'Location','northwest');
	set(gca,'Xtick',[1:numel(parsup)/2+1],'XtickLabel',{'NoStim', '25\muA', '50\muA', '\geq100\muA', ' '});
	xlabel('current level (\muA)')
	switch k
		case 1
			ylim([-.5 .5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<0.05;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.35,'k*');
			ind_signi2=pwil_off2<0.05;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.35,'k*');
			ind_signi1=pwil_off1<0.01;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.35,'m*');
			ind_signi2=pwil_off2<0.01;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.35,'m*');
			ind_signi1=pwil_off1<0.001;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.35,'r*');
			ind_signi2=pwil_off2<0.001;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.35,'r*');
			%
			ylabel('\rhostimON - \rhostimOFF');
		case 2
			ylim([-1.5 1.5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<0.05;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*1.1,'k*');
			ind_signi2=pwil_off2<0.05;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-1.1,'k*');
			ind_signi1=pwil_off1<0.01;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*1.1,'k*');
			ind_signi2=pwil_off2<0.01;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-1.1,'k*');
			ind_signi1=pwil_off1<0.001;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*1.1,'k*');
			ind_signi2=pwil_off2<0.001;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-1.1,'k*');
			%
			ylabel('\etastimON - \etastimOFF');
		case 3
			ylim([-.5 .5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<0.05;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.47,'k*');
			ind_signi2=pwil_off2<0.05;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.47,'k*');
			ind_signi1=pwil_off1<0.01;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.47,'k*');
			ind_signi2=pwil_off2<0.01;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.47,'k*');
			ind_signi1=pwil_off1<0.001;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.47,'k*');
			ind_signi2=pwil_off2<0.001;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.47,'k*');
			%
			ylabel('\epsilonstimON - \epsilonstimOFF');
	end
	xlim([0.5 4.5])
	axis square
	box off
end
suptitle('Figure 2')



% % %
% function 1
% % %
function [info] = probitLink_regression(alldata, info)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A) + b2 *(delta_AB - delta_BA)
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
	order = alldata.Order(ind_trials,:);
	%
	% remove forced choices
	ii = ~prod([QtyA,QtyB],2)==0;
	%
	iOF = stim==-1;
	iON = stim==1;
	%
	logB_A = log(abs(QtyB)./abs(QtyA));
	%
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),order(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1),order(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
	orderbiasOFF = betaOFF(3);
	
	% for stimON
	jj = ii & iON;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),order(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1),order(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	orderbiasON = betaON(3);
	
	
	% save data
	relative_value_OFF(isession,:) = rhoOFF;
	relative_value_ON(isession,:) = rhoON;
	steepness_OFF(isession,:) = steepnessOFF;
	steepness_ON(isession,:) = steepnessON;
	order_bias_OFF(isession,:) = orderbiasOFF;
	order_bias_ON(isession,:) = orderbiasON;
	
end % for isession

% save data to info
info.relative_value_OFF = relative_value_OFF;
info.relative_value_ON = relative_value_ON;
info.steepness_OFF = steepness_OFF;
info.steepness_ON = steepness_ON;
info.order_bias_OFF = order_bias_OFF;
info.order_bias_ON = order_bias_ON;
end


