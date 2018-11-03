% Creates the Visualizations of BiondiRighi (2018, JEIC) concerning the
% effect of  additive wealth accumulation process
% Data for these visualizations can be created through the file Experiment_EandF.m

close all
clear all
clc

times=0:250:5000; times(1)=1;

spacing=(3.71-1)/(21-10);
timesval=1;
timeslog=1:9;
for i=10:21
    timeslog(i)=round(10^timesval);
    timesval=timesval+spacing;
end
timeslog(end)=5000;    

load('baseline_setup.mat')
figure(1)
load('02InterestSimple.mat'); 
errorbar(times,mean(gini_avg(times,:),2),std(gini_avg(times,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
errorbar(times,mean(gini_avg(times,:),2),std(gini_avg(times,:),0,2),'k');
title('Gini Coefficient','FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('Gini Coefficient', 'FontSize', 19);
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
saveas(1,'Exp_EF_Gini','fig')
print -depsc Exp_EF_Gini.eps


load('baseline_setup.mat')
figure(2)
load('02InterestSimple.mat'); 
errorbar(times,mean(theil_avg(times,:),2),std(theil_avg(times,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
errorbar(times,mean(theil_avg(times,:),2),std(theil_avg(times,:),0,2),'k');
title('Theil Coefficient','FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('Theil Coefficient', 'FontSize', 19);
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
saveas(2,'Exp_EF_Theil','fig')
print -depsc Exp_EF_Theil.eps

load('baseline_setup.mat')
figure(3)
load('02InterestSimple.mat'); 
timesma=times; timesma(end)=length(mav_Weighted_movement_W);
errorbar(timesma,mean(mav_Weighted_movement_W(timesma,:),2),std(mav_Weighted_movement_W(timesma,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
errorbar(timesma,mean(mav_Weighted_movement_W(timesma,:),2),std(mav_Weighted_movement_W(timesma,:),0,2),'k');
title('Moving Average Weighted Movements Index (M_t)','FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('Weighted Movements Index (M_t)', 'FontSize', 19);
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
ylim([0.16 0.34])
saveas(3,'Exp_EF_MavWeightedMovementW','fig')
print -depsc Exp_EF_MavWeightedMovementW.eps



load('baseline_setup.mat')
figure1=figure(4);
axes1 = axes('Parent',figure1);
load('02InterestSimple.mat'); 
for i=10:tmax; Varricchezzeabs_MA(i-9,:)=sum(Varricchezzeabs(i-9:i,:),1)./10; end;
timesma=timeslog; timesma(end)=length(Varricchezzeabs_MA);
errorbar(timesma,mean(Varricchezzeabs_MA(timesma,:),2),std(Varricchezzeabs_MA(timesma,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
for i=10:tmax; Varricchezzeabs_MA(i-9,:)=sum(Varricchezzeabs(i-9:i,:),1)./10; end;
errorbar(timesma,mean(Varricchezzeabs_MA(timesma,:),2),std(Varricchezzeabs_MA(timesma,:),0,2),'k');
title('Absolute Mean Wealth Change Index (|V_t|)','FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('| V_t |', 'FontSize', 19);
set(axes1,'XMinorTick','on','XScale','log');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
saveas(4,'Exp_EF_Varricchezzeabs','fig')
print -depsc Exp_EF_Varricchezzeabs.eps

load('baseline_setup.mat')
figure1=figure(5);
axes1 = axes('Parent',figure1);
load('02InterestSimple.mat');  %semilogx
timesma(1)=2;
errorbar(timesma,mean(Varricchezze(timesma,:),2),std(Varricchezze(timesma,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
errorbar(timesma,mean(Varricchezze(timesma,:),2),std(Varricchezze(timesma,:),0,2),'k');
title('Mean Wealth Change Index (V_t)','FontSize',19);
xlabel('Log(Time)', 'FontSize', 19)
ylabel('V_t', 'FontSize', 19);
set(axes1,'XMinorTick','on','XScale','log');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
ylim([0 0.3])
saveas(5,'Exp_EF_Varricchezze','fig')
print -depsc Exp_EF_Varricchezze.eps


% % % load('baseline_setup.mat')
% % % figure(7)
% % % subplot(3,2,1)
% % % load('02InterestSimple.mat'); 
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % histfit(Final_wealth);
% % % title('Normal \mu=0.05, \sigma=0.05','FontSize',14);
% % % ylabel('Frequency');
% % % subplot(3,2,4)
% % % title('Normal \mu=0.025, \sigma=0.05','FontSize',14);
% % % load('02InterestSimple.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % histfit(Final_wealth);
% % % ylabel('Frequency');
% % % subplot(3,2,6)
% % % title('Normal \mu=0.075, \sigma=0.05','FontSize',14);
% % % load('03InterestSimpleGamma.mat'); 
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % hist(Final_wealth,71);
% % % ylabel('Frequency');
% % % subplot(3,2,3)
% % % title('Normal \mu=0.05, \sigma=0.025','FontSize',14);
% % % load('04mu005sigma0025.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % histfit(Final_wealth);
% % % ylabel('Frequency');
% % % subplot(3,2,5)
% % % title('Normal \mu=0.05, \sigma=0.075','FontSize',14);
% % % load('05mu005sigma0075.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % histfit(Final_wealth);
% % % ylabel('Frequency');
% % % subplot(3,2,2)
% % % title('Gamma a=0.25, b=0.2','FontSize',14);
% % % load('06a02b025Gamma.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % histfit(Final_wealth);
% % % ylabel('Frequency');
% % % saveas(7,'Exp_EF_FinalWealthDistro','fig')
% % % print -depsc Exp_EF_FinalWealthDistro.eps
% % 
% % 
% % % load('baseline_setup.mat')
% % % figure(8)
% % % load('02InterestSimple.mat'); 
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'b')
% % % load('02InterestSimple.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'k')
% % % load('03InterestSimpleGamma.mat'); 
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'g')
% % % load('04mu005sigma0025.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'r')
% % % load('05mu005sigma0075.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'m')
% % % load('06a02b025Gamma.mat');
% % % hold on
% % % Final_wealth_v=[]; for i=1:Niter; Final_wealth_v=[Final_wealth_v Final_wealth(i,:)]; end; clear Final_wealth; Final_wealth=Final_wealth_v;
% % % [Num,Xvals] = hist(Final_wealth,100);
% % % loglog(Xvals,Num/N,'Color',[255./255 196./255 0./255])
% % % title('Final Distribution of wealth','FontSize',19);
% % % ylabel('Frequency');
% % % legend('Normal \mu=0.05, \sigma=0.05','Normal \mu=0.025, \sigma=0.05','Normal \mu=0.075, \sigma=0.05',...
% % %     'Normal \mu=0.05, \sigma=0.025','Normal \mu=0.05, \sigma=0.075','Gamma a=0.25, b=0.2');
% % % saveas(8,'Exp_EF_FinalWealthDistro_alltogether','fig')
% % % print -depsc Exp_EF_FinalWealthDistro_alltogether.eps



%%%%
load('baseline_setup.mat')
figure(10)
subplot(1,2,1)
load('02InterestSimple.mat'); 
bar(mean(wealthquartiles,3)','histc'); 
ylabel('Share of Wealth', 'FontSize', 12);
xlabel('Quartiles', 'FontSize', 12)
title('Normal \mu=0.05, \sigma=0.05','FontSize',14);
subplot(1,2,2)
title('Normal \mu=0.025, \sigma=0.05','FontSize',14);
load('03InterestSimpleGamma.mat'); 
hold on
bar(mean(wealthquartiles,3)','histc'); 
ylabel('Share of Wealth', 'FontSize', 12);
xlabel('Quartiles', 'FontSize', 12)
legend('t=1','t=10','t=100','t=1000','t=tmax','Location','best');
saveas(10,'Exp_EF_WealthQuartiles','fig')
print -depsc Exp_EF_WealthQuartiles.eps

load('baseline_setup.mat')
figure(11)
subplot(1,2,1)
load('02InterestSimple.mat'); 
bar(mean(wealthdecile,3)','histc'); 
ylabel('Share of Wealth', 'FontSize', 12);
xlabel('Deciles', 'FontSize', 12)
title('Normal \mu=0.05, \sigma=0.05','FontSize',14);
xlim([0 11]);
subplot(1,2,2)
title('Normal \mu=0.025, \sigma=0.05','FontSize',14);
load('03InterestSimpleGamma.mat'); 
hold on
bar(mean(wealthdecile,3)','histc'); 
ylabel('Share of Wealth', 'FontSize', 12);
xlabel('Deciles', 'FontSize', 12)
xlim([0 11])
legend('t=1','t=10','t=100','t=1000','t=tmax','Location','best');
saveas(11,'Exp_EF_WealthDeciles','fig')
print -depsc Exp_EF_WealthDeciles.eps



figure1=figure(18)
axes1 = axes('Parent',figure1);
load('02InterestSimple.mat'); 
errorbar(times,mean(Proportion_total_wealth(times,:),2),std(Proportion_total_wealth(times,:),0,2),'b')
load('03InterestSimpleGamma.mat'); 
hold on
errorbar(times,mean(Proportion_total_wealth(times,:),2),std(Proportion_total_wealth(times,:),0,2),'k');
title({'Proportion of total wealth owned by top 1%';' '},'FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('Proportion of Wealth', 'FontSize', 19)
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
xlim([0 tmax])
ylim([0 1])
set(axes1,'XMinorTick','on','YScale','log');

saveas(18,'Exp_EF_proportion_wealth_top1perc','fig')
print -depsc Exp_EF_proportion_wealth_top1perc.eps



figure1=figure(19)
axes1 = axes('Parent',figure1);
load('02InterestSimple.mat'); 
stdval=std(Proportion_relative_wealth(timeslog,:),0,2);
lb=max(mean(Proportion_relative_wealth(timeslog,:),2),0.0001)-stdval;
hb=mean(Proportion_relative_wealth(timeslog,:),2)+stdval;
errorbar(timeslog,mean(Proportion_relative_wealth(timeslog,:),2),lb,hb,'b')
load('03InterestSimpleGamma.mat'); 
hold on
stdval=std(Proportion_relative_wealth(timeslog,:),0,2);
lb=max(mean(Proportion_relative_wealth(timeslog,:),2),0.0001)-stdval;
hb=mean(Proportion_relative_wealth(timeslog,:),2)+stdval;
errorbar(timeslog,mean(Proportion_relative_wealth(timeslog,:),2),lb,hb,'k');
set(gca,'yscale','log')
set(axes1,'XMinorTick','on','XScale','log','YScale','log');
ax = get(figure1,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution','Location','best');
title({'Proportion of wealth owned by top 1%'; 'as proportion of the bottom 50%'},'FontSize',19);
xlabel('Time', 'FontSize', 19)
ylabel('Proportion of Wealth', 'FontSize', 19)
xlim([0 tmax])
ylim([10^-5 10^0])
saveas(19,'Exp_EF_rel_proportion_top1bottom50','fig')
print -depsc Exp_EF_rel_proportion_top1bottom50.eps