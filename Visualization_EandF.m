% Creates the Visualizations of BiondiRighi (2018, JEIC) concerning the
% study of the effects additive accumulation process. Data for
% these visualizations can be created through the file Experiment_EandF.m


load('baseline_setup.mat')
figure(1)
load('02InterestSimple.mat');
plot(gini_avg,'b')
load('03InterestSimpleGamma.mat');
hold on
plot(gini_avg,'k');
title('Gini Coefficient','FontSize',19);
xlabel('Time')
ylabel('Gini Coefficient');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(1,'Exp_EF_Gini','fig')
print -depsc Exp_EF_Gini.eps

load('baseline_setup.mat')
figure(2)
load('02InterestSimple.mat');
semilogy(theil_avg,'b')
load('03InterestSimpleGamma.mat');
hold on
semilogy(theil_avg,'k');
title('Theil Coefficient','FontSize',19);
xlabel('Time')
ylabel('Theil Coefficient');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(2,'Exp_EF_Theil','fig')
print -depsc Exp_EF_Theil.eps

load('baseline_setup.mat')
figure(3)
load('02InterestSimple.mat');
plot(mav_Weighted_movement_W,'b')
load('03InterestSimpleGamma.mat');
hold on
plot(mav_Weighted_movement_W,'k');
title('Moving Average Weighted Movements Index (M_t)','FontSize',19);
xlabel('Time')
ylabel('Weighted Movements Index (M_t)');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(3,'Exp_EF_MavWeightedMovementW','fig')
print -depsc Exp_EF_MavWeightedMovementW.eps

load('baseline_setup.mat')
figure(4)
load('02InterestSimple.mat');
semilogy(Varricchezzeabs,'b')
load('03InterestSimpleGamma.mat');
hold on
semilogy(Varricchezzeabs,'k');
title('Absolute Mean Wealth Change Index (| V_t |)','FontSize',19);
xlabel('Time')
ylabel('Log(| V_t |)');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(4,'Exp_EF_Varricchezzeabs','fig')
print -depsc Exp_EF_Varricchezzeabs.eps

load('baseline_setup.mat')
figure(5)
load('02InterestSimple.mat');
loglog([2:1:tmax],Varricchezze(2:end),'b')
load('03InterestSimpleGamma.mat');
hold on
loglog([2:1:tmax],Varricchezze(2:end),'k');
title('Mean Wealth Change Index (V_t)','FontSize',19);
xlabel('Log(Time)')
ylabel('Log(V_t)');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(5,'Exp_EF_Varricchezze','fig')
print -depsc Exp_EF_Varricchezze.eps

load('baseline_setup.mat')
figure(6)
load('02InterestSimple.mat');
plot(share_WY,'b')
load('03InterestSimpleGamma.mat');
hold on
plot(share_WY,'k');
title('Share Wealth from Saved Income with reinvestment: Y_{i}^{t}','FontSize',19);
xlabel('Time')
ylabel('Share Wealth from Saved Income');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(6,'Exp_EF_share_WY','fig')
print -depsc Exp_EF_share_WY.eps

load('baseline_setup.mat')
figure(7)
subplot(1,2,1)
load('02InterestSimple.mat');
histfit(Final_wealth);
title('Simple Interest: Normal');
ylabel('Frequency');
subplot(1,2,2)
title('Simple Interest: Gamma','FontSize',19);
load('03InterestSimpleGamma.mat');
hold on
histfit(Final_wealth);
ylabel('Frequency');
saveas(7,'Exp_EF_FinalWealthDistro','fig')
print -depsc Exp_EF_FinalWealthDistro.eps

load('baseline_setup.mat')
figure(8)
load('02InterestSimple.mat');
[Num,Xvals] = hist(Final_wealth,100);
loglog(Xvals,Num/N,'b')
load('03InterestSimpleGamma.mat');
hold on
[Num,Xvals] = hist(Final_wealth,100);
loglog(Xvals,Num/N,'k')
title('Final Distribution of wealth','FontSize',19);
ylabel('Frequency');
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
saveas(8,'Exp_EF_FinalWealthDistro_alltogether','fig')
print -depsc Exp_EF_FinalWealthDistro_alltogether.eps


load('baseline_setup.mat')
figure(9)
load('02InterestSimple.mat');
loglog(B_sf,A_sf','b+')
%hold on
%plot(log(D_sf),Y_top,'k--')
load('03InterestSimpleGamma.mat');
hold on
loglog(B_sf,A_sf','k+')
%hold on
%plot(log(D_sf),Y_top,'k--')
xlabel('log(rank)')
ylabel('log(Wealth)')
title('log(wealth) over log(rank)','FontSize',19)
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
saveas(9,'Exp_EF_LogWealth_vs_LogRank','fig')
print -depsc Exp_EF_LogWealth_vs_LogRank.eps


%%%%
load('baseline_setup.mat')
figure(10)
subplot(1,2,1)
load('02InterestSimple.mat');
bar(wealthquartiles','histc'); 
ylabel('Share of Wealth');
xlabel('Quartiles')
ylim([0 1])
title('Simple Interest: Normal','FontSize',14);
subplot(1,2,2)
title('Simple Interest: Gamma','FontSize',14);
load('03InterestSimpleGamma.mat');
hold on
bar(wealthquartiles','histc'); 
ylabel('Share of Wealth');
xlabel('Quartiles')
ylim([0 1])
legend('t=1','t=10','t=100','t=1000','t=tmax',0);
saveas(10,'Exp_EF_WealthQuartiles','fig')
print -depsc Exp_EF_WealthQuartiles.eps
%%%%

load('baseline_setup.mat')
figure(11)
subplot(1,2,1)
load('02InterestSimple.mat');
bar(wealthdecile','histc'); 
ylabel('Share of Wealth');
xlabel('Deciles')
title('Simple Interest: Normal','FontSize',14);
xlim([0 11]);
ylim([0 1])
subplot(1,2,2)
title('Simple Interest: Gamma','FontSize',14);
load('03InterestSimpleGamma.mat');
hold on
bar(wealthdecile','histc'); 
ylabel('Share of Wealth');
xlabel('Deciles')
xlim([0 11])
ylim([0 1])
legend('t=1','t=10','t=100','t=1000','t=tmax',0);
saveas(11,'Exp_EF_WealthDeciles','fig')
print -depsc Exp_EF_WealthDeciles.eps



load('baseline_setup.mat')
figure(12)
load('02InterestSimple.mat');
plot(Proportion_total_wealth,'b')
load('03InterestSimpleGamma.mat');
hold on
plot(Proportion_total_wealth,'k');
title('Proportion of total wealth owned by top 1%');
xlabel('Time')
ylabel('Proportion of Wealth')
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(12,'Exp_EF_proportion_wealth_top1perc','fig')
print -depsc Exp_EF_proportion_wealth_top1perc.eps

load('baseline_setup.mat')
figure(13)
load('02InterestSimple.mat');
plot(Proportion_relative_wealth,'b')
load('03InterestSimpleGamma.mat');
hold on
plot(Proportion_relative_wealth,'k');
title('Proportion of wealth owned by top 1% as proportion of the bottom 50%');
xlabel('Time')
ylabel('Proportion of Wealth')
legend('Simple Interest: Normal Distribution','Simple Interest: Gamma Distribution',0);
xlim([0 tmax])
saveas(13,'Exp_EF_rel_proportion_top1bottom50','fig')
print -depsc Exp_EF_rel_proportion_top1bottom50.eps

