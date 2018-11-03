% Creates the Visualizations of BiondiRighi (2018, JEIC) concerning the
% study of different modalities of taxation and redistribution. Data for
% these visualizations can be created through the file Experiment_M.m
% Note that more visualizations on these simulations are available running
% the file Visualization_M_1.m 

load('baseline_setup.mat')
figure(25)
load('02TaxesProportionalLiberalist.mat');
hold on
plot(ratio_top1T_bottom10W,'k-o');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(ratio_top1T_bottom10W,'g-o');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(ratio_top1T_bottom10W,'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(ratio_top1T_bottom10W,'m');
title('Ratio between Tax Levy on Top 1% and Wealth of bottom 10%','FontSize',19);
xlabel('Time')
ylabel('Ratio');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(25,'Exp_M_ratioTop1T_bottom10W','fig')
print -depsc Exp_M_ratioTop1T_bottom10W.eps


load('baseline_setup.mat')
figure(26)
load('02TaxesProportionalLiberalist.mat');
hold on
plot(ratio_top10T_bottom50W,'k-o');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(ratio_top10T_bottom50W,'g-o');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(ratio_top10T_bottom50W,'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(ratio_top10T_bottom50W,'m');
title('Ratio between Tax Levy on Top 10% and Wealth of bottom 50%','FontSize',19);
xlabel('Time')
ylabel('Ratio');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(26,'Exp_M_ratioTop10T_bottom50W','fig')
print -depsc Exp_M_ratioTop10T_bottom50W.eps


load('baseline_setup.mat')
figure(27)
load('02TaxesProportionalLiberalist.mat');
hold on
plot(ratio_top50T_bottom50W,'k-o');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(ratio_top50T_bottom50W,'g-o');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(ratio_top50T_bottom50W,'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(ratio_top50T_bottom50W,'m');
title('Ratio between Tax Levy on Top 50% and Wealth of bottom 50%','FontSize',19);
xlabel('Time')
ylabel('Ratio');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(27,'Exp_M_ratioTop50T_bottom50W','fig')
print -depsc Exp_M_ratioTop50T_bottom50W.eps


load('baseline_setup.mat')
figure(12)
%load('01Baseline.mat');
%plot(MeanTaxRate(2:end),'b')
load('02TaxesProportionalLiberalist.mat');
hold on
plot(MeanTaxRate(2:end),'k');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(MeanTaxRate(2:end),'g');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(MeanTaxRate(2:end),'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(MeanTaxRate(2:end),'m');
title('Mean Tax Rate','FontSize',19);
xlabel('Time')
ylabel('Effective Mean Tax Rate');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(12,'Exp_M_MeanTaxRate','fig')
print -depsc Exp_M_MeanTaxRate.eps


load('baseline_setup.mat')
figure(13)
%load('01Baseline.mat');
%plot(MedianTaxRate(2:end),'b')
load('02TaxesProportionalLiberalist.mat');
hold on
plot(MedianTaxRate(2:end),'k');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(MedianTaxRate(2:end),'g');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(MedianTaxRate(2:end),'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(MedianTaxRate(2:end),'m');
title('Median Tax Rate','FontSize',19);
xlabel('Time')
ylabel('Effective Median Tax Rate');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(13,'Exp_M_MedianTaxRate','fig')
print -depsc Exp_M_MedianTaxRate.eps

load('baseline_setup.mat')
figure(14)
%load('01Baseline.mat');
%plot(MeanRedistributionRate(2:end),'b')
load('02TaxesProportionalLiberalist.mat');
hold on
plot(MeanRedistributionRate(2:end),'k');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(MeanRedistributionRate(2:end),'g');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(MeanRedistributionRate(2:end),'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(MeanRedistributionRate(2:end),'m');
title('Mean Redistribution Rate','FontSize',19);
xlabel('Time')
ylabel('Effective Mean Redistribution Rate');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(14,'Exp_M_MeanRedistributionRate','fig')
print -depsc Exp_M_MeanRedistributionRate.eps

load('baseline_setup.mat')
figure(15)
%load('01Baseline.mat');
%plot(MedianRedistributionRate(2:end),'b')
load('02TaxesProportionalLiberalist.mat');
hold on
plot(MedianRedistributionRate(2:end),'k');
load('03TaxesProgressiveLiberalist.mat');
hold on
plot(MedianRedistributionRate(2:end),'g');
load('04TaxesProportionaKeynesian.mat');
hold on
plot(MedianRedistributionRate(2:end),'r');
load('05TaxesProgressiveKeynesian.mat');
hold on
plot(MedianRedistributionRate(2:end),'m');
title('Median Redistribution Rate','FontSize',19);
xlabel('Time')
ylabel('Effective Median Redistribution Rate');
legend('Proportional & PS','Progressive & PS',...
    'Proportional & Welfare','Progressive & Welfare',0);
xlim([0 tmax])
saveas(15,'Exp_M_MedianRedistributionRate','fig')
print -depsc Exp_M_MedianRedistributionRate.eps
