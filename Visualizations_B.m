% Creates the Visualizations of BiondiRighi (2018, JEIC) concerning the
% exploration of different average and std in the returns of a
% multiplicative wealth accumulation process.
% Data for these visualizations can be created through the file Experiment_B.m

close all
clear all
clc

load('baseline_setup.mat')
load('Experiment_B')
figure(15)
imagesc(interesse_1,interesse_2,mean(WatZero_end,3)./N)
title('Proportion of Agents ending with W_{tmax}=0','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(15,'B_RuinedAgents','fig')
print -depsc B_RuinedAgents.eps



load('baseline_setup.mat')
load('Experiment_B')
figure(1)
imagesc(interesse_1,interesse_2,mean(P_logw_logr_top,3))
title('Slope of Power law (top 10%)','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(1,'Power_law_slope','fig')
print -depsc Power_law_slope.eps

load('baseline_setup.mat')
figure(2)
load('Experiment_B')
imagesc(interesse_1,interesse_2,mean(M_growth,3))
title('Average Aggregate Growth','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(2,'Mean_Growth','fig')
print -depsc Mean_Growth.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(3)
imagesc(interesse_1,interesse_2,mean(Std_growth,3))
title('Standard Deviation of Growth','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(3,'Std_Growth','fig')
print -depsc Std_Growth.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(4)
imagesc(interesse_1,interesse_2,mean(gini_tmax,3))
title('Gini Coefficient at t_{max}','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(4,'Gini_at_tmax','fig')
print -depsc Gini_at_tmax.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(5)
imagesc(interesse_1,interesse_2,mean(theil_tmax,3))
title('Theil Coefficient at t_{max}','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(5,'Theil_at_tmax','fig')
print -depsc Theil_at_tmax.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(6)
imagesc(interesse_1,interesse_2,mean(max_gini,3))
title('Gini Coefficient Maximum Value','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(6,'Gini_maxval','fig')
print -depsc Gini_maxval.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(7)
imagesc(interesse_1,interesse_2,mean(when_max_gini,3))
title('Gini Coefficient Time of Maximum Value','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(7,'Gini_time_maxval','fig')
print -depsc Gini_time_maxval.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(8)
imagesc(interesse_1,interesse_2,mean(max_theil,3))
title('Theil Coefficient Maximum Value','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(8,'Gini_maxval','fig')
print -depsc Gini_maxval.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(9)
imagesc(interesse_1,interesse_2,mean(when_max_theil,3))
title('Theil Coefficient Time of Maximum Value','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(9,'Theil_time_maxval','fig')
print -depsc Theil_time_maxval.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(10)
imagesc(interesse_1,interesse_2,mean(weighted_movements_tmax,3))
title('Weighted Movements Index (M_t) at t_{max}','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(10,'WeightedMovementsW_at_Tmax','fig')
print -depsc WeightedMovementsW_at_Tmax.eps


load('baseline_setup.mat')
load('Experiment_B')
figure(11)
imagesc(interesse_1,interesse_2,mean(meanVarRicchezzeabs,3))
title('Absolute Mean Wealth Change Index (| V_t |)','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(11,'MeanVarRicchezzeAbs','fig')
print -depsc MeanVarRicchezzeAbs.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(12)
imagesc(interesse_1,interesse_2,mean(meanVarRicchezze,3))
title('Time Average of Mean Wealth Change Index (V_t)','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(12,'MeanVarRicchezze','fig')
print -depsc MeanVarRicchezze.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(13)
imagesc(interesse_1,interesse_2,mean(FinalVarRicchezzeabs,3))
title('Absolute Mean Wealth Change Index (| V_t |) at tmax','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(13,'FinalVarRicchezzeAbs','fig')
print -depsc FinalVarRicchezzeAbs.eps

load('baseline_setup.mat')
load('Experiment_B')
figure(14)
imagesc(interesse_1,interesse_2,mean(FinalVarRicchezze,3))
title('Mean Wealth Change Index (V_t) at tmax','FontSize',19);
xlabel('Standard Deviation of Individual Returns','FontSize',19)
ylabel('Average Individual Return','FontSize',19);
set(gca,'YDir','normal');
colorbar
saveas(14,'FinalVarRicchezze','fig')
print -depsc FinalVarRicchezze.eps