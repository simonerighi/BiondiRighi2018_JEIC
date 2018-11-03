% This file creates data and visualizations for figure 3 in Biondi Righi
% (2018, JEIC).

clear all
close all
clc

load('baseline_setup.mat');


N=20000; % number of agents

interesse_1=0.05; %average interest   (average o a) 0.25
interesse_2=0.05; %std interest (sigma o b) 0.2

[gini_avg,theil_avg, mav_Weighted_movement_W, Varricchezzeabs,...
    Varricchezze, share_WY, Final_wealth,A_sf,B_sf,C_sf,D_sf, Y_top,...
    wealthquartiles, wealthdecile, P_logw_logr_top, M_growth, Std_growth,...
    gini_tmax,theil_tmax,when_max_gini, max_gini, when_max_theil, max_theil,...
    weighted_movements_tmax,MeanTaxRate,MedianTaxRate,MeanRedistributionRate,MedianRedistributionRate,...
    Proportion_total_wealth,Proportion_relative_wealth,W]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed);

%figure(12)
ntop=round(N/100);
timportant=10;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(12,'Exp_AC_ranktop20_t10','fig')
%print -depsc Exp_AC_ranktop20_t10.eps
avg_agent_rank(1,:)=avg_pos_i;
total_agents_rank(1,:)=avg_pos_cum;
time_avg(1,:)=avg_pos_t;
time_std(1,:)=std_pos_t;
 
%figure(13)
ntop=round(N/100);
timportant=100;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(13,'Exp_AC_ranktop20_t100','fig')
%print -depsc Exp_AC_ranktop20_t100.eps
avg_agent_rank(2,:)=avg_pos_i;
total_agents_rank(2,:)=avg_pos_cum;
time_avg(2,:)=avg_pos_t;
time_std(2,:)=std_pos_t;

%figure(14)
ntop=round(N/100);
timportant=1000;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(14,'Exp_AC_ranktop20_t1000','fig')
%print -depsc Exp_AC_ranktop20_t1000.eps
avg_agent_rank(3,:)=avg_pos_i;
total_agents_rank(3,:)=avg_pos_cum;
time_avg(3,:)=avg_pos_t;
time_std(3,:)=std_pos_t;

%figure(15)
ntop=round(N/100);
timportant=2000;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(15,'Exp_AC_ranktop20_t2000','fig')
%print -depsc Exp_AC_ranktop20_t2000.eps
avg_agent_rank(4,:)=avg_pos_i;
total_agents_rank(4,:)=avg_pos_cum;
time_avg(4,:)=avg_pos_t;
time_std(4,:)=std_pos_t;

%figure(16)
ntop=round(N/100);
timportant=3000;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(16,'Exp_AC_ranktop20_t3000','fig')
%print -depsc Exp_AC_ranktop20_t3000.eps
avg_agent_rank(5,:)=avg_pos_i;
total_agents_rank(5,:)=avg_pos_cum;
time_avg(5,:)=avg_pos_t;
time_std(5,:)=std_pos_t;

%figure(17)
ntop=round(N/100);
timportant=4000;
traccia_ricchezze
%title(['Rank of Wealth, top 1% richest agents at t=' num2str(timportant)]);
%saveas(17,'Exp_AC_ranktop20_t4000','fig')
%print -depsc Exp_AC_ranktop20_t4000.eps
avg_agent_rank(6,:)=avg_pos_i;
total_agents_rank(6,:)=avg_pos_cum;
time_avg(6,:)=avg_pos_t;
time_std(6,:)=std_pos_t;

close all

% media e varianza del top 1% nel tempo (time_avg e time_std errorbar vs tempo)
figure(20)
errorbar(1:1:1000,time_avg(1,:),time_std(1,:),'r')
hold on
errorbar(1:1:1000,time_avg(2,:),time_std(2,:),'b')
hold on
errorbar(1:1:1000,time_avg(3,:),time_std(3,:),'k')
hold on
errorbar(1:1:1000,time_avg(4,:),time_std(4,:),'m')
hold on
errorbar(1:1:1000,time_avg(5,:),time_std(5,:),'c')
hold on
errorbar(1:1:1000,time_avg(6,:),time_std(6,:),'g')
hold on
title('Evolution of average rank of 1%');
ylabel('Average Rank')
xlabel('Time from measurement')
legend('Top 1% of t=10','Top 1% of t=100','Top 1% of t=1000','Top 1% of t=2000','Top 1% of t=3000','Top 1% of t=4000',0);
saveas(20,'Exp_AC_ranktop20_errobar_in_time','fig')
print -depsc Exp_AC_ranktop20_errobar_in_time.eps



% total_agents_rank distribuzioni
figure(18)
[x,y]=hist(total_agents_rank(1,:),20);
plot(y,x,'r');
hold on
[x,y]=hist(total_agents_rank(2,:),20);
plot(y,x,'b');
hold on
[x,y]=hist(total_agents_rank(3,:),20);
plot(y,x,'k');
hold on
[x,y]=hist(total_agents_rank(4,:),20);
plot(y,x,'m');
hold on
[x,y]=hist(total_agents_rank(5,:),20);
plot(y,x,'c');
hold on
[x,y]=hist(total_agents_rank(6,:),20);
plot(y,x,'g');
hold on
title('Distribution of average rank in the next 1000 steps, top 1% richest agents');
xlabel('Average Rank')
ylabel('Frequency')
legend('mean rank in t \in [10-1010]','mean rank in t \in [100-1100]','mean rank in t \in [1000-2000]','mean rank in t \in [2000-3000]',...
    'mean rank in t \in [3000-4000]','mean rank in t \in [4000-5000]',0);
saveas(18,'Exp_AC_ranktop20_distroavgposition(all)','fig')
print -depsc Exp_AC_ranktop20_distroavgposition(all).eps



% avg_agent_rank distribuzioni

figure(19)
[x,y]=hist(avg_agent_rank(1,:),20);
semilogx(y/N,x/ntop,'r');
hold on
[x,y]=hist(avg_agent_rank(2,:),20);
semilogx(y/N,x/ntop,'b');
hold on
[x,y]=hist(avg_agent_rank(3,:),20);
semilogx(y/N,x/ntop,'k');
hold on
[x,y]=hist(avg_agent_rank(4,:),20);
semilogx(y/N,x/ntop,'m');
hold on
[x,y]=hist(avg_agent_rank(5,:),20);
semilogx(y/N,x/ntop,'c');
hold on
[x,y]=hist(avg_agent_rank(6,:),20);
semilogx(y/N,x/ntop,'g');
hold on
title('PDF of Average rank over the next 1000 time steps, top 1% richest agents at period t','FontSize',19);
xlabel('Normalized Average Rank')
ylabel('Density')
legend('mean rank in t \in [10-1010]','mean rank in t \in [100-1100]','mean rank in t \in [1000-2000]','mean rank in t \in [2000-3000]',...
    'mean rank in t \in [3000-4000]','mean rank in t \in [4000-5000]',0);
saveas(19,'Exp_AC_ranktop20_distroavgposition','fig')
print -depsc Exp_AC_ranktop20_distroavgposition.eps


figure(21)
plot(1:1:1000,time_avg(1,:),'r')
hold on
plot(1:1:1000,time_avg(2,:),'b')
hold on
plot(1:1:1000,time_avg(3,:),'k')
hold on
plot(1:1:1000,time_avg(4,:),'m')
hold on
plot(1:1:1000,time_avg(5,:),'c')
hold on
plot(1:1:1000,time_avg(6,:),'g')
hold on
title('Evolution of average rank of 1%');
ylabel('Average Rank')
xlabel('Time from measurement')
legend('Top 1% of t=10','Top 1% of t=100','Top 1% of t=1000','Top 1% of t=2000','Top 1% of t=3000','Top 1% of t=4000',0);
saveas(21,'Exp_AC_ranktop20_plot_in_time','fig')
print -depsc Exp_AC_ranktop20_plot_in_time.eps

% figure(22)
% semilogx(1:1:1000,time_avg(1,:),'r')
% hold on
% semilogx(1:1:1000,time_avg(2,:),'b')
% hold on
% semilogx(1:1:1000,time_avg(3,:),'k')
% hold on
% semilogx(1:1:1000,time_avg(4,:),'m')
% hold on
% semilogx(1:1:1000,time_avg(5,:),'c')
% hold on
% semilogx(1:1:1000,time_avg(6,:),'g')
% hold on
% title('Evolution of average rank of 1%');
% ylabel('Average Rank')
% xlabel('Log(Time from measurement)')
% legend('Top 1% of t=10','Top 1% of t=100','Top 1% of t=1000','Top 1% of t=2000','Top 1% of t=3000','Top 1% of t=4000',0);
% saveas(22,'Exp_AC_ranktop20_logx_in_time','fig')
% print -depsc Exp_AC_ranktop20_logx_in_time.eps
% 
% figure(23)
% semilogy(1:1:1000,time_avg(1,:),'r')
% hold on
% semilogy(1:1:1000,time_avg(2,:),'b')
% hold on
% semilogy(1:1:1000,time_avg(3,:),'k')
% hold on
% semilogy(1:1:1000,time_avg(4,:),'m')
% hold on
% semilogy(1:1:1000,time_avg(5,:),'c')
% hold on
% semilogy(1:1:1000,time_avg(6,:),'g')
% hold on
% title('Evolution of average rank of 1%');
% ylabel('Log(Average Rank)')
% xlabel('Time from measurement')
% legend('Top 1% of t=10','Top 1% of t=100','Top 1% of t=1000','Top 1% of t=2000','Top 1% of t=3000','Top 1% of t=4000',0);
% saveas(23,'Exp_AC_ranktop20_logy_in_time','fig')
% print -depsc Exp_AC_ranktop20_logy_in_time.eps

vct_1(1:1000,1)=ntop;
vct_2(1:1000,1)=ntop*10;

figure(24)
loglog(1:1:1000,time_avg(1,:)/N,'r')
hold on
loglog(1:1:1000,time_avg(2,:)/N,'b')
hold on
loglog(1:1:1000,time_avg(3,:)/N,'k')
hold on
loglog(1:1:1000,time_avg(4,:)/N,'m')
hold on
loglog(1:1:1000,time_avg(5,:)/N,'c')
hold on
loglog(1:1:1000,time_avg(6,:)/N,'g')
hold on
loglog(1:1:1000,vct_1/N,'Color',[0.8 0.8 0.8]);
hold on
loglog(1:1:1000,vct_2/N,'Color',[0.7 0.7 0.7]);
title('Evolution of average rank of 1%','FontSize',19);
ylabel('Log(Normalized Rank)')
xlabel('Log(Time from measurement)')
legend('Top 1% of t=10','Top 1% of t=100','Top 1% of t=1000','Top 1% of t=2000','Top 1% of t=3000','Top 1% of t=4000',...
    'First Percentile','First Decile',0);
saveas(24,'Exp_AC_ranktop20_loglog_in_time','fig')
print -depsc Exp_AC_ranktop20_loglog_in_time.eps

clear W