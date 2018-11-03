function [P_logw_logr, P_logr_logw, P_r_w, P_logw_r, P_w_r, P_logw_logr_top, P_logr_logw_top, P_r_w_top, P_logw_r_top, P_w_r_top] = visualizzazioni(growth,YT,share_WY,WY,W,N,tmax,time_steps_length,risparmio, reddito)
% Sample visualizations for single simulation runs.
% Biondi Righi (2018, JEIC)


P_logr_logw_top=NaN;
P_r_w_top=NaN;

countfigs=1;
figure(countfigs)
plot(growth(:));
title('Total Wealth Growth Rate over time')

if reddito==1
    countfigs=countfigs+1;
    figure(countfigs)
    bar(YT(:));
    title('Total Income: Y_{i}^t')
    
    countfigs=countfigs+1;
    figure(countfigs)
    histfit(YT(:),100)
    title('Total Income Dist baseline')

    countfigs=countfigs+1;
    figure(countfigs)
    plot(share_WY(:));
    title('Share Wealth from Saved Income with reinvestment: Y_{i}^t')
end

if risparmio==1
    
    countfigs=countfigs+1;
    figure(countfigs)
    %subplot(1,2,1);
    hist([WY(:,tmax) W(:,tmax)-WY(:,tmax)],10);
    title('Total Wealth from (Saved and Capital) Income: WY_{i} at tmax')
    legend('Saved','Capital',0)
    %subplot(1,2,2);
    %figure(countfigs)
    %hold on
    %hist(W(:,tmax)-WY(:,tmax),10,'r');
    %title('Total Wealth from (Capital) Income: W_{i}-WY_{i} at tmax')
    
end



countfigs=countfigs+1;
figure(countfigs)
for ag=1:N
    plot([1:1:tmax],W(ag,:));
    hold on
end
title('Wealths: W_{i}^t')


[gini_avg]=ginicoeff(W',2,1);
%[theil_avg] = theilt(W);
%theil_avg=theil_avg./log(N);
theil_avg=zeros(tmax,1);
for t=1:tmax
    data=W(:,t);
    val=(data./mean(data)).*(log(data./mean(data)));
    num=sum(val);
    den=N*log(N);
    theil_avg(t)=num/den;
end



countfigs=countfigs+1;
figure(countfigs)
plot([1:1:tmax],gini_avg)
hold on
title('Gini index')
legend('Gini')

countfigs=countfigs+1;
figure(countfigs)
plot([1:1:tmax],theil_avg)
hold on
title('Theil index')
legend('Theil')

countfigs=countfigs+1;
figure(countfigs)
histfit(W(:,tmax),100)
title('Wealth Distribution baseline')

% countfigs=countfigs+1;
% figure(countfigs)
% aa=(W(:,tmax)~=0);
% histfit(log(W(aa,tmax)),100)
% title('Log Wealth Dist baseline')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% histfit(W(:,tmax)/mean(W(:,tmax)),100)
% title('Standardised Wealth Dist baseline')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% aa=(W(:,tmax)~=0);
% histfit(log(W(aa,tmax))/mean(W(:,tmax)),100)
% title('Standardised Log Wealth Dist baseline')


P_logw_logr=0;
P_logr_logw=0;
P_r_w=0;
P_w_r=0;
P_logw_r=0;

countfigs=countfigs+1;
figure(countfigs)
A=sort(W(:,tmax), 'descend');
B=[1:1:N];
plot(log(B),log(A'),'r+')
%P_logw_logr=polyfit(log(B),log(A'),1);
%Y=polyval(P_logw_logr,log(B));
%hold on 
%plot(log(B),Y)
C=A(1:round(N/10));
D=[1:1:round(N/10)];
P_logw_logr_top=polyfit(log(D),log(C'),1); % aggiungi top, cambia lettere
Y_top=polyval(P_logw_logr_top,log(D)); %aggiungi top x2 cambia lettere
hold on 
plot(log(D),Y_top,'k--') % top e lettere
xlabel('log(rank)')
ylabel('log(Wealth)')
title('log(wealth) over log(rank)')
%legend('data','fit','fit top 10%')
legend('data','fit top 10%')

% countfigs=countfigs+1;
% figure(countfigs)
% A=sort(W(:,tmax), 'descend');
% B=[1:1:N];
% plot(log(A),log(B'),'r+');
% %P_logr_logw=polyfit(log(A),log(B'),1);
% %Y=polyval(P_logr_logw,log(A));
% %hold on 
% %plot(log(A),Y)
% C=A(1:round(N/10));
% D=[1:1:round(N/10)];
% P_logr_logw_top=polyfit(log(C),log(D'),1);
% Y=polyval(P_logr_logw_top,log(C));
% hold on 
% plot(log(C),Y_top,'k--')
% xlabel('log(Wealth)')
% ylabel('log(rank)')
% title('log(rank) over log(wealth)')
% %legend('data','fit','fit top 10%')
% legend('data','fit top 10%')

% countfigs=countfigs+1;
% figure(countfigs)
% A=sort(W(:,tmax), 'descend');
% B=[1:1:N];
% plot(A,B','r+')
% %P_r_w=polyfit(A,B',1);
% %Y=polyval(P_r_w,A);
% %hold on 
% %plot(A,Y)
% C=A(1:round(N/10));
% D=[1:1:round(N/10)];
% P_r_w_top=polyfit(C,D',1);
% Y_top=polyval(P_r_w_top,C);
% hold on 
% plot(C,Y_top,'k--');
% xlabel('wealth')
% ylabel('rank')
% title('(rank) over (wealth)')
% %legend('data','fit','fit top 10%')
% legend('data','fit top 10%')


countfigs=countfigs+1;
figure(countfigs)
A=sort(W(:,tmax), 'descend');
B=[1:1:N];
plot(B,A','r+');
%P_w_r=polyfit(B,A',1);
%Y=polyval(P_w_r,B);
%hold on 
%plot(B,Y)
C=A(1:round(N/10));
D=[1:1:round(N/10)];
P_w_r_top=polyfit(D,C',1);
Y_top=polyval(P_w_r_top,D);
hold on 
plot(D,Y_top,'k--')
xlabel('rank')
ylabel('wealth')
title('(wealth) over (rank)')
%legend('data','fit','fit top 10%')
legend('data','fit top 10%')


countfigs=countfigs+1;
figure(countfigs)
A=sort(W(:,tmax), 'descend');
B=[1:1:N];
plot((B),log(A'),'r+');
%P_logw_r=polyfit(B,log(A'),1);
%Y=polyval(P_logw_r,B);
%hold on 
%plot(B,Y);
C=A(1:round(N/10));
D=[1:1:round(N/10)];
P_logw_r_top=polyfit(D,log(C'),1);
Y_top=polyval(P_logw_r_top,D);
hold on 
plot(D,Y_top,'k--');
xlabel('rank')
ylabel('log(wealth)')
title('log(wealth) over (rank)');
%legend('data','fit','fit top 10%')
legend('data','fit top 10%')

if tmax>=1000
countfigs=countfigs+1;
figure(countfigs)
% t=1
C=sort(W(:,1), 'descend');
wealthquartiles(1,4)=sum(C(1:round(N/4)))/sum(C); %top25
wealthquartiles(1,3)=sum(C(round(N/4)+1:round(N/4)*2))/sum(C); %second25
wealthquartiles(1,2)=sum(C((round(N/4)*2)+1:round(N/4)*3))/sum(C); % third25
wealthquartiles(1,1)=sum(C(round(N/4)*3+1:end))/sum(C); % low25

%t=10
C=sort(W(:,10), 'descend');
wealthquartiles(2,4)=sum(C(1:round(N/4)))/sum(C); %top25
wealthquartiles(2,3)=sum(C(round(N/4)+1:round(N/4)*2))/sum(C); %second25
wealthquartiles(2,2)=sum(C((round(N/4)*2)+1:round(N/4)*3))/sum(C); % third25
wealthquartiles(2,1)=sum(C(round(N/4)*3+1:end))/sum(C); % low25

%t=100
C=sort(W(:,100), 'descend');
wealthquartiles(3,4)=sum(C(1:round(N/4)))/sum(C); %top25
wealthquartiles(3,3)=sum(C(round(N/4)+1:round(N/4)*2))/sum(C); %second25
wealthquartiles(3,2)=sum(C((round(N/4)*2)+1:round(N/4)*3))/sum(C); % third25
wealthquartiles(3,1)=sum(C(round(N/4)*3+1:end))/sum(C); % low25

% t=1000
C=sort(W(:,1000), 'descend');
wealthquartiles(4,4)=sum(C(1:round(N/4)))/sum(C); %top25
wealthquartiles(4,3)=sum(C(round(N/4)+1:round(N/4)*2))/sum(C); %second25
wealthquartiles(4,2)=sum(C((round(N/4)*2)+1:round(N/4)*3))/sum(C); % third25
wealthquartiles(4,1)=sum(C(round(N/4)*3:end))/sum(C); % low25

% tmax
wealthquartiles(5,4)=sum(A(1:round(N/4)))/sum(A); %top25
wealthquartiles(5,3)=sum(A(round(N/4)+1:round(N/4)*2))/sum(A); %second25
wealthquartiles(5,2)=sum(A((round(N/4)*2)+1:round(N/4)*3))/sum(A); % third25
wealthquartiles(5,1)=sum(A(round(N/4)*3:end))/sum(A); % low25


bar(wealthquartiles','histc'); 
title('Share fo Wealth for each quartile (poorest to richest)')
legend('t=1','t=10','t=100','t=1000','t=tmax',0);



tint=[1, 10, 100, 1000, tmax];
for j=1:length(tint)
    C=sort(W(:,tint(j)), 'ascend');
    for ii=0:9
        wealthdecile(j,ii+1)=sum(C((round(N/10)*ii)+1:(round(N/10))*(ii+1)))/sum(C);
    end
end

countfigs=countfigs+1;
figure(countfigs)
bar(wealthdecile','histc'); 
title('Share fo Wealth for each decile (poorest to richest)')
legend('t=1','t=10','t=100','t=1000','tmax',0);

end

% Plot of transition between deciles.

countfigs=countfigs+1;
changes_mtx_dec=zeros(10,10,tmax-1);
Weighted_changes_mtx_dec=zeros(10,10,tmax/time_steps_length);
Weighted_changes_mtx_dec_w=zeros(10,10,tmax/time_steps_length);
Weighted_changes_mtx_dec_w_noabs=zeros(10,10,tmax/time_steps_length);


for i=1:tmax-1 % for every time step
  [~, B]=sort(W(:,i), 'descend');
  [WW, E]=sort(W(:,i+1), 'descend');
  DecW=zeros(10,1); % create the deciles averages
  Declength=N/10;
  for q=1:10
      DecW(q)=median(WW(((q-1)*Declength+1):(q*Declength),1));
  end
  for j=1:N
      D(j,i)=abs(B(j)-E(j))./N;
      dec_init=floor(((B(j)-1)/N)*10)+1;
      if dec_init>10
          dec_init=10;
      end
      dec_end=floor(((E(j)-1)/N)*10)+1;
      if dec_end>10
          dec_end=10;
      end
      changes_mtx_dec(dec_init,dec_end,i)=changes_mtx_dec(dec_init,dec_end,i)+1;
      Weighted_changes_mtx_dec(dec_init,dec_end,i)=Weighted_changes_mtx_dec(dec_init,dec_end,i)+(abs(dec_init-dec_end)/9);

      %Weighted_changes_mtx_dec_w(dec_init,dec_end,i)=Weighted_changes_mtx_dec_w(dec_init,dec_end,i)+(abs(DecW(dec_init)-DecW(dec_end)))/sum(DecW);
      Weighted_changes_mtx_dec_w(dec_init,dec_end,i)=Weighted_changes_mtx_dec_w(dec_init,dec_end,i)+(abs(DecW(dec_init)-DecW(dec_end)))/(DecW(1)-DecW(10));
      Weighted_changes_mtx_dec_w_noabs(dec_init,dec_end,i)=Weighted_changes_mtx_dec_w_noabs(dec_init,dec_end,i)+((DecW(dec_init)-DecW(dec_end)))/(DecW(1)-DecW(10));

  end
  
end
changes_mtx_dec_norm=changes_mtx_dec./N;
Weighted_changes_mtx_dec_w=Weighted_changes_mtx_dec_w./N;
Weighted_changes_mtx_dec_w_noabs=Weighted_changes_mtx_dec_w_noabs./N;

countfigs=countfigs+1;
plot(mean(D,1))
title('progression of differences')

countfigs=countfigs+1;
figure(countfigs)
pcolor(Weighted_changes_mtx_dec_w(:,:,2))
colorbar
title('2')

countfigs=countfigs+1;
figure(countfigs)
pcolor(Weighted_changes_mtx_dec_w(:,:,tmax-1))
colorbar
title('tmax')

Varricchezzeabs=zeros(tmax,1);
Varricchezze=zeros(tmax,1);

for t=2:tmax
    Varricchezzeabs(t)=sum(abs((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t)))))./N;
    Varricchezze(t)=sum((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t))))./N;
end

countfigs=countfigs+1;
figure(countfigs) 
plot(Varricchezzeabs)
title('Varricchezzeabs')

countfigs=countfigs+1;
figure(countfigs) 
plot(Varricchezze)
title('Varricchezze')


% countfigs=countfigs+1;
% figure(countfigs)
% [~, ind]=sort(W(:,round(tmax/2)), 'descend');
% for j=1:N % for every agent
%     ix=ind(j); % take note of the index
%     for t=round(tmax/2):tmax 
%         [~, E]=sort(W(:,t), 'descend');
%         rank_richest(t)=find(E==ix);
%     end
%     hold on
%     if ix~=1 && ix~=N
%         plot(round(tmax/2):1:tmax,rank_richest(round(tmax/2):tmax));
%     end
%     if ix==1
%         plot(round(tmax/2):1:tmax,rank_richest(round(tmax/2):tmax),'r');
%     end
%     if ix==N
%         plot(round(tmax/2):1:tmax,rank_richest(round(tmax/2):tmax),'k');
%     end    
% end
% xlabel('time')
% ylabel('rank')
% title('Evolution of Ranking')

[x,y]=size(changes_mtx_dec(:,:,1));
T=(tmax-1);
sum_movement=zeros(x,T);
sum_movement_up=zeros(x,T);
sum_movement_down=zeros(x,T);
sum_movement_same=zeros(x,T);

sum_movement_tot=zeros(T,1);
sum_movement_up_tot=zeros(T,1);
sum_movement_down_tot=zeros(T,1);
sum_movement_same_tot=zeros(T,1);

weighted_movements=zeros(T,1);
weighted_movements_up=zeros(T,1);
weighted_movements_down=zeros(T,1);

Weighted_movement_W=zeros(T,1);
Weighted_movement_Wnoabs=zeros(T,1);

for t=1:T
    for i=1:x
        for j=1:y
            weighted_movements(t)=weighted_movements(t)+Weighted_changes_mtx_dec(i,j,t);
            if i~=j % movimenti 
                sum_movement(i,t)=sum_movement(i,t)+changes_mtx_dec_norm(i,j,t);
                sum_movement_tot(t)=sum_movement_tot(t)+changes_mtx_dec_norm(i,j,t);
                Weighted_movement_W(t)=Weighted_movement_W(t)+Weighted_changes_mtx_dec_w(i,j,t);
                Weighted_movement_Wnoabs(t)=Weighted_movement_Wnoabs(t)+Weighted_changes_mtx_dec_w_noabs(i,j,t);
            end
            if i>j % movimenti su
                sum_movement_up(i,t)=sum_movement_up(i,t)+changes_mtx_dec_norm(i,j,t);
                sum_movement_up_tot(t)=sum_movement_up_tot(t)+changes_mtx_dec_norm(i,j,t);
                weighted_movements_up(t)=weighted_movements_up(t)+Weighted_changes_mtx_dec(i,j,t);
            end
            if i<j % movimenti giu
            	sum_movement_down(i,t)=sum_movement_down(i,t)+changes_mtx_dec_norm(i,j,t);
                sum_movement_down_tot(t)=sum_movement_down_tot(t)+changes_mtx_dec_norm(i,j,t);
                weighted_movements_down(t)=weighted_movements_down(t)+Weighted_changes_mtx_dec(i,j,t);
            end
            if i==j % non movimenti
             	sum_movement_same(i,t)=sum_movement_same(i,t)+changes_mtx_dec_norm(i,j,t);
                sum_movement_same_tot(t)=sum_movement_same_tot(t)+changes_mtx_dec_norm(i,j,t);
            end
        end
    end
end
countfigs=countfigs+1;
figure(countfigs)
plot(Weighted_movement_W)
title('Weighted_movement_WW');

countfigs=countfigs+1;
figure(countfigs)
plot(Weighted_movement_Wnoabs)
title('Weighted_movement_Wnoabs')

% countfigs=countfigs+1;
% figure(countfigs)
% plot(sum_movement_tot)
% hold on
% plot(sum_movement_up_tot,'r')
% hold on
% plot(sum_movement_down_tot,'k')
% hold on
% plot(sum_movement_same_tot,'m')
% title('movements')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% legend('total','upward','downward','not moving')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(t,sum_movement'.*10)
% title('movements (each decile)')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(t,sum_movement_up'.*10)
% title('movements UP (each decile)')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(t,sum_movement_down'.*10)
% title('movements DOWN (each decile)')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(t,sum_movement_same'.*10)
% title('NON movements (each decile)')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(t,sum_movement_same'.*N)
% xlabel('time')
% ylabel('Number')
% ylim([0 N/10]);
% title('NON movements (each decile): NUMBER')
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)

window_of_mav=10;
for i=1:tmax-window_of_mav-1 % construct the moving averages
    mav_sum_movement_tot(i)=sum(sum_movement_tot(i:i+window_of_mav-1))./(window_of_mav);
    mav_sum_movement_up_tot(i)=sum(sum_movement_up_tot(i:i+window_of_mav-1))./((window_of_mav));
    mav_sum_movement_down_tot(i)=sum(sum_movement_down_tot(i:i+window_of_mav-1))./((window_of_mav));
    mav_sum_movement_same_tot(i)=sum(sum_movement_same_tot(i:i+window_of_mav-1))./((window_of_mav));
    for j=1:10
        mav_sum_movement(j,i)=sum(sum_movement(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_up(j,i)=sum(sum_movement_up(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_down(j,i)=sum(sum_movement_down(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_same(j,i)=sum(sum_movement_same(j,i:i+window_of_mav-1))./((window_of_mav));
    end      
end
countfigs=countfigs+1;
figure(countfigs)
plot(mav_sum_movement_tot)
hold on
plot(mav_sum_movement_up_tot,'r')
hold on
plot(mav_sum_movement_down_tot,'k')
hold on
plot(mav_sum_movement_same_tot,'m')
xlabel('time')
ylabel('proportion')
ylim([0 1]);
title('movements (moving averages)')
legend('total','upward','downward','not moving')


% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(mav_sum_movement'.*10)
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% title('movements (each decile): moving averages (10 steps)')
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(mav_sum_movement_up'.*10)
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% title('movements UP (each decile): moving averages (10 steps)')
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(mav_sum_movement_down'.*10)
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% title('movements DOWN (each decile): moving averages (10 steps)')
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=[1:time_steps_length:tmax-1];
% plot(mav_sum_movement_same'.*10)
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% title('NON movements (each decile): moving averages (10 steps)')
% legend('1st Decile (richest)','2nd Decile','3rd Decile','4th Decile','5th Decile','6th Decile','7th Decile','8th Decile','9th Decile','10th Decile (poorest)',0)

countfigs=countfigs+1;
figure(countfigs)
t=1:time_steps_length:tmax-1;
plot(weighted_movements'./N,'k')
hold on
for i=1:tmax-10-1
    mav_weighted(i)=sum(weighted_movements(i:i+10)./(N*10));
end
plot(mav_weighted','r+')
xlabel('time')
ylabel('Weighted Movements')
title('Weigthed movements across deciles (per agent per period \in [0,1])')

% countfigs=countfigs+1;
% figure(countfigs)
% t=1:time_steps_length:tmax-1;
% plot(weighted_movements_up'./N,'k')
% hold on
% for i=1:tmax-10-1
%     mav_weighted_up(i)=sum(weighted_movements_up(i:i+10)./(N*10));
% end
% plot(mav_weighted_up','r+')
% xlabel('time')
% ylabel('Weighted Movements UP')
% title('Weigthed movements UP across deciles ')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% t=1:time_steps_length:tmax-1;
% plot(weighted_movements_down'./N,'k')
% for i=1:tmax-10-1
%     mav_weighted_down(i)=sum(weighted_movements_down(i:i+10)./(N*10));
% end
% plot(mav_weighted_down','r+')
% xlabel('time')
% ylabel('Weighted Movements DOWN')
% title('Weigthed movements DOWN across deciles ')

% countfigs=countfigs=1
% figure(countfigs)
% for t=1:tmax
%     aaa(t)=(max(r(:,t))*(median(W(:,t))))/(quantile(W(:,t),0.75)-median(W(:,t)));
% end
% plot(aaa(2:end))
% title('Potential social ascension')
