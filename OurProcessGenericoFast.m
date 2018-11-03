% This file produces results of BiondiRighi (2018, JEIC) for a single
% simulation.

clear all
close all
clc

filename='prova';
filename=[filename '.mat'];
% PARAMETRI TIPO SIMULAZIONE
reddito=0;
risparmio=0;
reduction_of_interest='n'; % 'n' no reduction of interest in time
                          % 't' reduction along time (or log of time, just change the line below)
                          % 'w' reduction on the basis of total wealth.
type_of_interest='c'; % 'c' compound interest over time
                      % 's' simple interest over time
taxes_yes='n';  % yes: 'y' or no: 'n'

% PARAMETRI TASSE
tax_base='i'; % 'w' tassa su ricchezza (W_(t), 'i' tassa su income (W_t-W_(t-1))
tax_type='prog'; % 'prop': proportional 'prog': progessive,'lump': lumpsum
redistribution_type='key'; % 'key': invesely proportional to tax base, 'lib': T/N per tutti
taxrate=0.1; % for 'lump': percentage of W(1) extracted from all
             % for 'prop': proportion of the current tax base extracted
             %             from each agent
             % for 'prog': Maximum proportion of tax base extracted to the richest (for
             %             everybody else is a proportion of it).

% PARAMETRI POPOLAZIONE E TEMPO
N=5000; % number of agents
tmax=5000; % number of steps

% PARAMETRI DI RICCHEZZA E INCOME
interesse=0.05; %interesse di riferimento annuale
init_wealth_avg=10; % average of initial wealths
initial_wealth_type='equ'; % 'equ'
income_avg=1; % average of periodic income

% PARAMETRI DI DISEGNO
time_steps_length=1; % every how many step should i generate a transition plot?

interesse_giornaliero=interesse; %tasso interesse convertito in giorni con t = giorni
pd = makedist('Normal', 'mu',  0.05, 'sigma', 0.05);
%pd=makedist('Gamma','a',0.25,'b',0.2); % mean = a * b = 0.05 ; variance = a * b^2 = 0.05


%%%%%%%%%%%%%%%%%%%%%
% END OF PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%

% random seed
randseed=round(rand*1000000000); % The random seed is authomatically changed each run
randn('seed',randseed);
rand('seed',randseed);

% controllo ricchezze iniziali minina > 0
initialmin=init_wealth_avg/N;

% controllO se risparmio ? 1 allora lo ? anche reddito x definizione
if risparmio==1; reddito=1; end;

tic

% INIZIALIZZAZIONI
W=zeros(N,tmax);
if reddito==1
    Y=zeros(N,tmax);
    YT=zeros(N,1);
    share_WY=zeros(N,1);
    growth=zeros(N,1);
end
if risparmio==1
    WY=zeros(N,tmax);
end

if strcmp(initial_wealth_type,'lin')==1
    W(:,1)=linspace(initialmin,init_wealth_avg,N); % agent's initial wealth is split around the average above
else
    W(:,1)=init_wealth_avg;
end
            


if reddito==1
    % NOTE: if initial incomes are different than this is not ok... but why
    % should it be?
    %Y(:,1)=income_avg; % agent's initial periodic income is the average above 
    Y=income_avg.*rand(N,tmax); % agent's saving rate 'rand' is applied to each period j
    YT=sum(Y,2); % sum7 is the sume across time of Y 
    if risparmio==1
        WY(:,1)=Y(:,1);
    end
end

% CREATION OF THE VECTOR OF RETURNS
r=max(-1,random(pd,[N,tmax]));
if reduction_of_interest=='t'
    tt=zeros(N,tmax);
    for i=1:N
        tt(i,:)=[1:tmax];
    end
    r=r./(log(1+tt)); % either .^(1/t) or .^(1/(log(t+1))) or /log(1+t) se ^
end


if type_of_interest=='c'  % interesse composto
    for t=2:tmax % for each period t
        if reduction_of_interest=='w'
            if t>2
                r(:,t-1)=r(:,t-1)./log(1+sum(W(:,t-1),1));
            end
        end
        W(:,t)=max(0,W(:,t-1).*((1+r(:,t-1))));
        W(W(:,t-1)<=0,t)=0;
        if risparmio==1
            WY(:,t)= Y(:,t) + max(0,WY(:,t-1).*((1+r(:,t-1)))); % is it right, i was not summing for the W?i,t-1)<=0
        end
        if strcmp(taxes_yes,'y')
            W=taxes(W,t,N,tax_base,tax_type,redistribution_type,taxrate);
        end
    end
end

if type_of_interest=='s' % interesse semplice
    if risparmio==1
        sum3=zeros(N,1);
    end
    for t=2:tmax
        %t
        if reduction_of_interest=='w'
            if t>2
                if risparmio==0
                    r(:,t-1)=r(:,t-1)./log(1+sum(W(:,t-1),1));
                else
                    r(:,t-1)=r(:,t-1)./log(1+sum(W(:,t-1),1)+sum(WY(:,t-1),1));
                end
            end
         end
         if risparmio==1
             for n=1:t
                %sum4=sum(r(:,n:t),2);
                sum3=sum3+Y(:,n).*(1+sum(r(:,n:t),2));
             end
             WY(:,t)=max(0,sum3);
             WY(WY(:,t-1)<=0,t)=0;

         end
         sum2=sum(r(:,1:t-1),2);
         W(:,t)= max(0,W(:,1).*(1+sum2));
         W(W(:,t-1)<=0,t)= 0;
         if strcmp(taxes_yes,'y')
            W=taxes(W,t,N,tax_base,tax_type,redistribution_type,taxrate);
         end
    end
    if risparmio==1
        for t=1:tmax
            W(:,t)=W(:,t)+WY(:,t);
        end
    end
end

if reddito==1 && risparmio==0
    W=W+Y;     % aggiungo income come variabile autonoma additiva
end

for t=2:tmax
    growth(t)=(sum(W(:,t),1)-sum(W(:,t-1),1))/sum(W(:,t-1),1);
end

if reddito==1 && risparmio==0
    share_WY=sum(Y,1)./(sum(W,1)+sum(Y,1));
end
if risparmio==1
    share_WY=sum(WY,1)./(sum(WY,1)+sum(W,1));
end


M_growth=mean(growth);
Std_growth=std(growth);



% PRODUZIONE GRAFICI
% if reddito==0 && risparmio==0 % se non c'? il reddito procedi con i grafici
%     [P_logw_logr, P_logr_logw, P_r_w, P_logw_r, P_w_r, P_logw_logr_top, P_logr_logw_top, P_r_w_top, P_logw_r_top, P_w_r_top]=visualizzazioni(growth,NaN,NaN,NaN,W,N,tmax,time_steps_length,risparmio,reddito);
% end
% if reddito==1 && risparmio==0
%     [P_logw_logr, P_logr_logw, P_r_w, P_logw_r, P_w_r, P_logw_logr_top, P_logr_logw_top, P_r_w_top, P_logw_r_top, P_w_r_top]=visualizzazioni(growth,YT,share_WY,0,W,N,tmax,time_steps_length,risparmio,reddito);
% end
% if reddito==1 && risparmio==1
%     [P_logw_logr, P_logr_logw, P_r_w, P_logw_r, P_w_r, P_logw_logr_top, P_logr_logw_top, P_r_w_top, P_logw_r_top, P_w_r_top]=visualizzazioni(growth,YT,share_WY,WY,W,N,tmax,time_steps_length,risparmio,reddito);
% end


%countfigs=1;
%figure(countfigs)
%plot(growth(:));
%title('Total Wealth Growth Rate over time')

% if reddito==1
%     countfigs=countfigs+1;
%     figure(countfigs)
%     %bar(YT(:));
%     title('Total Income: Y_{i}^t')
%     
%     countfigs=countfigs+1;
%     figure(countfigs)
%     %histfit(YT(:),100)
%     title('Total Income Dist baseline')
% 
%     countfigs=countfigs+1;
%     figure(countfigs)
%     %plot(share_WY(:));
%     title('Share Wealth from Saved Income with reinvestment: Y_{i}^t')
% end

% if risparmio==1
%     
%     countfigs=countfigs+1;
%     figure(countfigs)
%     %subplot(1,2,1);
%     %hist([WY(:,tmax) W(:,tmax)-WY(:,tmax)],10);
%     title('Total Wealth from (Saved and Capital) Income: WY_{i} at tmax')
%     legend('Saved','Capital',0)
%     %subplot(1,2,2);
%     %figure(countfigs)
%     %hold on
%     %hist(W(:,tmax)-WY(:,tmax),10,'r');
%     %title('Total Wealth from (Capital) Income: W_{i}-WY_{i} at tmax')
%     
% end

% countfigs=countfigs+1;
% figure(countfigs)
% for ag=1:N
%     %plot([1:1:tmax],W(ag,:));
%     hold on
% end
% title('Wealths: W_{i}^t')


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



% countfigs=countfigs+1;
% figure(countfigs)
% %plot([1:1:tmax],gini_avg)
% hold on
% title('Gini index')
% legend('Gini')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% %plot([1:1:tmax],theil_avg)
% hold on
% title('Theil index')
% legend('Theil')
% 
% countfigs=countfigs+1;
% figure(countfigs)
% %histfit(W(:,tmax))
% title('Wealth Distribution baseline')


%countfigs=countfigs+1;
%figure(countfigs)
A=sort(W(:,tmax), 'descend');
B=[1:1:N];
%plot(log(B),log(A'),'r+')
C=A(1:round(N/10));
D=[1:1:round(N/10)];
P_logw_logr_top=polyfit(log(D),log(C'),1); % aggiungi top, cambia lettere
Y_top=polyval(P_logw_logr_top,log(D)); %aggiungi top x2 cambia lettere
%hold on 
%plot(log(D),Y_top,'k--') % top e lettere
% xlabel('log(rank)')
% ylabel('log(Wealth)')
% title('log(wealth) over log(rank)')
% legend('data','fit top 10%')
% 


if tmax>=1000
    %countfigs=countfigs+1;
    %figure(countfigs)
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


    %bar(wealthquartiles','histc'); 
%     title('Share fo Wealth for each quartile (poorest to richest)')
%     legend('t=1','t=10','t=100','t=1000','t=tmax',0);



    tint=[1, 10, 100, 1000, tmax];
    for j=1:length(tint)
        C=sort(W(:,tint(j)), 'ascend');
        for ii=0:9
            wealthdecile(j,ii+1)=sum(C((round(N/10)*ii)+1:(round(N/10))*(ii+1)))/sum(C);
        end
    end

%     countfigs=countfigs+1;
%     figure(countfigs)
%     %bar(wealthdecile','histc'); 
%     title('Share fo Wealth for each decile (poorest to richest)')
%     legend('t=1','t=10','t=100','t=1000','tmax',0);


end

% Plot of transition between deciles.

%countfigs=countfigs+1;
%figure(countfigs)
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

%countfigs=countfigs+1;
%plot(mean(D,1))
%title('progression of differences')

% countfigs=countfigs+1;
% figure(countfigs)
% pcolor(Weighted_changes_mtx_dec_w(:,:,2))
% colorbar
% title('2')

% countfigs=countfigs+1;
% figure(countfigs)
% pcolor(Weighted_changes_mtx_dec_w(:,:,tmax-1))
% colorbar
% title('tmax')

Varricchezzeabs=zeros(tmax,1);
Varricchezze=zeros(tmax,1);

for t=2:tmax
    Varricchezzeabs(t)=sum(abs((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t)))))./N;
    Varricchezze(t)=sum((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t))))./N;
end

% countfigs=countfigs+1;
% figure(countfigs) 
% %plot(Varricchezzeabs)
% title('Varricchezzeabs')
% 
% countfigs=countfigs+1;
% %figure(countfigs) 
 plot(Varricchezze)
% title('Varricchezze')



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
% countfigs=countfigs+1;
% figure(countfigs)
%plot(Weighted_movement_W)
% title('Weighted_movement_WW');
% 
% countfigs=countfigs+1;
% figure(countfigs)
% %plot(Weighted_movement_Wnoabs)
% title('Weighted_movement_Wnoabs')


window_of_mav=10;
for i=1:tmax-window_of_mav-1 % construct the moving averages
    mav_sum_movement_tot(i)=sum(sum_movement_tot(i:i+window_of_mav-1))./(window_of_mav);
    mav_sum_movement_up_tot(i)=sum(sum_movement_up_tot(i:i+window_of_mav-1))./((window_of_mav));
    mav_sum_movement_down_tot(i)=sum(sum_movement_down_tot(i:i+window_of_mav-1))./((window_of_mav));
    mav_sum_movement_same_tot(i)=sum(sum_movement_same_tot(i:i+window_of_mav-1))./((window_of_mav));
    mav_Weighted_movement_W(i)=sum(Weighted_movement_W(i:i+window_of_mav-1))./((window_of_mav));
    for j=1:10
        mav_sum_movement(j,i)=sum(sum_movement(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_up(j,i)=sum(sum_movement_up(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_down(j,i)=sum(sum_movement_down(j,i:i+window_of_mav-1))./((window_of_mav));
        mav_sum_movement_same(j,i)=sum(sum_movement_same(j,i:i+window_of_mav-1))./((window_of_mav));
    end      
end
% countfigs=countfigs+1;
% figure(countfigs)
% %plot(mav_sum_movement_tot)
% hold on
% %plot(mav_sum_movement_up_tot,'r')
% hold on
% %plot(mav_sum_movement_down_tot,'k')
% hold on
% %plot(mav_sum_movement_same_tot,'m')
% xlabel('time')
% ylabel('proportion')
% ylim([0 1]);
% title('movements (moving averages)')
% legend('total','upward','downward','not moving')

% countfigs=countfigs+1;
% figure(countfigs)
% t=1:time_steps_length:tmax-1;
% plot(weighted_movements'./N,'k')
% hold on
for i=1:tmax-10-1
    mav_weighted(i)=sum(weighted_movements(i:i+10)./(N*10));
end
%plot(mav_weighted','r+')
% xlabel('time')
% ylabel('Weighted Movements')
% title('Weigthed movements across deciles (per agent per period \in [0,1])')



%save(filename, 'gini_avg','theil_avg','W',...
%    'tmax','N','A','B','P_logw_logr_top','wealthdecile','Varricchezze',...
%    'Varricchezzeabs','Weighted_movement_W','WY','share_WY','M_growth','Std_growth','YT')

P_logw_logr_top
M_growth
Std_growth
gini_tmax
theil_tmax
when_max_gini
max_gini
when_max_theil
max_theil
weighted_movements_tmax

%save(filename)
toc
close all