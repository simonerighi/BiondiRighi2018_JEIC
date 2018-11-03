function [gini_avg,theil_avg, mav_Weighted_movement_W, Varricchezzeabs,...
    Varricchezze, share_WY, Final_wealth,A_sf,B_sf,C_sf,D_sf, Y_top,...
    wealthquartiles, wealthdecile, P_logw_logr_top, M_growth, Std_growth,...
    gini_tmax,theil_tmax,when_max_gini, max_gini, when_max_theil, max_theil,...
    weighted_movements_tmax,MeanTaxRate,MedianTaxRate,MeanRedistributionRate,MedianRedistributionRate,Proportion_total_wealth,...
    Proportion_relative_wealth,W,ratio_top1T_bottom10W,ratio_top10T_bottom50W,ratio_top50T_bottom50W,WatZero_end]=...
    OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed)

% This function generates the results for single iterations of the
% simuations of Biondi and Righi (2018, JEIC): Inequality, Mobility and the
% financial accumulation process: A computational economic analysis.
% with different input paramters it is possible to obtain results
% concerning different processes of wealth accumulation.



% Design parameters
time_steps_length=1; % every how many step should i generate a transition plot?

if tipo_interesse=='Norma';pd = makedist('Normal', 'mu',  interesse_1, 'sigma', interesse_2); end;
if tipo_interesse=='Gamma';pd = makedist('Gamma','a',interesse_1,'b',interesse_2); end; 
    % mean = a * b = 0.05 ; variance = a * b^2 = 0.05

MeanTaxRate=0;
MedianTaxRate=0;
MeanRedistributionRate=0;
MedianRedistributionRate=0;

%%%%%%%%%%%%%%%%%%%%%
% END OF PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%

% random seed
%randseed=round(rand*1000000000); % The random seed is authomatically changed each run
randn('seed',randseed);
rand('seed',randseed);

% checks if the initial wealth is larger than zero. 
initialmin=init_wealth_avg/N;

% check if the simulation requires the use of savings. if yes it requires by definition also the income 
if risparmio==1; reddito=1; end;

T=zeros(N,tmax);
Subsidy=zeros(N,tmax);

% Initializations
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
    %Y(:,1)=income_avg; % agent's initial periodic income is the average above 
    Y=income_avg.*rand(N,tmax); % agent's saving rate 'rand' is applied to each period j
    YT=sum(Y,2); % sum7 is the sume across time of Y 
    if risparmio==1
        WY(:,1)=Y(:,1);
    end
end

% Creation of the return vector
r=max(-1,random(pd,[N,tmax]));
if reduction_of_interest=='t'
    tt=zeros(N,tmax);
    for i=1:N
        tt(i,:)=[1:tmax];
    end
    r=r./(log(1+tt)); 
end


if type_of_interest=='c'  % compounded interest
    for t=2:tmax % for each period t
        if reduction_of_interest=='w'
            if t>2
                r(:,t-1)=r(:,t-1)./log(1+sum(W(:,t-1),1));
            end
        end
        W(:,t)=max(0,W(:,t-1).*((1+r(:,t-1))));
        W(W(:,t-1)<=0,t)=0;
        if risparmio==1
            WY(:,t)= Y(:,t) + max(0,WY(:,t-1).*((1+r(:,t-1)))); % is it right, i was not summing for the W(i,t-1)<=0
        end
        if strcmp(taxes_yes,'y')
            [W,MeanTaxRate(t),MedianTaxRate(t),MeanRedistributionRate(t),MedianRedistributionRate(t),Taxes(:,t),Subsidy(:,t)]=taxes(W,t,N,tax_base,tax_type,redistribution_type,taxrate);
        end
    end
end

if type_of_interest=='s' % simple interest
    if risparmio==1
        sum3=zeros(N,1);
    end
    for t=2:tmax
        
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
             clc
             t
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
            [W,MeanTaxRate(t),MedianTaxRate(t),MeanRedistributionRate(t),MedianRedistributionRate(t),Taxes(:,t),Subsidy(:,t)]=taxes(W,t,N,tax_base,tax_type,redistribution_type,taxrate);
         end
    end
    if risparmio==1
        for t=1:tmax
            W(:,t)=W(:,t)+WY(:,t);
        end
    end
end

%%%%%
% creation of the output variables


if reddito==1 && risparmio==0
    W=W+Y;     % add income as and additive variable
end

for t=2:tmax
    growth(t)=(sum(W(:,t),1)-sum(W(:,t-1),1))/sum(W(:,t-1),1);
end

if reddito==1 && risparmio==0
    share_WY=sum(Y,1)./(sum(W,1)+sum(Y,1));
else
    if risparmio==1
        share_WY=sum(WY,1)./(sum(WY,1)+sum(W,1));
    else
        share_WY=NaN;
    end
end
    



M_growth=mean(growth);
Std_growth=std(growth);




[gini_avg]=ginicoeff(W',2,1);
theil_avg=zeros(tmax,1);
for t=1:tmax
    data=W(:,t);
    val=(data./mean(data)).*(log(data./mean(data)));
    num=sum(val);
    den=N*log(N);
    theil_avg(t)=num/den;
end

gini_tmax=gini_avg(end);
theil_tmax=theil_avg(end);
[max_theil,when_max_theil]=max(theil_avg);
[max_gini,when_max_gini]=max(gini_avg);
	

A_sf=sort(W(:,tmax), 'descend');
B_sf=[1:1:N];
C_sf=A_sf(1:round(N/10));
D_sf=[1:1:round(N/10)];
P_logw_logr_top=polyfit(log(D_sf),log(C_sf'),1); % aggiungi top, cambia lettere
Y_top=polyval(P_logw_logr_top,log(D_sf)); %aggiungi top x2 cambia lettere

   
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
C=sort(W(:,end), 'descend');
wealthquartiles(5,4)=sum(C(1:round(N/4)))/sum(C); %top25
wealthquartiles(5,3)=sum(C(round(N/4)+1:round(N/4)*2))/sum(C); %second25
wealthquartiles(5,2)=sum(C((round(N/4)*2)+1:round(N/4)*3))/sum(C); % third25
wealthquartiles(5,1)=sum(C(round(N/4)*3:end))/sum(C); % low25


tint=[1, 10, 100, 1000, tmax];
wealthdecile=zeros(length(tint),10);
for j=1:length(tint)
    C=sort(W(:,tint(j)), 'ascend');
    for ii=0:9
        wealthdecile(j,ii+1)=sum(C((round(N/10)*ii)+1:(round(N/10))*(ii+1)))/sum(C);
    end
end


% Plot of transition between deciles.
Weighted_changes_mtx_dec_w=zeros(10,10,tmax/time_steps_length);
for i=1:tmax-1 % for every time step
  [~, B]=sort(W(:,i), 'descend');
  [WW, E]=sort(W(:,i+1), 'descend');
  DecW=zeros(10,1); % create the deciles averages
  Declength=N/10;
  for q=1:10
      DecW(q)=median(WW(((q-1)*Declength+1):(q*Declength),1));
  end
  for j=1:N
      dec_init=floor(((B(j)-1)/N)*10)+1;
      if dec_init>10
          dec_init=10;
      end
      dec_end=floor(((E(j)-1)/N)*10)+1;
      if dec_end>10
          dec_end=10;
      end
      Weighted_changes_mtx_dec_w(dec_init,dec_end,i)=Weighted_changes_mtx_dec_w(dec_init,dec_end,i)+(abs(DecW(dec_init)-DecW(dec_end)))/(DecW(1)-DecW(10));
  end
  
end
Weighted_changes_mtx_dec_w=Weighted_changes_mtx_dec_w./N;

Varricchezzeabs=zeros(tmax,1);
Varricchezze=zeros(tmax,1);
for t=2:tmax
    Varricchezzeabs(t)=sum(abs((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t)))))./N;
    Varricchezze(t)=sum((W(:,t)-W(:,t-1))./(max(W(:,t))-min(W(:,t))))./N;
end
[x,y]=size(Weighted_changes_mtx_dec_w(:,:,1));
T=(tmax-1);
Weighted_movement_W=zeros(T,1);
for t=1:T
    for i=1:x
        for j=1:y
            if i~=j % movimenti 
                Weighted_movement_W(t)=Weighted_movement_W(t)+Weighted_changes_mtx_dec_w(i,j,t);
            end
        end
    end
end
weighted_movements_tmax=Weighted_movement_W(end);
window_of_mav=10;
mav_Weighted_movement_W=zeros(tmax-window_of_mav-1,1);
for i=1:tmax-window_of_mav-1 % construct the moving averages
   mav_Weighted_movement_W(i)=sum(Weighted_movement_W(i:i+window_of_mav-1))./((window_of_mav));   
end
Final_wealth=W(:,tmax)';

P_logw_logr_top=P_logw_logr_top(1);

Proportion_total_wealth=zeros(tmax,1); % proportion of total wealth owned by top 1%
Proportion_relative_wealth=zeros(tmax,1); % proportion of total wealth owned by top 1% wrt to bottom 50%
for t=1:tmax
    [Wlist, ~]=sort(W(:,t), 'descend');
    NUM=sum(sum(Wlist(1:50)));
    Wbottom=sum(sum(Wlist(2500:5000)));
    WTot=sum(sum(W(:,t))); 
    Proportion_total_wealth(t)=NUM./WTot;
    Proportion_relative_wealth(t)=NUM./Wbottom;
end


ratio_top1T_bottom10W=NaN;
ratio_top10T_bottom50W=NaN;
if taxes_yes=='y'
    create_ratio_Tax_Wealth
end


WatZero_end=sum(squeeze(W(:,tmax))<=0)

ratio_top50T_bottom50W=NaN;

