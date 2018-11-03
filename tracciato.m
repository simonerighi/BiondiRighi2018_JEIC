% traks the wealth of the 20 richest individuals at time t

figure(12)
subplot(1,1,1)
load('')
ntop=20;
timportant=2000;
[~, Pos]=sort(W(:,timportant), 'descend');
Pos=Pos(1:ntop);
Pos_i_at_t=zeros(ntop,tmax);
for t=1:tmax
    [~, PositionInOriginalVct_at_t]=sort(W(:,t), 'descend');
    for i=1:ntop
        Pos_i_at_t(i,t)=find(PositionInOriginalVct_at_t==Pos(i));
    end
end
for i=1:ntop
    plot(Pos_i_at_t(i,:),'b-');
    hold on;
end
title(['Wealth Rank of top 20 richest agents at t=' num2str(timportant)])
xlabel('time')
ylabel('Ranking Wealth')





%%%%%

Proportion_wealth=zeros(tmax);
figure(13)
for t=1:tmax
    [Wlist, ~]=sort(W(:,t), 'descend');
    Wlist=Wlist(1:50);
    WTot=sum(sum(W(:,t)));
    Proportion_wealth(t)=Wlist/WTot;
end