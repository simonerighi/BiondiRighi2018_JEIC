% This file creates data and visualizations for figure 3 in Biondi Righi
% (2018, JEIC). It is called by Figure3_data.m and cannot be
% run alone.

[~, Pos]=sort(W(:,timportant), 'descend'); % ordina agenti in base a ricchezza al tempo timportant
Pos=Pos(1:ntop); % e si segna le posizioni nel vettore degli agenti (non ordinato)
Pos_i_at_t=zeros(ntop,1000); % azzera una matrice che traccia le posizioni da t a tmax)
avg_pos_cum=[];
tt=1;
for t=timportant+1:timportant+1000
    [~, PositionInOriginalVct_at_t]=sort(W(:,t), 'descend'); % ordina il vettore al tempo t
    for i=1:ntop
        Pos_i_at_t(i,tt)=find(PositionInOriginalVct_at_t==Pos(i)); % e si segna la posizione nella matrice di tracciamento
    end
    avg_pos_cum=[avg_pos_cum Pos_i_at_t(i,tt)]; % vettore con tutte le posizoni dei top agents.
    tt=tt+1;
end
% evolutione temporale
%for i=1:ntop
%    plot(timportant+1:timportant+1000,Pos_i_at_t(i,:),'b-');
%    hold on;
%end
%xlabel('time')
%ylabel('Ranking Wealth')

% media posizioni
avg_pos_i=zeros(ntop,1);

avg_pos_i=mean(Pos_i_at_t,2); % posizione media

avg_pos_t=mean(Pos_i_at_t,1); % posizione media at t
std_pos_t=std(Pos_i_at_t,0,1); % posizione std at t



