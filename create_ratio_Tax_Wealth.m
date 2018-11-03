% This file pre-treats data for the visualizations concerning taxation and
% redistribution. It cannot be run alone.



ntop=round(N/100);
ntop10=round((N/100)*10);
nbottom=round((N/100)*10);
nbottom50=round((N/100)*50);
ratio_top1T_bottom10W=zeros(tmax,1);
ratio_top10T_bottom50W=zeros(tmax,1);
ratio_top50T_bottom50W=zeros(tmax,1);
for t=1:tmax
    [~, Pos]=sort(W(:,t), 'descend'); % ordina agenti in base a ricchezza al tempo timportant
    Pos1top=Pos(1:ntop); % e si segna le posizioni nel vettore degli agenti (non ordinato)
    Pos10top=Pos(1:ntop10);
    Pos50top=Pos(1:nbottom50-1);
    Tax_top1=Taxes(Pos1top,t); % mi segno le tasse pagate dal top 1%
    Tax_top10=Taxes(Pos10top,t);
    Tax_top50=Taxes(Pos50top,t);
    Pos10bottom=Pos(end-nbottom:end);
    Pos50bottom=Pos(end-nbottom50:end);
    W_bottom10=W(Pos10bottom,t);
    W_bottom50=W(Pos50bottom,t);
    ratio_top1T_bottom10W(t)=sum(Tax_top1)/sum(W_bottom10);
    ratio_top10T_bottom50W(t)=sum(Tax_top10)/sum(W_bottom50);
    ratio_top50T_bottom50W(t)=sum(Tax_top50)/sum(W_bottom50);
end
