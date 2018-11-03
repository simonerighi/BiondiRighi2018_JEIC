function [W,MeanTaxRate,MedianTaxRate,MeanRedistributionRate,MedianRedistributionRate,T,Subsidy]=taxes(W,t,N,tax_WR,tax_type,redistribution_type,taxrate)
% this function implements the various modes of taxation and redistribution
% of wealth discussed in Biondi and Righi (2018, JEIC).

MeanTaxRate=-1;
MedianTaxRate=-1;
MeanRedistributionRate=-1;
MedianRedistributionRate=-1;

Vct_wealths=W(:,t); % I estract the vct of wealths
if tax_WR=='w'
    Tax_base=Vct_wealths;
else
    Tax_base=max(0,W(:,t)-W(:,t-1)); % se la tassa e' sul reddito, quelli con reddito negativo non vengono tassati :)
end

% vector of taxes
if strcmp(tax_type,'lump')==1 % takes a certain (taxrate) percentage of the average initial income
    T=mean(W(:,1))*taxrate;
    T((Vct_wealths-T)<0)=Vct_wealths((Vct_wealths-T)<0);
end
if strcmp(tax_type,'prop')==1
    T=Tax_base*taxrate; % takes a certain proportion of the current tax base
    MeanTaxRate=taxrate;
    MedianTaxRate=taxrate;
    
end
if strcmp(tax_type,'prog')==1 % tax progressive
    %TR=zeros(N,1);
    %for i=1:N
    TR=taxrate.*((Tax_base-min(Tax_base))/(max(Tax_base)-min(Tax_base)));
    T=Tax_base.*TR;
    MeanTaxRate=mean(TR);
    MedianTaxRate=median(TR);
    %end
end
Total_taxes=sum(T); % total taxes gathered by the state
if strcmp(redistribution_type,'lib')==1
    Subsidy=Total_taxes/N;
    MeanRedistributionRate=1/N;
    MedianRedistributionRate=1/N;
end
if strcmp(redistribution_type,'key')==1
    alpha=1-(Tax_base/sum(Tax_base));
    index=alpha./sum(alpha);
    Subsidy=Total_taxes.*index;
    MeanRedistributionRate=mean(index);
    MedianRedistributionRate=median(index);
end

Vct_wealths=Vct_wealths-T+Subsidy;
size(Vct_wealths);
W(:,t)=Vct_wealths; % the new we


