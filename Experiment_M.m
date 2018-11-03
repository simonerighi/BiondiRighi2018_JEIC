% Creates the DATA of BiondiRighi (2018, JEIC) concerning the
% study of different modalities of taxation and redistribution.
% Visualizations for these data can be created through the file
% Visualizations_M_1.m and Visualizations_M_2.m


clear all
close all
clc

load('baseline_setup.mat');

filename='01Baseline';
filename=[filename '.mat'];


p=parpool;
share_WY=[];

parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i),share_WY_one,Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),MeanTaxRate(i),MedianTaxRate(i),MeanRedistributionRate(i),MedianRedistributionRate(i),...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,~,~,~,WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));

    share_WY=[share_WY share_WY_one'];

end
delete(p)

save(filename)

filename=[filename '.mat'];


%%%%%%
warning off
clear all
close all
clc

load('baseline_setup.mat');

filename='02TaxesProportionalLiberalist';
filename=[filename '.mat'];

taxes_yes='y'; 
tax_base='i'; % 'w' tassa su ricchezza (W_(t), 'i' tassa su income (W_t-W_(t-1))
tax_type='prop'; % 'prop': proportional 'prog': progessive,'lump': lumpsum
redistribution_type='lib'; % 'key': invesely proportional to tax base, 'lib': T/N per tutti
taxrate=0.05;

%matlabpool open %p=parpool;
share_WY=[];
%par

MeanTaxRate=zeros(tmax,Niter);
MedianTaxRate=zeros(tmax,Niter);
MeanRedistributionRate=zeros(tmax,Niter);
MedianRedistributionRate=zeros(tmax,Niter);

p=parpool
parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i),share_WY_one,Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),...
    MeanTaxRate_one,...
    MedianTaxRate_one,...
    MeanRedistributionRate_one,...
    MedianRedistributionRate_one,...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,ratio_top1T_bottom10W(:,i),ratio_top10T_bottom50W(:,i),ratio_top50T_bottom50W(:,i),WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));


    MeanTaxRate(:,i)=MeanTaxRate_one';
    MedianTaxRate(:,i)=MedianTaxRate_one';
    MeanRedistributionRate(:,i)=MeanRedistributionRate_one';
    MedianRedistributionRate(:,i)=MedianRedistributionRate_one';
    %clear MeanTaxRate_one
    share_WY=[share_WY share_WY_one'];

end
delete(p)
save(filename)

filename=[filename '.mat'];


%%%%
clear all
close all
clc

load('baseline_setup.mat');

filename='03TaxesProgressiveLiberalist';
filename=[filename '.mat'];

taxes_yes='y'; 
tax_base='i'; % 'w' tassa su ricchezza (W_(t), 'i' tassa su income (W_t-W_(t-1))
tax_type='prog'; % 'prop': proportional 'prog': progessive,'lump': lumpsum
redistribution_type='lib'; % 'key': invesely proportional to tax base, 'lib': T/N per tutti
taxrate=0.10;

%matlabpool open %p=parpool;
share_WY=[];
%par
MeanTaxRate=zeros(tmax,Niter);
MedianTaxRate=zeros(tmax,Niter);
MeanRedistributionRate=zeros(tmax,Niter);
MedianRedistributionRate=zeros(tmax,Niter);

p=parpool
parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i),share_WY_one,Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),...
    MeanTaxRate_one,...
    MedianTaxRate_one,...
    MeanRedistributionRate_one,...
    MedianRedistributionRate_one,...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,ratio_top1T_bottom10W(:,i),ratio_top10T_bottom50W(:,i),ratio_top50T_bottom50W(:,i),WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));


    MeanTaxRate(:,i)=MeanTaxRate_one';
    MedianTaxRate(:,i)=MedianTaxRate_one';
    MeanRedistributionRate(:,i)=MeanRedistributionRate_one';
    MedianRedistributionRate(:,i)=MedianRedistributionRate_one';
    %clear MeanTaxRate_one
    share_WY=[share_WY share_WY_one'];

end
delete(p)
save(filename)



%%%%
clear all
close all
clc

load('baseline_setup.mat');

filename='04TaxesProportionaKeynesian';
filename=[filename '.mat'];

taxes_yes='y'; 
tax_base='i'; % 'w' tassa su ricchezza (W_(t), 'i' tassa su income (W_t-W_(t-1))
tax_type='prop'; % 'prop': proportional 'prog': progessive,'lump': lumpsum
redistribution_type='key'; % 'key': invesely proportional to tax base, 'lib': T/N per tutti
taxrate=0.05;


%matlabpool open %p=parpool;
share_WY=[];
%par
MeanTaxRate=zeros(tmax,Niter);
MedianTaxRate=zeros(tmax,Niter);
MeanRedistributionRate=zeros(tmax,Niter);
MedianRedistributionRate=zeros(tmax,Niter);

p=parpool
parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i),share_WY_one,Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),...
    MeanTaxRate_one,...
    MedianTaxRate_one,...
    MeanRedistributionRate_one,...
    MedianRedistributionRate_one,...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,ratio_top1T_bottom10W(:,i),ratio_top10T_bottom50W(:,i),ratio_top50T_bottom50W(:,i),WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));


    MeanTaxRate(:,i)=MeanTaxRate_one';
    MedianTaxRate(:,i)=MedianTaxRate_one';
    MeanRedistributionRate(:,i)=MeanRedistributionRate_one';
    MedianRedistributionRate(:,i)=MedianRedistributionRate_one';
    %clear MeanTaxRate_one
    share_WY=[share_WY share_WY_one'];

end
delete(p)
save(filename)

filename=[filename '.mat'];


%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='05TaxesProgressiveKeynesian';
filename=[filename '.mat'];

taxes_yes='y'; 
tax_base='i'; % 'w' tassa su ricchezza (W_(t), 'i' tassa su income (W_t-W_(t-1))
tax_type='prog'; % 'prop': proportional 'prog': progessive,'lump': lumpsum
redistribution_type='key'; % 'key': invesely proportional to tax base, 'lib': T/N per tutti
taxrate=0.10;

%matlabpool open %p=parpool;
share_WY=[];
%par
MeanTaxRate=zeros(tmax,Niter);
MedianTaxRate=zeros(tmax,Niter);
MeanRedistributionRate=zeros(tmax,Niter);
MedianRedistributionRate=zeros(tmax,Niter);

p=parpool
parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i),share_WY_one,Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),...
    MeanTaxRate_one,...
    MedianTaxRate_one,...
    MeanRedistributionRate_one,...
    MedianRedistributionRate_one,...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,ratio_top1T_bottom10W(:,i),ratio_top10T_bottom50W(:,i),ratio_top50T_bottom50W(:,i),WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));


    MeanTaxRate(:,i)=MeanTaxRate_one';
    MedianTaxRate(:,i)=MedianTaxRate_one';
    MeanRedistributionRate(:,i)=MeanRedistributionRate_one';
    MedianRedistributionRate(:,i)=MedianRedistributionRate_one';
    %clear MeanTaxRate_one
    share_WY=[share_WY share_WY_one'];

end
delete(p)
save(filename)



name_experiment='Experiment_M/';
mkdir(name_experiment);
%visualizzazioni_M
%Visualizzazioni_nuove_M

%movefile('*.eps',name_experiment)
%movefile('*.fig',name_experiment)

movefile('baseline_setup.mat','baseline_setup.cip') % rinomino il creatore di basline.
movefile('*.mat',name_experiment)
movefile('baseline_setup.cip','baseline_setup.mat') % rinomino il creatore di basline.
copyfile('baseline_setup.mat',name_experiment)