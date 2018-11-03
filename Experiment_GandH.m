% Creates the DATA of BiondiRighi (2018, JEIC) concerning the
% study of the effects of intrudcing Income and Savings in the wealth accumulation process. 
% Visualizations for these data can be created through the file Visualizations_GandH.m


clear all
close all
clc

load('baseline_setup.mat');

filename='01Baseline';
filename=[filename '.mat'];


p=parpool;

parfor i=1:Niter
    display([filename ' '  num2str(i)])
[gini_avg(:,i),theil_avg(:,i), mav_Weighted_movement_W(:,i), Varricchezzeabs(:,i),...
    Varricchezze(:,i), share_WY(i), Final_wealth(i,:),A_sf(:,i),B_sf(i,:),C_sf(:,i),D_sf(i,:), Y_top(i,:),...
    wealthquartiles(:,:,i), wealthdecile(:,:,i), P_logw_logr_top(i), M_growth(i), Std_growth(i),...
    gini_tmax(i),theil_tmax(i),when_max_gini(i), max_gini(i), when_max_theil(i), max_theil(i),...
    weighted_movements_tmax(i),MeanTaxRate(i),MedianTaxRate(i),MeanRedistributionRate(i),MedianRedistributionRate(i),...
    Proportion_total_wealth(:,i),Proportion_relative_wealth(:,i),~,~,~,~,WatZero_end(i)]=OurProcessGenericoFast_fct(reddito,risparmio,...
    reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
    taxrate,N,tmax,tipo_interesse,interesse_1,interesse_2,init_wealth_avg,...
    initial_wealth_type,income_avg,randseed(i));
end

delete(p)


save(filename)

filename=[filename '.mat'];


%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='02WithReddito';
filename=[filename '.mat'];

reddito=1;

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


%%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='03WithRedditoERisparmio';
filename=[filename '.mat'];

reddito=1;
risparmio=1;

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



name_experiment='Experiment_GandH/';
mkdir(name_experiment)
%visualizzazioni_GandH

%movefile('*.eps',name_experiment)
%movefile('*.fig',name_experiment)

movefile('baseline_setup.mat','baseline_setup.cip') % rinomino il creatore di baseline.
movefile('*.mat',name_experiment)
movefile('baseline_setup.cip','baseline_setup.mat') % rinomino il creatore di baseline.
copyfile('baseline_setup.mat',name_experiment)
