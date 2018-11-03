% Creates the DATA of BiondiRighi (2018, JEIC) concerning the
% study of multiplicative wealth accumulation process (baseline) under different setup concerning the distribution of returns. 
% Visualizations for these data can be created through the file Visualizations_AandC.m

clear all
close all
clc

%baseline_setup_creator

load('baseline_setup.mat');

filename='01Baseline';
filename=[filename '.mat'];

interesse_1=0.05; %interesse di riferimento  (media o a) 0.25
interesse_2=0.05; %interesse di riferimento (sigma o b) 0.2


p=parpool
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

clear all
close all
clc

load('baseline_setup.mat');

filename='02mu0025sigma005';
filename=[filename '.mat'];

interesse_1=0.025; %interesse di riferimento  (media o a) 0.25
interesse_2=0.05; %interesse di riferimento (sigma o b) 0.2

p=parpool
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


%%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='03mu0075sigma005';
filename=[filename '.mat'];

interesse_1=0.075; %interesse di riferimento  (media o a) 0.25
interesse_2=0.05; %interesse di riferimento (sigma o b) 0.2

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


clear all
close all
clc

load('baseline_setup.mat');

filename='04mu005sigma0025';
filename=[filename '.mat'];

interesse_1=0.05; %interesse di riferimento  (media o a) 0.25
interesse_2=0.025; %interesse di riferimento (sigma o b) 0.2

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


%%%%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='05mu005sigma0075';
filename=[filename '.mat'];

interesse_1=0.05; %interesse di riferimento  (media o a) 0.25
interesse_2=0.075; %interesse di riferimento (sigma o b) 0.2


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

%%%%%

clear all
close all
clc

load('baseline_setup.mat');

filename='06a02b025Gamma';
filename=[filename '.mat'];

tipo_interesse='Gamma';
interesse_1=0.25; %interesse di riferimento  (media o a) 0.25
interesse_2=0.20; %interesse di riferimento (sigma o b) 0.2


p=parpool
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

mkdir('Experiment_AandC')
%visualizzazioniAandC

%additional_run_for_baseline

%movefile('*.eps','Experiment_AandC')
%smovefile('*.fig','Experiment_AandC')

movefile('baseline_setup.mat','baseline_setup.cip') % rinomino il creatore di basline.
movefile('*.mat','Experiment_AandC')
movefile('baseline_setup.cip','baseline_setup.mat') % rinomino il creatore di basline.
copyfile('baseline_setup.mat','Experiment_AandC')


