% Creates the DATA of BiondiRighi (2018, JEIC) concerning the
% exploration of different average and std in the returns of a multiplicative wealth accumulation process. 
% Visualizations for these data can be created through the file Visualizations_AandC.m


clear all
close all
clc
 
load('baseline_setup.mat');

filename='Experiment_B';
filename=[filename '.mat'];

interesse_1=[0.02:0.01:0.08];
interesse_2=[0.02:0.01:0.08];
NIter=100;
randseed=round(rand(NIter,1)*1000000000);


for i=1:length(interesse_1)
    for j=1:length(interesse_2)
        clc
        display(num2str(interesse_1(i)))
        display(num2str(interesse_2(j)))
        int_1=interesse_1(i);
        int_2=interesse_2(j);
        for iter=1:NIter
            display(num2str(iter))
            [~,~, ~, Varricchezzeabs,...
            Varricchezze, ~, ~,~,~,~,~, ~, ~, ~,...
            P_logw_logr_top(i,j,iter), M_growth(i,j,iter), Std_growth(i,j,iter),...
            gini_tmax(i,j,iter),theil_tmax(i,j,iter),when_max_gini(i,j,iter), max_gini(i,j,iter),...
            when_max_theil(i,j,iter), max_theil(i,j,iter),...
            weighted_movements_tmax(i,j,iter),~,~,~,~,~,...
            ~,~,~,~,~,WatZero_end(i,j,iter)]=OurProcessGenericoFast_fct(reddito,risparmio,...
            reduction_of_interest,type_of_interest,taxes_yes,tax_base,tax_type,redistribution_type,...
            taxrate,N,tmax,tipo_interesse,int_1,int_2,init_wealth_avg,...
            initial_wealth_type,income_avg,randseed(iter));
           
            meanVarRicchezzeabs(i,j,iter)=mean(Varricchezzeabs);
            meanVarRicchezze(i,j,iter)=mean(Varricchezze);
            FinalVarRicchezzeabs(i,j,iter)=Varricchezzeabs(end);
            FinalVarRicchezze(i,j,iter)=Varricchezze(end);
            save(filename)
        end
        
    end
end



visualizzazioni_B

name_experiment='Experiment_B/';
mkdir(name_experiment);

movefile('*.eps',name_experiment)
movefile('*.fig',name_experiment)

movefile('baseline_setup.mat','baseline_setup.cip') % rinomino il creatore di basline.
movefile('*.mat',name_experiment)
movefile('baseline_setup.cip','baseline_setup.mat') % rinomino il creatore di basline.
copyfile('baseline_setup.mat',name_experiment)