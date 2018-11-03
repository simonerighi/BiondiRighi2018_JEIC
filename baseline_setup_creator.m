% This file creates the baseline_setup.mat file containing the value of
% alla variales in the most baseline case. Such file is then loaded and
% only the relevant variables are changed to obtain different configurations.
% Using this method ensures the same initial randseed sequence for all
% simulations

Niter=100;


randseed=round(rand(Niter,1)*1000000000);

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
tipo_interesse='Norma'; % 'Gamma'
interesse_1=0.05; %interesse di riferimento  (media o a) 0.25
interesse_2=0.05; %interesse di riferimento (sigma o b) 0.2
init_wealth_avg=10; % average of initial wealths
initial_wealth_type='equ'; % 'equ'
income_avg=1; % average of periodic income



save('baseline_setup.mat')
