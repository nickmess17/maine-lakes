%script to run GA on MyLake
clc,clear,close all
fun=@NJM_MyLake;
nvars=6;%number of parameters to calibrate
lb=[0.01 0.5 1 1 0.25 0.8];%lower bounds of parameters
ub=[0.03 0.99 2 1.5 5 1.3];%upper bounds of parameters
population_size = 72;  % Populations size for each generation of the genetic algorithm
max_generations = 15;  % How many generations to run the genetic algorithm for
parallelize     = true; % 15 generation takes 12 hours on 24 cores
options=optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);
[par,fval]=ga(fun,nvars,[],[],[],[],lb,ub,[],[],options);