% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
% Requirement: Dynare 4.4.3
% Add to path the following folders and subfolders: GMMtoolbox, HOSScriptFiles and HOSUtils
addpath(genpath(pwd))

%% Replication of table 1
dynare usmodel_Gaussian
dynare usmodel_Student_t

%% Replication of table 2
dynare AnSchorfheide_Gaussian
dynare AnSchorfheide_Student_t

%% Replication of table 3
dynare neoclassical_Gaussian
dynare neoclassical_Student_t

%% GMM estimation example
% For GMM estimation run
dynare RBCmodel
% The used seeds for the Monte Carlo Simulation study are saved in seedsNR.mat
