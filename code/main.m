%% clean
clear
close all
clc

% Running this script reproduces the plots shown in the report in the
% same order as they appear in the latter (for the chosen matrix name)

%% import utilities and fix seed
addpath("utils")
rng("default")
rng(0)

%% set matrix name and load matrix
matname = "mesh3em5";
matfile = sprintf("../matrices/%s.mat", matname);
datastruct = load(matfile);
M = datastruct.Problem.A;

%% MC estimator
N = 10000;    % number of Monte Carlo samples
avgs = 10;    % number of averages
alpha = 0.05; % width of confidence interval
namefile = sprintf("../figures/%s/mc_estimator", matname);
mc(M, N, avgs, alpha, namefile);

%% Lanczos estimator
d = dictionary(["nos3", "mesh3em5", "mhdb416"], [500, 200, 400]);
k = lookup(d, matname); % number of Lanczos iterations
namefile = sprintf("../figures/%s/lanczos_estimator", matname);
lanczos(M, k, namefile);

%% Lanczos-MC estimator
N = 1000;                   % number of Monte Carlo samples
base_ks = [10; 50; 100; 200];
d = dictionary(["nos3", "mesh3em5", "mhdb416"], [{[base_ks; 500]}, {[base_ks; 250]}, {[base_ks; 400]}]);
ks = cell2mat(lookup(d, matname));    % number of Lanczos iterations
avgs = 10;                  % number of averages
alpha = 0.05;               % width of confidence interval
namefile = sprintf("../figures/%s/lanczos_mc_estimator", matname);
lanczos_mc(M, N, ks, avgs, alpha, namefile);

%% comparison
k = 100;        % number of Lanczos iterations
N = 1000;       % number of Monte Carlo samples
avgs = 10;      % number of averages
alpha = 0.05;   % width of confidence interval
namefile = sprintf("../figures/%s/comparison", matname);
comparison(M, k, N, avgs, alpha, namefile);