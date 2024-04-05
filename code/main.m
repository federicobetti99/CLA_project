%% clean
clear
close all
clc

%% import utilities and fix seed
addpath("utils")
rng("default")
rng(0)

% Running the three subsections below reproduces the plots shown in
% the report in the same order as they appear in the latter

%% run MC estimator
mc

%% run Lanczos estimator
lanczos

%% run Lanczos-MC estimator
lanczos_mc