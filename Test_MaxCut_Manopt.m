maxcut_data = 'G67';


rng(0,'twister');
addpath utils;
addpath solver;

%% Load data

A = readSMAT(['~/Gset/',maxcut_data, '.smat']);

n = size(A,1);
C = spdiags(A*ones(n,1),0,n,n) - A;
C = 0.5*(C+C'); % symmetrize if not symmetric

clearvars A;

r = 10;

maxcut(C, r)

%% Last edit: Alp Yurtsever - July 24, 2020