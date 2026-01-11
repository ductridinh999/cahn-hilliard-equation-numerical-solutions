clear all; close all; clc;

% Parameters
N = 100;        % Grid points
T = 5000;       % Time steps
ep = 0.01;      % Interface width
dt = 1e-5;      % Time step
seed = 42;      % Random seed
k = 1;          % Random initial condition
visualize = true;

[cvecs, u_vecs, time] = CH1D_Coupled_Live(N, T, ep, dt, seed, k, visualize);
