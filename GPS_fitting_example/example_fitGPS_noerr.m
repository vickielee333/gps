%% Example to fit GPS data without considering error
% function used: build_G.m
% Vickie Lee @ VT
% 2021/03/15 v1

clear; clc;

file = 'LNCH.neu';
[t, N, E, U, Sn, Se, Su] = textread(file,'%f%f%f%f%f%f%f');
ENU = [E,N,U];
Senu = [Se,Sn,Su];

%% Detect bad data (large sigma & sigma = 0)
ic = find( Senu(:,1)>25 | Senu(:,2)>25 | Senu(:,3)>50 | Senu(:,1)==0 | Senu(:,2)==0 | Senu(:,3)==0);
t(ic,:) = []; ENU(ic,:) = []; Senu(ic,:) = [];

%% Make G matrix
G = build_G(t, 1, 1, 1, [2016.0979], [2016.0979], [1], [0], []);

%% Regression

for i = 1:3
    d = ENU(:,1);          

    m = (G'*G)\G'*d;       % regression
    dhat(:,i) = G*m;       % forwarding
    r(:,i) = d-dhat(:,i);  % residual
end

%% PLOT

figure; hold on;
subplot(311); hold on; grid on; box on;
plot(t,ENU(:,1),'k.')
plot(t,dhat(:,1),'r-')

subplot(312); hold on; grid on; box on;
plot(t,ENU(:,2),'k.')
plot(t,dhat(:,2),'r-')

subplot(313); hold on; grid on; box on;
plot(t,ENU(:,3),'k.')
plot(t,dhat(:,3),'r-')
