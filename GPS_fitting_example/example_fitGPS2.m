%% Example to fit GPS data considering error using 2 iterations
% 1. 2 iterations to remove the bad data (residual > 3 rms)
% 2. plot the data with error
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

%% Regression twice

 % first iteration of regression
for i = 1:3
    d = ENU(:,i);          
    s = Senu(:,i);
    Cd = diag(s.^2);       % covariance matrix

    G = build_G(t, 1, 1, 1, [2016.0979], [2016.0979], [1], [0], []);
    
    m = (G'/Cd*G)\G'/Cd*d; % regression with covariance matrix
    dhat(:,i) = G*m;            % forwarding
    r(:,i) = d-dhat(:,i);  % residual
end
    
% Detect bad data (r > 3*wrms)
wrms = sqrt( sum(r.^2)/(length(t)) );

[id1,id2] = find( r>3*wrms);
t(id1,:) = []; ENU(id1,:) = []; Senu(id1,:) = [];
    
 % second iteration of regression
for i = 1:3    
   
    d = ENU(:,i);          
    s = Senu(:,i);
    Cd = diag(s.^2);       % covariance matrix
    
    G = build_G(t, 1, 1, 1, [2016.0979], [2016.0979], [1], [0], []);
    
    m = (G'/Cd*G)\G'/Cd*d;
    dhat2(:,i) = G*m;
    r2(:,i) = d-dhat2(:,i);
    
end


%% PLOT

figure; hold on;
subplot(311); hold on; grid on; box on;
plot([t t]',[ENU(:,1)+2*Senu(:,1), ENU(:,1)-2*Senu(:,1)]','color',[.5 .5 .5])
plot(t,ENU(:,1),'k.')
plot(t,dhat2(:,1),'r-')

subplot(312); hold on; grid on; box on;
plot([t t]',[ENU(:,2)+2*Senu(:,2), ENU(:,2)-2*Senu(:,2)]','color',[.5 .5 .5])
plot(t,ENU(:,2),'k.')
plot(t,dhat2(:,2),'r-')

subplot(313); hold on; grid on; box on;
plot([t t]',[ENU(:,3)+2*Senu(:,3), ENU(:,3)-2*Senu(:,3)]','color',[.5 .5 .5])
plot(t,ENU(:,3),'k.')
plot(t,dhat2(:,3),'r-')