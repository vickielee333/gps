function [G] = build_G(time, rate_flag, yr_flag, syr_flag, ...
                       offsets, epoch_tau, TAUs, k_value, VBREAKs)

% Input:
%     time        nx1 vector of x axis values (time for instance)
%     rate_flag   [0 1] estimate a rate
%     yr_flag     [0 1] estimate a yearly signal (phase wrt Jan.1)
%                (adds 2 parameters to the model)
%     syr_flag    [0 1] estimate a semi-annual signal (phase wrt Jan.1)
%                (adds 2 parameters to the model)
%     offsets     (nx1) vector of offset times.
%     epoch_tau   (dx1) vector of epochs for exponential decays
%     TAUs        (dx1) vector of exponential decay time constants
%     k_value     (dx1) vector of the k value that control the curvature of
%                       exponential decays (k=1~7, k=0 for unchanged curve.)
%     VBREAKs     (vx1) vector of velocity break times.  A new rate 
%                will be estimated starting at each break.
%
% Output:
%     G = [(y-axis-1) (linear-1) (yr-2) (semi-yr-2) (offset-?) (poseismic-?)]

% build the G matrix
n=length(time);
G=[ones(n,1)];

if (rate_flag==1)
    G=[G, time];
    nvb=length(VBREAKs);
    if nvb>0
        mask=time<VBREAKs(1);
        G(~mask,2)=VBREAKs(1)*ones(sum(~mask),1);
        for i=1:nvb-1
            mask1=time<VBREAKs(i);
            mask2=time>VBREAKs(i+1);
            G=[G [zeros(sum(mask1),1);time(~mask1&(~mask2))-VBREAKs(i); ...
          ones(sum(mask2),1)*VBREAKs(i+1)-VBREAKs(i)]];
        end
        lastvb = VBREAKs(nvb);
        mask=time>lastvb;
        G=[G mask.*(time-lastvb)];    
    end
end

% remove periodic effect
if (yr_flag==1)
    G=[G, sin(2*pi*time), cos(2*pi*time)];
end

if (syr_flag==1)
    G=[G, sin(4*pi*time), cos(4*pi*time)];
end

% estimate magnitude of offsets
for i=1:length(offsets)
    G=[G, time>=offsets(i)];
end

% estimate magnitude of post-seismic deformation
for i=1:length(TAUs)
    if TAUs(i)~=0
        if k_value(i)~=0
        G=[G, (time>=epoch_tau(i)).*(1-(2/k_value(i)).*acoth(exp((time-epoch_tau(i))./TAUs(i)).*coth(k_value(i)/2)))];    
        else
        G=[G, (time>=epoch_tau(i)).*(1-exp(-(time-epoch_tau(i))./TAUs(i)))];
        end
    end
end


