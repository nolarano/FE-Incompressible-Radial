function radial_set_parameters(param)
global NewDir
beta = param(1);
k = param(2);
gammac = param(3);
L = param(4);
c_H = param(5);
pBar = param(6);
eta = param(7);
R0 = param(8);
tumorID = param(9);
TimeRelease = param(10);
delB = param(11); 
sen = param(12); 


scale_v = 0; %alpha
mu=1; %elastic coefficient on stresses

Nr = 1001; dt = (1e-3);
record_every = (1e3); % record every this many time steps 

%Compare with the experimental data
switch tumorID
    case -2 % parameter study of Morgan. load Free data regardless of pBar
        load morgan_fig1c_free1s.mat f0
        f0 = f0(2:end,:);
        tspan = [f0(1,1) 100];
    case -1
        isFree = pBar==0;
        if isFree
            load morgan_fig1c_free1s.mat f0
            f0 = f0(2:end,:);
        else
            load morgan_fig1c_pressures.mat f0
            f0 = f0(2:end,:);
        end
%         tspan = [0 1];
        tspan = [f0(1,1) 100];
    case 0
        load helminger_fig1a.mat f0
%         tspan = [f0(1,1) 200];
        tspan = [f0(1,1) 100];
    case {1,7}
        load helminger_fig1a.mat f07
        f0 = f07;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 2
        load helminger_fig1a.mat f10
        f0 = f10;
        tspan = [f0(1,1) 200];
    case 3
        load helminger_fig1a.mat f03
        f0 = f03(1:end,:);
        tspan = [f0(1,1) 200];
    case 5
        load helminger_fig1a.mat f05
        f0 = f05;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 8
        load helminger_fig1a.mat f08
        f0 = f08;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 9
        load helminger_fig1a.mat f09
        f0 = f09;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 10
        load helminger_fig1b.mat f0
        f0 = f0(2:end,:);
        tspan = [f0(1,1) 150];
    case 11
        load helminger_fig1b.mat f07
        f0 = f07(3:end,:);
        tspan = [f0(1,1) 150];
    case 12
        load helminger_fig1b.mat f10
        f0 = f10;
        tspan = [f0(1,1) 150];
    case 20
        load prl_fig1as.mat f0
        tspan = [f0(1,1) 100];
    case 21 % Morgan pressure release
        load prl_fig1a_500pa_releases.mat f0
        f0 = f0(2:end,:); % skip first data point
        tspan = [f0(1,1) 100];
    case 221 % Morgan pressure release, pre-release data
        load prl_fig1a_500pa_releases.mat f0
        f0 = f0(2:end-1,:);
        tspan = [f0(1,1) 100];
    case 31
        load prl_fig1b.mat f0
        f0(:,2) = f0(:,2)/2;
    case 32 
        load prl_fig1b_500pa.mat f0
        f0(:,2) = f0(:,2)/2;
    case 33 
        load prl_fig1b_2000pa.mat f0
        f0(:,2) = f0(:,2)/2;
    case 34 
        load prl_fig1b_5000pa.mat f0
        f0(:,2) = f0(:,2)/2;

end

if exist('radial_grid_search.lock','file')
    tspan = [f0(1,1) f0(end,1)];
end

R0 = f0(1,2);
%R0 =param(5);
tspan = [0 2];


save(fullfile(NewDir,'parameters.mat'), ...
    'param','f0','R0','dt','tspan','Nr','scale_v','mu','tumorID',...
    'record_every');

end
