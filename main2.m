function main2(beta,k,gammac,L,ch,eta,R0,type)%,delB,sen)
global NewDir
addpath('solver')
addpath('data_as_diameter')
addpath('plotting')

%% parameters for Fig 1
param = [
    beta;    % beta
    k;      % k chemical coefficient
    gammac;      %gammac
    L;      % L
    0;      % c_H
    ch;      % hydrostatic pressure in surrounding material, \bar p
    eta;    % Eta in Gamma 
    R0;     % gamma_lambda eq. 8 (lambda_A in the code)
    type;      % tumor ID for experimental data, see radial_set_parameters, line 47
    0;      % Time Release
    0;      %delB;  
    0;      %sen;    
 ];

NewDir = ['./Fit',num2str(type),'P',num2str(param(6)),'eta',num2str(param(7)),...
          'k',num2str(param(2)),'gamc',num2str(param(3)),'L',num2str(param(4)),...
          'Beta',num2str(param(1)),'CH',num2str(param(5)),'R0',num2str(param(8))];

%,'release',num2str(param(11)),'DelB',num2str(param(11)),'Phi',num2str(param(12))];

% if ~isfile([NewDir,'/solution.mat'])
 mkdir(NewDir(3:end))
% else
%     return
% end

%% run simulation
radial_time_evolution2(param);

% solution.mat structure:
% 'R' is tumor radius
% solved field functions are Nr*Nt arrays. For example, p(:,1) is the
%   initial condition of pressure, p(end,:) is the boundary solution of
%   pressure for all times.
%   Note: 'radial' and 'hoop' are the total radial/hoop stress.

end