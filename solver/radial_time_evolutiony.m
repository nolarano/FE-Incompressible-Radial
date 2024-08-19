function [err,R,max_hoop] = radial_time_evolutiony(param,param_rhoc)
global NewDir
if nargin>2, disp('Input should be a vector of parameters.'); return; end
if nargin>0 && ~isempty(param), radial_set_parametersy(param); end


load([NewDir '/parameters.mat']);

%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,VT,YR,C,B,LA,LAMBDA,TMP,lamBs,VA,PRA,CA,PD,residual] = solve_radial_fvp_bugcheck2(0);

%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual] = solve_radial_fvp_old(0);
%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual,Gam] = solve_radial_2D_DisN(0);
%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual,Gam] = solve_radial_DisN3(0);
%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual,Gam] = solve_radial_DisN_23(0);
[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual,Gam] = solve_radial_y(0);
%[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,~,~,C,~,LA,LAMBDA,~,~,~,~,~,PD,residual,Gam] = solve_radial_2023_witholdc(0);
err = get_error(R,f0);
max_hoop = max(max(abs(hoop)));



if ~exist('radial_grid_search.lock','file')
%     save 'solution.mat' r y p v R radial hoop VT YR C B LA LAMBDA TMP lamBs VA PRA CA PD RHOC residual
%    save(fullfile(NewDir,'solution.mat'), 'r', 'fer', 'fet', 'Fr', 'Ft', 'Gr', 'Gt', 'p' ,'v', 'R', 'radial', 'hoop', 'C', 'LA', 'LAMBDA', 'PD', 'residual','Gam')
%   save(fullfile(NewDir,'solution.mat'), 'r', 'fer', 'fet', 'Fr', 'Ft', 'Gr', 'Gt', 'p' ,'v', 'R', 'radial', 'hoop', 'C', 'LA', 'LAMBDA', 'PD', 'residual','err','Gam')
   save(fullfile(NewDir,'solution.mat'), 'r','C','R','Gam','radial','hoop','v','p','err','fer','fet')
%   save(fullfile(NewDir,'solution.mat'), 'R','err')

%     disp(['Error = ' num2str(err)]);
    % radial_make_movies;
end

if exist(fullfile(NewDir,'solution.mat'),'file')
%    save(fullfile(NewDir,'solution.mat'),'v','fer','fet','Fr','Ft','Gr','Gt','r','R','hoop','err','C','-append')
    save(fullfile(NewDir,'solution.mat'), 'r','C','R','Gam','radial','hoop','v','p','err','fer','fet','-append')
%    save(fullfile(NewDir,'solution.mat'), 'R','err','-append')
else
%    save(fullfile(NewDir,'solution.mat'),'v','fer','fet','Fr','Ft','Gr','Gt','r','R','hoop','err','C')
    save(fullfile(NewDir,'solution.mat'), 'r','C','R','Gam','radial','hoop','v','p','err','fer','fet')
% save(fullfile(NewDir,'solution.mat'), 'R','err')
end

%  if exist('param_rhoc','var')
%      solve_rhoc_fast(param_rhoc);
%  else
%      solve_rhoc_fast;
%  end

%  if ~exist('radial_grid_search.lock','file')
%      radial_plot_single_run;
%  end
end

function err = get_error(R,f0)
global NewDir
load([NewDir,'/parameters.mat']);

err = 0;
for i=1:length(f0(:,1))
    weight = 1;% + (f0(i,1)>4) * 10; % last few points have more weights
    err = err + weight*(R(min((f0(i,1)-f0(1,1))/dt+1,length(R)))-f0(i,2))^2;
end
err = sqrt(err);
end
