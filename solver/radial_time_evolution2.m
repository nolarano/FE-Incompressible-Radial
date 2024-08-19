function [err,R,max_hoop] = radial_time_evolution2(param,param_rhoc)
global NewDir

if nargin>2, disp('Input should be a vector of parameters.'); return; end
if nargin>0 && ~isempty(param), radial_set_parameters(param); end

load([NewDir '/parameters.mat']);

[r,fer,fet,Fr,Ft,Gr,Gt,p,v,R,radial,hoop,VT,C,residual,Gam] = solve_radial_final2023(0);

err = get_error(R,f0);
max_hoop = max(max(abs(hoop)));

%save(fullfile(NewDir,'solution.mat'), 'R','err')
save(fullfile(NewDir,'solution.mat'), 'r','C','R','Gam','radial','hoop','v','p','err','fer','fet','Fr','Ft','Gr','Gt')

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
