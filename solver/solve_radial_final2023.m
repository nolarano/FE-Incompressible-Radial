%implicit nutrient equation
function [r,fe_r,fe_t,F_r,F_t,G_r,G_t,P,V,R,radial,hoop,VT,C,residual,Gam] = solve_radial_final2023(~)
global NewDir

if nargin==0, radial_time_evolution; return; end % call radial_time_evolution when directly executed

load([NewDir,'/parameters.mat'],'param','f0','R0','dt','tspan','Nr','scale_v','mu','tumorID',...
    'record_every');

R0_ref = R0;
dim = 3;
beta = param(1); k = param(2); gammac = param(3); L = param(4);
cH = param(5); pBar = param(6); eta = param(7); R0 = param(8);
tumorID = param(9); TimeRelease = param(10);delB = param(11);
sen = param(12);

if ~exist('record_every','var'), record_every = 1; end

% Helminger Fig 1b pressure release time. 0 to disable.
if ~ismember('match_fig1b', who('-file', [NewDir,'/parameters.mat']))
    match_fig1b = 0;
else
    load([NewDir,'/parameters.mat'],'match_fig1b');
    if match_fig1b>0 && ~exist('radial_grid_search.lock','file')
        %         disp(['Pressure release at T=' num2str(match_fig1b)]);
        R0_ref = R0;
    end
end


hoop_warned = false;
debug_output = 0;

% p, v iteration
use_nonlinear = 0;
max_vp_ite = 100; % max iterations of v,p
tol = 1e-6;
Drho = L^2;
fsolve_ops = optimoptions(@fsolve,'Display','off','algorithm','trust-region');

%% nonlinear solver
    function F = fvp(x)
        % x = [v; p]
        F = ones(size(x));

        % v
        gamman = get_Gamma(fer,fet,x(Nr+1:2*Nr));
        F(1) = x(1); % v(0) = 0
        F(2:Nr-1) = r(2:Nr-1).*x(3:Nr)+4.*dx.*x(2:Nr-1)-r(2:Nr-1).*x(1:Nr-2)...
            -2.*dx.*r(2:Nr-1).*R(n-1).*gamman(2:Nr-1);
        F(Nr) = ((2*dx+r(Nr)).*x(Nr) - r(Nr).*x(Nr-1) ...
            - dx*R(n-1)*r(Nr).*(gamman(Nr)));
        % p
        F(Nr+1:2*Nr-2) = -3/2*x(Nr+1:2*Nr-2)+4/2*x(Nr+2:2*Nr-1)-1/2*x(Nr+3:2*Nr) ...
            + dx*scale_v*R(n-1)*x(1:Nr-2) - dx*fr(1:Nr-2);
        F(2*Nr-1) = -1/2*x(2*Nr-2)+1/2*x(2*Nr) ...
            + dx*scale_v*R(n-1)*x(Nr-1) - dx*fr(Nr-1);
        F(2*Nr) = x(2*Nr) - pEnd; % p(1) = pEnd

    end

%% feedback functions
    function Gamma = get_Gamma(fer,fet,p)

        Ptil = -mu*(fer.^2+2*fet.^2)/dim+p;
        W = 1/2*mu*(fer.^2+2*fet.^2-dim);

        Gamma = eta*c.*(1/2*k*c.^2-(W+Ptil));
        %Gamma = lambda_base*c.*(1/2*c.^2-(W+Ptil)/k); FEEDBACK COEFFICIENT

        %Huaming's witout feedback
        %       Gamma = 1.1*c-0.2;

        %Huaming with feedback
        %         cir_sig = fet.^2-p;
        %         Gamma = 1.1 ...%./(1+0.1.*(cir_sig).^2.*(cir_sig < -tol)).*1 ...
        %          -0.2-0.9.*(0.1.*(cir_sig).^2.*(cir_sig < -tol))./(1+0.1.*(cir_sig).^2.*(cir_sig < -tol));

    end

%% initialize or restart
r = linspace(0,1,Nr)';
Nt = round((tspan(2)-tspan(1))/dt);
dx = r(2)-r(1);

R = zeros(Nt+1,1); R(1) = R0;
fer = ones(length(r),1);
fet = ones(length(r),1);
Fr = ones(length(r),1);
Ft = ones(length(r),1);
Gr = ones(length(r),1);
Gt = ones(length(r),1);
fe_r(:,1) = fer;
fe_t(:,1) = fet;
F_r(:,1) = Fr;
F_t(:,1) = Ft;
G_r(:,1) = Gr;
G_t(:,1) = Gt;
c = ones(length(r),1);
%c = (R(1)*sinh(r))./(r*sinh(R(1)));
C(:,1) = c;
Gam(:,1) = get_Gamma(ones(length(r),1),ones(length(r),1),0);
nStart = 2;
residual = zeros(Nt,1); residual(1) = NaN;


%% finite difference loop
for n=nStart:Nt+1
    currTime = tspan(1)+dt*(n-1);
    %if disp_progress && mod(currTime-tspan(1),disp_progress)==0, disp(['T=' num2str(currTime)]); end
    if (TimeRelease>0 && abs(currTime-TimeRelease)<dt) % Helminger Fig 1b, pressure release
        if cH>0
            cH = 0;
        end
        if pBar>0
            pBar = 0;
        end
    end

    %% calculate stress from lagged p. Note: fet is current
    %fe_rr = gradient(fer, dx);
    %fe_tr = gradient(fet, dx);
    fe_tr = [(-3*fet(1)+4*fet(2)-fet(3))/2/dx;(fet(3:end)-fet(1:end-2))/2/dx;(fet(end-2)-4*fet(end-1)+3*fet(end))/2/dx];

    if ~exist('p','var'), p = 0; end % assuming no feedback at T=0

        sr = mu*fet.^(-4);
        st = mu*fet.^(2);

    gamman = get_Gamma(fer,fet,p);

    %% update v and p togetheR
    gel_stress = -cH/2*(5-R0_ref*(R0_ref^3+4*R(n-1)^3)/(R(n-1)^4)); %...
        %- cH2*(1-(2*R(n-1)^3-R0_ref^3)/(R(n-1)^2*R0_ref));

    sigma_r3 = -4*mu*(fet).^(-5).*fe_tr;
    pEnd = mu*fet(end).^(-4)- gel_stress + pBar;

    fr = (sigma_r3 + 2./r.*(sr - st));
    fr(1) = 3*sigma_r3(1)+-4*mu*fet(1)*fe_tr(1);
    N = 2*Nr; % 2x large matrix
    diagV = [1; 4*dx.*ones(Nr-2,1); 2*dx+1]; subDiagV = -r(2:end); % central diff
    supDiagV = [0;r(2:end-1);0];
    %diagP = [1;zeros(Nr-2,1);1]; supDiagP = -ones(Nr-1,1); subDiagP = [0;ones(Nr-2,1);0]; % central difference
    diagP = [-3/2*ones(Nr-2,1);0;1]; supDiagP = [4/2*ones(Nr-2,1);1/2]; subDiagP = [0;zeros(Nr-3,1);-1/2;0]; supsupDiagP = -1/2*ones(Nr-2,1);% forward difference
    VPcoupled = [ones(Nr-1,1)*scale_v*R(n-1)*dx; 0]; % lower left submatrix. the p_r+R*v bit
    diag = [diagV; diagP];
    supDiag = [supDiagV; supDiagP];
    subDiag = [subDiagV; subDiagP];

    A = sparse(1:N,1:N,diag) + sparse(1:N-1,2:N,supDiag,N,N) + sparse(2:N,1:N-1,subDiag,N,N) ...
        + sparse(Nr+1:N,1:Nr,VPcoupled,N,N)+sparse(Nr+1:N-2,Nr+3:N,supsupDiagP,N,N);

    %     diagV = [1; -1+2*dx./r(2); 3+4*dx./r(3:end)];
    %     subDiagV = [-1;-4*ones(Nr-2,1)];
    %     subsubDiagV = ones(Nr-2,1);
    %     supDiagV = [0;0;zeros(Nr-2,1)];
    %     diagP = [-3*ones(Nr-2,1);-1;1];
    %     subDiagP = [0;zeros(Nr-3,1);0;0];
    %     supDiagP = [4*ones(Nr-2,1);1];
    %     supsupDiagP = -ones(Nr-2,1);
    %     VPcoupled = [2*ones(Nr-2,1)*scale_v*R(n-1)*dx;scale_v*R(n-1)*dx; 0]; % lower left submatrix. the p_r+R*v bit
    %     diag = [diagV; diagP];
    %     supDiag = [supDiagV; supDiagP];
    %     subDiag = [subDiagV; subDiagP];
    %     A = sparse(1:N,1:N,diag);
    %     A = A + sparse(1:N-1,2:N,supDiag,N,N);
    %     A = A + sparse(2:N,1:N-1,subDiag,N,N);
    %     A = A + sparse(Nr+1:N,1:Nr,VPcoupled,N,N);
    %     A = A + sparse(Nr+1:N-2,Nr+3:N,supsupDiagP,N,N);
    %     A = A+ sparse(3:Nr,1:Nr-2,subsubDiagV,N,N);

    if ~use_nonlinear
        for i_vp=1:max_vp_ite % iterate a few times to eliminate oscillation
            % RHS of v and p

            vr = R(n-1)*gamman.*r;
            rhsV = [0; 2*dx*vr(2:end-1); dx*vr(end)];
            rhsP = [dx*fr(1:Nr-2); dx*fr(Nr-1); pEnd];
            %rhsV = [0; 2*dx*vr(2); 2*dx*vr(3:end)];
            %rhsP = [2*dx*fr(1:Nr-2); dx*fr(Nr-1); pEnd];

            rhs = [rhsV; rhsP];
            % solution is [v; p]
            sol = A\rhs;
            v = sol(1:Nr);
            p = sol(Nr+1:end);
            gamman = get_Gamma(fer,fet,p);

            res = norm(fvp([v; p]));
            if res<tol
                break
            else
                if i_vp==max_vp_ite, use_nonlinear = 1; end
            end
        end
    else
        x0 = [v; p];
        sol = fsolve(@fvp, x0, fsolve_ops);
        v = sol(1:Nr);
        p = sol((Nr+1):(2*Nr));
        gamman = get_Gamma(fer,fet,p);
    end

    residual(n) = norm(fvp([v; p]));
    if debug_output
        fprintf('T=%.2f, residual=%e, toc=%f\n',currTime,residual(n),toc);
    end

    % anatical solution of v, p_r
    va = R(n-1)./r.^2.*cumsimps(r,(gamman).*r.^2);
    pra = -scale_v*R(n-1)*va + R(n-1)*sigma_r3 + 2./r.*(sr - st);

    % pD is the pressure from Darcy's law
    pD = zeros(Nr,1);
    pD(end) = 2*gamman(end)/R(n-1) - gel_stress + pBar;
    pDrhs = -scale_v*R(n-1)*va;
    for i_pD=Nr-1:-1:1
        pD(i_pD) = pD(i_pD+1) - dx*pDrhs(i_pD);
    end

    %% update R
    R(n) = R(n-1) + dt*v(end);
    if R(n) < tol
        break;
    end

    %% solve y
    drdt = v(end);
    vt = (v-r*drdt)/R(n);

    v_r = zeros(Nr,1);
    v_r(1) = gamman(1)/3;
    v_r(2:end) = gamman(2:end)-v(2:end).*(2./(r(2:end).*R(n-1)));

    %Siginv = (2*(st-p)+(sr-p))/3;
    %b = beta * ones(size(r)).*(1+delB*sen.*Siginv./(1+sen.*Siginv)); %beta remodeling. 
    
    b = beta * ones(size(r));
    bbar = 0 * ones(size(r));
    fetnew = zeros(Nr,1);
    fetnew(1)   = fet(1);% + dt.*fet(1) .*(-1/3*b(1).*(fet(1)^2-fet(1)^-4));
    %   fetnew(end) = fet(end)+ dt.*fet(end).*(v(end)./(r(end).*R(n-1))-gamman(end)/3-1/3*b(end).*(fet(end)^2-fet(end)^-4));

    %setup for nutrient solver
    upwc = zeros(Nr,1);
    %upwc(1) = 0;
    %upwc(end) = 0;

    % upwind
    vp = max(vt, 0); vn = min(vt, 0);
    %     for i=2:Nr-1
    %         upwfet = (vp(i)*(fet(i)-fet(i-1))+vn(i)*(fet(i+1)-fet(i)))/dx; %first order upwind
    %         upwc(i) = (vp(i)*(c(i)-c(i-1))+vn(i)*(c(i+1)-c(i)))/dx;
    %         fetnew(i) = fet(i)-dt.*upwfet+dt.*fet(i).*(v(i)./(r(i).*R(n-1))-gamman(i)/3-1/3*b(i).*(fet(i)^2-fet(i)^-4));
    %     end
    %       for i=2:Nr-1
    %         if i > 2 && i < Nr-1
    %             upwfet = (vp(i)*(3*fet(i)-4*fet(i-1)+fet(i-2))+vn(i)*(-fet(i+2)+4*fet(i+1)-3*fet(i)))/(2*dx); %second order upwind
    %             %upwc(i) = (vp(i)*(3*c(i)-4*c(i-1)+c(i-2))+vn(i)*(-c(i+2)+4*c(i+1)-3*c(i)))/(2*dx);
    %         else
    %             upwfet = (vp(i)*(fet(i)-fet(i-1))+vn(i)*(fet(i+1)-fet(i)))/dx; %first order upwind
    %             %upwc(i) = (vp(i)*(c(i)-c(i-1))+vn(i)*(c(i+1)-c(i)))/dx;
    %         end
    %         %upwfet = (vp(i)*(fet(i+1)-fet(i-1))+vn(i)*(fet(i+1)-fet(i-1)))/(2*dx); %second order upwind
    %         %upwc(i) = (vp(i)*(c(i)-c(i-1))+vn(i)*(c(i+1)-c(i)))/dx;
    %         fetnew(i) = fet(i)-dt.*upwfet+dt.*fet(i).*(v(i)./(r(i).*R(n-1))-gamman(i)/3-1/3*b(i).*(fet(i)^2-fet(i)^-4));
    %       end

    Fet_r_p = zeros(Nr,1);
    Fet_r_p(1:end-2) = (-3.*fet(1:end-2)+4.*fet(2:end-1)-fet(3:end))./2./dx;
    Fet_r_p(end-1) = (fet(end)-fet(end-1))./dx;
    Fet_r_m = zeros(Nr,1);
    Fet_r_m(3:end) = (3.*fet(3:end)-4.*fet(2:end-1)+fet(1:end-2))./2./dx;
    Fet_r_m(2) = (fet(2)-fet(1))./dx;
    D2 = gamman/3+1/3*b.*(fet.^2-fet.^(-4));
    fetnew(2:end) = fet(2:end)-dt*(max(vt(2:end),0).*Fet_r_m(2:end)+min(vt(2:end),0).*Fet_r_p(2:end))+dt*(v(2:end)./r(2:end)./R(n-1)-D2(2:end)).*fet(2:end);
    fetnew(1) = fet(1)+dt*(v_r(1)-D2(1)).*fet(1);


    c_r_p = zeros(Nr,1);
    c_r_p(1:end-2) = (-3.*c(1:end-2)+4.*c(2:end-1)-c(3:end))./2./dx;
    c_r_p(end-1) = (c(end)-c(end-1))./dx;
    c_r_m = zeros(Nr,1);
    c_r_m(3:end) = (3.*c(3:end)-4.*c(2:end-1)+c(1:end-2))./2./dx;
    c_r_m(2) = (c(2)-c(1))./dx;
    upwc(2:end-1) = (max(vt(2:end-1),0).*c_r_m(2:end-1)+min(vt(2:end-1),0).*c_r_p(2:end-1));

    if n == 2
        max_inter = 100;
    else
        max_inter = 1;
    end
    for i = 1:max_inter
        diagC    = [-3/2;(1+dt*(gammac+gamman(2:end-1))+dt*2*Drho/(dx)^2/R(n-1)^2);1];
        %diagC    = [1;dt*2*Drho/(dx)^2/R(n-1)^2);1];%Huaming c
        supDiagC = [2;-dt*Drho/((dx)*R(n-1)^2).*(1/dx+1./r(2:Nr-1))];
        subDiagC = [-dt*Drho/(dx*R(n-1)^2).*(1/dx-1./r(2:Nr-1));0];
        Amatrix_c = sparse(1:Nr,1:Nr,diagC,Nr,Nr)...
            + sparse(1:Nr-1,2:Nr,supDiagC,Nr,Nr)...
            + sparse(2:Nr,1:Nr-1,subDiagC,Nr,Nr)...
            +sparse(1,3,-1/2,Nr,Nr);
        Amatrix_c = full(Amatrix_c);
        Bmatrix_c = [0;c(2:end-1);1]-dt*upwc;
        cnew = Amatrix_c\Bmatrix_c;
        c = cnew;
    end

    %    c = sinh(R(n-1)*r./Lbase)./(r.*sinh(R(n-1)./Lbase));
    %    c(1) = R(n-1)/(sinh(R(n-1)/Lbase)*Lbase); %Huaming c in explicit form
    %c = ones(size(c));
    fernew = fetnew.^(-2);


    % Calculate F and G
    Frnew = zeros(Nr,1);
    Ftnew = zeros(Nr,1);
    Grnew = zeros(Nr,1);
    Gtnew = zeros(Nr,1);
    %     % Boundary condition
    %     %RHS F
    Frnew(1)   = Fr(1)  +dt.*Fr(1).*(v_r(1)+bbar(1)*(1-Fr(1)^2));
    Frnew(end) = Fr(end)+dt.*Fr(end).*(v_r(end)+bbar(end)*(1-Fr(end)^2));
    Ftnew(1)   = Ft(1)  +dt.*Ft(1).*(gamman(1)/3+bbar(1)*(1-Ft(1)^2));
    Ftnew(end) = Ft(end)+dt.*Ft(end).*(v(end)./(r(end).*R(n-1))+bbar(end)*(1-Ft(end)^2));


    Grnew(1)   = Gr(1)+dt.*(gamman(1)/3+2/3*b(1).*(-fet(1)^2+fet(1)^(-4)) + bbar(1)*(1-(Fr(1))^2)).*Gr(1) ;
    Grnew(end) = Gr(end)+dt.*(gamman(end)/3+2/3*b(end).*(-fet(end)^2+fet(end)^(-4)) + bbar(end)*(1-(Fr(end))^2)).*Gr(end);
    Gtnew(1)   = Gt(1)+dt.*(gamman(1)/3+1/3*b(1).*(-fet(1)^(-4)+fet(1)^2) + bbar(1)*(1-(Ft(1))^2)).*Gt(1);
    Gtnew(end) = Gt(end)+dt.*(gamman(end)/3+1/3*b(end).*(-fet(end)^(-4)+fet(end)^2) + bbar(end)*(1-(Ft(end))^2)).*Gt(end);

    for i=2:Nr-1
        upwFr = (vp(i)*(Fr(i)-Fr(i-1))+vn(i)*(Fr(i+1)-Fr(i)))/dx;
        upwFt = (vp(i)*(Ft(i)-Ft(i-1))+vn(i)*(Ft(i+1)-Ft(i)))/dx;
        upwGr = (vp(i)*(Gr(i)-Gr(i-1))+vn(i)*(Gr(i+1)-Gr(i)))/dx;
        upwGt = (vp(i)*(Gt(i)-Gt(i-1))+vn(i)*(Gt(i+1)-Gt(i)))/dx;

        Frnew(i) = Fr(i)-dt.*upwFr+dt.*Fr(i).*(v_r(i)+bbar(i)*(1-Fr(i)^2));
        Ftnew(i) = Ft(i)-dt.*upwFt+dt.*Ft(i).*(v(i)./(r(i).*R(n-1))+bbar(i)*(1-Ft(i)^2));

        Grnew(i) = Gr(i)-dt.*upwGr+dt.*(gamman(i)/3+2/3*b(i).*(fet(i)^(-4)-fet(i)^2) + bbar(i)*(1-(Fr(i))^2)).*Gr(i);
        Gtnew(i) = Gt(i)-dt.*upwGt+dt.*(gamman(i)/3+1/3*b(i).*(fet(i)^2-fet(i)^(-4)) + bbar(i)*(1-(Ft(i))^2)).*Gt(i);

    end

    fet = fetnew;
    fer = fernew;
    Ft =  Ftnew;%fetnew.*Gtnew;%Ftnew;
    Fr =  Frnew;%fernew.*Grnew;%Frnew;
    Gt = Gtnew;%fetnew.^(-1).*Ft;%Gtnew;
    Gr = Grnew;%fernew.^(-1).*Fr;%Grnew;

    %% bookkeeping
    if ~hoop_warned && max(abs(sr))>1e4 % do not display warning when fitting
        warning(['T=' num2str(currTime) ', sigma_rr-p too large = ' num2str(max(abs(sr))),' fet=',num2str(min((fet)))]);
        hoop_warned = true;
        break
        %     save('fetbug.mat','fe_t','gamman','v','r','R')
        %     error('stop')
    end

    if mod(n-1,record_every)==0
        nn = (n-1)/record_every+1;
        fe_r(:,nn) = fer;
        fe_t(:,nn) = fet;
        F_r(:,nn) = Fr;
        F_t(:,nn) = Ft;
        G_r(:,nn) = Gr;
        G_t(:,nn) = Gt;
        P(:,nn) = p;
        V(:,nn) = v;
        sr2(:,nn) = sr;
        st2(:,nn) = st;
        VT(:,nn) = vt;
        C(:,nn) = c;
        VA(:,nn) = va;
        PRA(:,nn) = pra;
        Gam(:,nn) = gamman;
    end
end
radial = sr2 - P;
hoop = st2 - P;
end

function m = mina(a,b)
if abs(a)<abs(b)
    m = a;
else
    m = b;
end
end
