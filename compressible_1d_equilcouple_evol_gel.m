%  clear
% clc;

clear citer gamma_iter Fer_iter Fet_iter str_rr_iter str_tt_iter viter div_str_iter
%% parameters
global dr alpha beta mu K Nt N dt tol iter_max L R0 D cH
N = 1000;
R1 = 0;
R2 = 1;
r = linspace(R1, R2, N + 1)'; 
dr = (R2-R1)/ N;

alpha = 0;
beta = 0.3; 
beta_bar1 = 100;
beta_bar2 = 100;
cH = 0;

rho1 = 1000;
rho2 = 1000;
K = 2;%rho1 + (rho2-rho1)./(R2-R1).*(r-R1);
mu = 1.*ones(size(r));
dmu = 0;
dK = 0;
pext=0;

L = 0.3;
eta = 1;
d = 3;
lambda = 1/3;
k=1;
gamma_c = 0.6;
% myGamma = @(lambda,C,Je,W) eta./k.*C.*((K.*(Je-1)-W./Je).*lambda.*d+k./2./Je.*C.^2.*(d*lambda+2));
myGamma = @(C,Je,W) eta.*(k./2./Je.*C.^2+(K.*(Je-1)-W./Je)).*C;
myW = @(Je,Fer,Fet) 1./2.*mu.*(Je.^(-2/d).*(Fer.^2+2.*Fet.^2)-d)+1./2.*K.*(Je-1).^2;

D = L^2;

%% initial data
R0= 1.3;
R = R0;
Fer = ones(size(r));
Fet = ones(size(r));
Jgr = ones(size(r));
Jgt = ones(size(r));

Jr = Fer.*Jgr;
Jt = Fet.*Jgt;



v = zeros(size(r));
vr = zeros(size(r));
vv = (v-v(end).*r)./R;
div_v = zeros(size(r));
rr = r.*R;
c = R/sinh(R/L).*sinh(rr./L)./rr;
c(1) = R/sinh(R/L)*cosh(rr(1)/L)/L; 
gamma = zeros(size(r));  

% error = zeros(iter_max,1);

dt = 1e-3;
Nt = round(100/dt);
tol = 1e-8;
iter_max = 100000;
iter_number = zeros(Nt,1);


Riter = zeros(Nt+1,1);
str_rr = zeros(N+1,1);
str_tt = zeros(N+1,1);
str_rr_r = zeros(N+1,1);

Riter(1) = R;        
        
%% 
tic;
for i = 1:Nt
    Fer_i = Fer;
    Fet_i = Fet;
    Jgr_i = Jgr;
    Jgt_i = Jgt;
    Jr_i = Jr;
    Jt_i = Jt;
    R_i = R;
    v_i = v;
    vv_i = vv;
    c_i = c;

%     iter = 0;
%     error = 100;    
%     while error> 1e-6
%     iter = iter+1;
%     Fer_copy = Fer;
%     Fet_copy = Fet;
%     R_copy = R;
%     v_copy = v;
%     c_copy = c;
%     gamma_copy = gamma;    
    
        rr = r.*R;
    Je = Fer.*Fet.^2;
    rho = 1./Je;
    W = myW(Je,Fer,Fet);  
    fb_m = K.*(Je-1)-W./Je;

    gamma = myGamma(c,Je,W);
    
    
%     Fext = -cH/2*(5-R0*(R0^3+4*R^3)/R^4);
%     dFext_coef = cH*2*R0*(R0^3/R^5+1/R^2);
    Fext = 0;
    dFext_coef = 0;

    
    
%     Fer_r = zeros(size(r));    
%     Fer_r(1) = (-3*Fer(1)+4*Fer(2)-Fer(3))/2/dr;
%     Fer_r(2:end-1) = (Fer(3:end)-Fer(1:end-2))/2/dr;
%     Fer_r(end) = (Fer(end-2)-4*Fer(end-1)+3*Fer(end))/2/dr;
%     Fet_r = zeros(size(r));    
%     Fet_r(1) = (-3*Fet(1)+4*Fet(2)-Fet(3))/2/dr;
%     Fet_r(2:end-1) = (Fet(3:end)-Fet(1:end-2))/2/dr;
%     Fet_r(end) = (Fet(end-2)-4*Fet(end-1)+3*Fet(end))/2/dr;

    Fer_r = zeros(size(r)); 
    Fer_r(1) = Fer(1:5)'*lagrange_p_der(r(1),r(1:5)',Fer(1:5)');
    Fer_r(2) = Fer(1:5)'*lagrange_p_der(r(2),r(1:5)',Fer(1:5)');
    Fer_r(3:end-2) = (Fer(1:end-4)-Fer(5:end))/12/dr+(Fer(4:end-1)-Fer(2:end-3))*2/3/dr;
    Fer_r(end-1) = Fer(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',Fer(end-4:end)');
    Fer_r(end) = Fer(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',Fer(end-4:end)');
     
    Fet_r = zeros(size(r)); 
    Fet_r(1) = Fet(1:5)'*lagrange_p_der(r(1),r(1:5)',Fet(1:5)');
    Fet_r(2) = Fet(1:5)'*lagrange_p_der(r(2),r(1:5)',Fet(1:5)');
    Fet_r(3:end-2) = (Fet(1:end-4)-Fet(5:end))/12/dr+(Fet(4:end-1)-Fet(2:end-3))*2/3/dr;
    Fet_r(end-1) = Fet(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',Fet(end-4:end)');
    Fet_r(end) = Fet(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',Fet(end-4:end)');
    
    AA = mu.*(Fer.*Fet.^2).^(-5/3).*(Fer.^2-Fet.^2);
    BB = K.*(Fer.*Fet.^2-1);
    

    AA_er = (-5/3).*mu.*(Fer.*Fet.^2).^(-8/3).*Fet.^2.*(Fer.^2-Fet.^2)+mu.*(Fer.*Fet.^2).^(-5/3).*2.*Fer;
    AA_erer = -2./9.*mu.*(Fer.*Fet.^2).^(-5/3)-40./9.*mu.*(Fer.*Fet.^2).^(-11/3).*Fet.^6;
    AA_eret = -10./9.*mu.*(Fer.*Fet.^2).^(-8/3).*Fet.*(Fer.^2+2.*Fet.^2);
    
    AA_et = (-5/3).*mu.*(Fer.*Fet.^2).^(-8/3).*2.*Fer.*Fet.*(Fer.^2-Fet.^2)...
            + mu.*(Fer.*Fet.^2).^(-5/3).*(-2.*Fet);
    AA_eter = -10./9.*mu.*(Fer.*Fet.^2).^(-8/3).*Fet.*(Fer.^2+2.*Fet.^2);
    AA_etet = 130./9.*mu.*(Fer.*Fet.^2).^(-8/3).*Fer.^3-28./9.*(Fer.*Fet.^2).^(-5/3);
    
    BB_er = K.*Fet.^2;
    BB_et = 2.*K.*Fer.*Fet;
    BB_erer = 0; BB_eret = 2.*K.*Fet; BB_eter = 2.*K.*Fet; BB_etet = 2.*K.*Fer;

    A = 2/3.*AA_er+BB_er;
    B = 2/3.*AA_et+BB_et;
    
    C1 = (2/3.*AA_erer+BB_erer).*Fer_r./R+(2/3.*AA_eret+BB_eret).*Fet_r./R+2./(r.*R).*AA_er;
    C2 = (2/3.*AA_eter+BB_eter).*Fer_r./R+(2/3.*AA_etet+BB_etet).*Fet_r./R+2./(r.*R).*AA_et;
    
    D1 = gamma/3+2./3.*beta.*(Fer.^2-Fet.^2);
    D2 = gamma/3+1./3.*beta.*(Fet.^2-Fer.^2);
    
    
    gamma_r = zeros(size(r));
    gamma_r(1) = gamma(1:5)'*lagrange_p_der(r(1),r(1:5)',gamma(1:5)');
    gamma_r(2) = gamma(1:5)'*lagrange_p_der(r(2),r(1:5)',gamma(1:5)');
    gamma_r(3:end-2) = (gamma(1:end-4)-gamma(5:end))/12/dr+(gamma(4:end-1)-gamma(2:end-3))*2/3/dr;
    gamma_r(end-1) = gamma(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',gamma(end-4:end)');
    gamma_r(end) = gamma(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',gamma(end-4:end)');
    
    

    D1_r = gamma_r./3+2./3.*beta.*(2.*Fer.*Fer_r-2.*Fet.*Fet_r);
    D2_r = gamma_r./3+1./3.*beta.*(2.*Fet.*Fet_r-2.*Fer.*Fer_r);
    
%     D1_r = zeros(size(r));    
%     D1_r(1) = D1(1:5)'*lagrange_p_der(r(1),r(1:5)',D1(1:5)');
%     D1_r(2) = D1(1:5)'*lagrange_p_der(r(2),r(1:5)',D1(1:5)');
%     D1_r(3:end-2) = (D1(1:end-4)-D1(5:end))/12/dr+(D1(4:end-1)-D1(2:end-3))*2/3/dr;
%     D1_r(end-1) = D1(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',D1(end-4:end)');
%     D1_r(end) = D1(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',D1(end-4:end)');
%     
%     D2_r = zeros(size(r));
%     D2_r(1) = D2(1:5)'*lagrange_p_der(r(1),r(1:5)',D2(1:5)');
%     D2_r(2) = D2(1:5)'*lagrange_p_der(r(2),r(1:5)',D2(1:5)');
%     D2_r(3:end-2) = (D2(1:end-4)-D2(5:end))/12/dr+(D2(4:end-1)-D2(2:end-3))*2/3/dr;
%     D2_r(end-1) = D2(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',D2(end-4:end)');
%     D2_r(end) = D2(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',D2(end-4:end)');
    
    
    
    a1 = C2.*Fet./(r.*R)-B.*Fet./(r.*R).^2+B.*Fet_r./R./r./R-2.*AA./(r.*R).^2;
    a1 = a1(2:end-1);
    a2 = C1.*Fer+B.*Fet./(r.*R)-B.*Fet_r/R;
    a2 = a2(2:end-1);
    a3 = A.*Fer;
    a3 = a3(2:end-1);
    a4 = -C1.*Fer.*D1-C2.*Fet.*D2;
    a4 = a4(2:end-1)-A(2:end-1).*(D1_r(2:end-1).*Fer(2:end-1)+D1(2:end-1).*Fer_r(2:end-1))./R...
        -B(2:end-1).*(D2_r(2:end-1).*Fet(2:end-1)+D2(2:end-1).*Fet_r(2:end-1))./R;
    

    str_rr = 2/3*AA+BB;
    str_tt = -1/3*AA+BB;    
    
%     str_rr_r(1) = (-3*str_rr(1)+4*str_rr(2)-str_rr(3))/2/dr;
%     str_rr_r(2:end-1)=(str_rr(3:end)-str_rr(1:end-2))/2/dr;
%     str_rr_r(end) = (str_rr(end-2)-4*str_rr(end-1)+3*str_rr(end))/2/dr;
    
    str_rr_r(1) = str_rr(1:5)'*lagrange_p_der(r(1),r(1:5)',str_rr(1:5)');
    str_rr_r(2) = str_rr(1:5)'*lagrange_p_der(r(2),r(1:5)',str_rr(1:5)');
    str_rr_r(3:end-2) = (str_rr(1:end-4)-str_rr(5:end))/12/dr+(str_rr(4:end-1)-str_rr(2:end-3))*2/3/dr;
    str_rr_r(end-1) = str_rr(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',str_rr(end-4:end)');
    str_rr_r(end) = str_rr(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',str_rr(end-4:end)');
  
    div_str = A.*Fer_r./R+B.*Fet_r./R+2./(r.*R).*AA;
    
    a4 = a4+alpha/dt*v_i(2:end-1)+beta_bar1.*div_str(2:end-1);

    v_r_p = zeros(N+1,1);
    v_r_p(1:end-2) = (-3.*v(1:end-2)+4.*v(2:end-1)-v(3:end))./2./dr;
    v_r_p(end-1) = (v(end)-v(end-1))./dr;
    v_r_m = zeros(N+1,1);
    v_r_m(3:end) = (3.*v(3:end)-4.*v(2:end-1)+v(1:end-2))./2./dr;
    v_r_m(2) = (v(2)-v(1))./dr;
    

    %%% check plus or minus sign
%     a4 = a4 + alpha*(max(vv_i(2:end-1),0).*v_r_m(2:end-1)+min(vv_i(2:end-1),0).*v_r_p(2:end-1));
    a4 = a4 - alpha*(max(vv_i(2:end-1),0).*v_r_m(2:end-1)+min(vv_i(2:end-1),0).*v_r_p(2:end-1));
   
    
    
    
%     myA = sparse(1, 1, 1,N+1,N+1)-sparse(2:N,2:N,a1,N+1,N+1)...
%             -sparse(2,1:5,[-3,-10,18,-6,1]/12/dr/R.*a2(1)+[11,-20,6,4,-1]/12/dr^2/R^2.*a3(1),N+1,N+1)...            
%             -sparse(3:N-1,1:N-3,1/12/dr/R.*a2(2:end-1)-1/12/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
%             -sparse(3:N-1,2:N-2,-2/3/dr/R.*a2(2:end-1)+4/3/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
%             -sparse(3:N-1,3:N-1,-5/2/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
%             -sparse(3:N-1,4:N,2/3/dr/R*a2(2:end-1)+4/3/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
%             -sparse(3:N-1,5:N+1,-1/12/dr/R.*a2(2:end-1)-1/12/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
%             -sparse(N,N-3:N+1,[-1,6,-18,10,3]./12./dr./R.*a2(end)+[-1,4,6,-20,11]./12./(dr*R)^2.*a3(end),N+1,N+1)...
%             +sparse(N+1,N+1,B(end)*Fet(end)/r(end)/R-dFext_coef,N+1,N+1)...
%             +sparse(N+1,N-3:N+1,A(end)*Fer(end).*[3,-16,36,-48,25]./12./dr./R,N+1,N+1);         
%     myb = [0;a4;A(end)*Fer(end)*D1(end)+B(end)*Fet(end)*D2(end)-beta_bar2.*(str_rr(end)-Fext)];
    
    myA = sparse(1, 1, 1,N+1,N+1)+sparse(2:N,2:N,alpha*(1/dt+beta_bar1),N+1,N+1)...
            -sparse(2:N,2:N,a1,N+1,N+1)...
            -sparse(2,1:5,[-3,-10,18,-6,1]/12/dr/R.*a2(1)+[11,-20,6,4,-1]/12/dr^2/R^2.*a3(1),N+1,N+1)...            
            -sparse(3:N-1,1:N-3,1/12/dr/R.*a2(2:end-1)-1/12/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
            -sparse(3:N-1,2:N-2,-2/3/dr/R.*a2(2:end-1)+4/3/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
            -sparse(3:N-1,3:N-1,-5/2/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
            -sparse(3:N-1,4:N,2/3/dr/R*a2(2:end-1)+4/3/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
            -sparse(3:N-1,5:N+1,-1/12/dr/R.*a2(2:end-1)-1/12/(dr*R)^2.*a3(2:end-1),N+1,N+1)...
            -sparse(N,N-3:N+1,[-1,6,-18,10,3]./12./dr./R.*a2(end)+[-1,4,6,-20,11]./12./(dr*R)^2.*a3(end),N+1,N+1)...
            +sparse(N+1,N+1,B(end)*Fet(end)/r(end)/R-dFext_coef,N+1,N+1)...
            +sparse(N+1,N-3:N+1,A(end)*Fer(end).*[3,-16,36,-48,25]./12./dr./R,N+1,N+1);         
    myb = [0;a4;A(end)*Fer(end)*D1(end)+B(end)*Fet(end)*D2(end)-beta_bar2.*(str_rr(end)-Fext)];   


    
    
    v = myA\myb;
    v(1) = 0;
    
%     R = R_i+dt*v(end); 
    v_r = zeros(N+1,1);
%     v_r(1) = (-3*v(1)+4*v(2)-v(3))/2/dr;
%     v_r(2:end-1)=(v(3:end)-v(1:end-2))/2/dr;
%     v_r(end) = (v(end-2)-4*v(end-1)+3*v(end))/2/dr;
    
    v_r(1) = v(1:5)'*lagrange_p_der(r(1),r(1:5)',v(1:5)');
    v_r(2) = v(1:5)'*lagrange_p_der(r(2),r(1:5)',v(1:5)');
    v_r(3:end-2) = (v(1:end-4)-v(5:end))/12/dr+(v(4:end-1)-v(2:end-3))*2/3/dr;
    v_r(end-1) = v(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',v(end-4:end)');
    v_r(end) = v(end-4:end)'*lagrange_p_der(r(end),r(end-4:end)',v(end-4:end)');
    
    vv = (v-v(end).*r)./R;
    
  
%     Fer_r_p = zeros(N+1,1);
%     Fer_r_p(2:end-3) = (-3.*Fer(1:end-4)-10.*Fer(2:end-3)+18.*Fer(3:end-2)-6.*Fer(4:end-1)+Fer(5:end))./12./dr;
%     Fer_r_p(end-2) = Fer(end-4:end)'*lagrange_p_der(r(end-2),r(end-4:end)',Fer(end-4:end)');
%     Fer_r_p(end-1) = Fer(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',Fer(end-4:end)');
%     Fer_r_m = zeros(N+1,1);
%     Fer_r_m(4:end-1) = (-Fer(1:end-4)+6.*Fer(2:end-3)-18.*Fer(3:end-2)+10.*Fer(4:end-1)+3.*Fer(5:end))./12./dr;
%     Fer_r_m(3) = Fer(1:5)'*lagrange_p_der(r(3),r(1:5)',Fer(1:5)');    
%     Fer_r_m(2) = Fer(1:5)'*lagrange_p_der(r(2),r(1:5)',Fer(1:5)');
% 
%     Fet_r_p = zeros(N+1,1);
%     Fet_r_p(2:end-3) = (-3.*Fet(1:end-4)-10.*Fet(2:end-3)+18.*Fet(3:end-2)-6.*Fet(4:end-1)+Fet(5:end))./12./dr;
%     Fet_r_p(end-2) = Fet(end-4:end)'*lagrange_p_der(r(end-2),r(end-4:end)',Fet(end-4:end)');
%     Fet_r_p(end-1) = Fet(end-4:end)'*lagrange_p_der(r(end-1),r(end-4:end)',Fet(end-4:end)');
%     Fet_r_m = zeros(N+1,1);
%     Fet_r_m(4:end-1) = (-Fet(1:end-4)+6.*Fet(2:end-3)-18.*Fet(3:end-2)+10.*Fet(4:end-1)+3.*Fet(5:end))./12./dr;
%     Fet_r_m(3) = Fet(1:5)'*lagrange_p_der(r(3),r(1:5)',Fet(1:5)');    
%     Fet_r_m(2) = Fet(1:5)'*lagrange_p_der(r(2),r(1:5)',Fet(1:5)');

    Fer_r_p = zeros(N+1,1);
    Fer_r_p(1:end-2) = (-3.*Fer(1:end-2)+4.*Fer(2:end-1)-Fer(3:end))./2./dr;
    Fer_r_p(end-1) = (Fer(end)-Fer(end-1))./dr;
    Fer_r_m = zeros(N+1,1);
    Fer_r_m(3:end) = (3.*Fer(3:end)-4.*Fer(2:end-1)+Fer(1:end-2))./2./dr;
    Fer_r_m(2) = (Fer(2)-Fer(1))./dr;
    
    Fet_r_p = zeros(N+1,1);
    Fet_r_p(1:end-2) = (-3.*Fet(1:end-2)+4.*Fet(2:end-1)-Fet(3:end))./2./dr;
    Fet_r_p(end-1) = (Fet(end)-Fet(end-1))./dr;
    Fet_r_m = zeros(N+1,1);
    Fet_r_m(3:end) = (3.*Fet(3:end)-4.*Fet(2:end-1)+Fet(1:end-2))./2./dr;
    Fet_r_m(2) = (Fet(2)-Fet(1))./dr;
    
     
    Fer = Fer_i-dt*(max(vv,0).*Fer_r_m+min(vv,0).*Fer_r_p)+dt*(v_r./R-D1).*Fer;
    Fet(2:end) = Fet_i(2:end)-dt*(max(vv(2:end),0).*Fet_r_m(2:end)+min(vv(2:end),0).*Fet_r_p(2:end))+dt*(v(2:end)./r(2:end)./R-D2(2:end)).*Fet(2:end);
    Fet(1) = Fet_i(1)+dt*(v_r(1)./R-D2(1)).*Fet(1);
    
    c_r_p = zeros(N+1,1);
    c_r_p(1:end-2) = (-3.*c(1:end-2)+4.*c(2:end-1)-c(3:end))./2./dr;
    c_r_p(end-1) = (c(end)-c(end-1))./dr;
    c_r_m = zeros(N+1,1);
    c_r_m(3:end) = (3.*c(3:end)-4.*c(2:end-1)+c(1:end-2))./2./dr;
    c_r_m(2) = (c(2)-c(1))./dr;
    
    
    r2rho = r.^2.*rho;
    Aa = sparse(1,1:3,[1,-4/3,1/3],N+1,N+1)+sparse(N+1,N+1,1,N+1,N+1)+...
        sparse(2:N,2:N,1+dt.*(gamma_c+gamma(2:end-1))+dt.*D./rho(2:end-1)./r(2:end-1).^2./2./(dr*R)^2.*(r2rho(1:end-2)+2.*r2rho(2:end-1)+r2rho(3:end)),N+1,N+1)+...
        sparse(2:N,1:N-1,-dt.*D./rho(2:end-1)./r(2:end-1).^2./2./(dr*R)^2.*(r2rho(1:end-2)+r2rho(2:end-1)),N+1,N+1)+...
        sparse(2:N,3:N+1,-dt.*D./rho(2:end-1)./r(2:end-1).^2./2./(dr*R)^2.*(r2rho(3:end)+r2rho(2:end-1)),N+1,N+1);

    bb = zeros(N+1,1);
    bb(2:end-1) = c_i(2:end-1)-dt*(max(vv(2:end-1),0).*c_r_m(2:end-1)+min(vv(2:end-1),0).*c_r_p(2:end-1));
    bb(N+1)=1;
    c = Aa\bb;

    R = R_i+dt*v(end); 
    Riter(i+1) = R;


    if i==1
        citer = c;
        gamma_iter=gamma;
        Fer_iter = Fer;
        Fet_iter = Fet;
        str_rr_iter = str_rr;
        str_tt_iter = str_tt;
        viter = v;
        div_str_iter = div_str;
    elseif rem(i,1/dt)==0 
%     else
        citer(:,end+1) = c;
        gamma_iter(:,end+1)=gamma;
        Fer_iter(:,end+1) = Fer;
        Fet_iter(:,end+1) = Fet;
        str_rr_iter(:,end+1) = str_rr;
        str_tt_iter(:,end+1) = str_tt;
        viter(:,end+1) = v;
        div_str_iter(:,end+1) = div_str;
    end

end
time = toc;
%%
% filename = sprintf('Sphere1d_relax_nofriction_data unifrom-growth gamma=%1.e K=%d beta=%1.e R0=%1.1e Nt=%d t=%d N=%d.mat',gamma(1),K(1),beta,R0,Nt,dt*Nt,N);
% filename = sprintf('Sphere1d_relax_rd_nofb_damp_data K=%d beta=%1.e beta_bar1=%1.e beta_bar2=%1.e R0=%1.1e Nt=%d t=%d N=%d tol=%1.e.mat',K(1),beta,beta_bar1,beta_bar2,R0,Nt,dt*Nt,N,tol);
filename = sprintf('Sphere1d_gel_crd_oldfb_3rd_damp_data_conserv_evol alpha=%1.2e cH=%1.2e L=%1.e K=%d beta=%1.2e R0=%1.2e Nt=%d t=%d N=%d tol=%1.e.mat',alpha,cH,L,K(1),beta,R0,Nt,dt*Nt,N,tol);
% save(filename,'Fer_iter','Fet_iter','Riter','Eiter','viter','citer','str_tt_iter','str_rr_iter','beta','K','r','dr','dt','Nt','iter_number','time','-v7.3');
save(filename,'-v7.3');
