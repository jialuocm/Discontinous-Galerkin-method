% This code solves a pure advective problems ; 
clc
clear 
close all
format long

% Define the physical problem
%% 
%advective coefficent 
u = 1.0;
% 
f = @(x) cos(2*pi*x);
g0 = 1;
% physical domain
omega_l = 0.0; omega_r = 1.0; 
% exact result 
exact = @(x) (1./(2*pi*u) * sin(2*pi*x)+g0);
exact_x = @(x) (1./(u) * cos(2*pi*x));
%-------------------------------------------------------------------
nElem_x = 3: 3: 36;
pp_x = 1:1:3;
cn = length(nElem_x);
p_max = length(pp_x);

hh_x = zeros(cn,1);
error_l2_x = zeros(cn,p_max);
error_l2s_x = zeros(cn,p_max);
error_h1_x = zeros(cn,p_max);

for pn = 1: p_max
    for num = 1 : cn 
% number of elements
nElem = nElem_x(num);
% choose polynomial (FEM basis function) degree
pp = pp_x(pn);
%-------------------------------------------------------------------
%DG_1D_advection equation solver
%-------------------------------------------------------------------
n_en = pp + 1; % number of nodal points of element
%n_an = pp * nElem + 1; % number of nodal points of all domain

%mesh uniform size
LL = omega_r - omega_l;
hh = LL./nElem;
tau = hh / (2* abs(u));

%discretization phiscal domain
x_coor = omega_l: hh/pp: omega_r;
% x_en = omega_l: hh: omega_r;

% location Matrix
IEN = zeros(pp+1, nElem);
for ee = 1 : nElem
    for aa = 1 : pp+1
        IEN(aa, ee) = (ee - 1) * pp + aa;
    end
end

%preallocate memory for the problem
Fai_h = zeros(nElem,n_en);


% quadrature rule
nqp = pp + 1;
[qp, wq] = Gauss( nqp, -1, 1 );
figure
for ee = 1: nElem
    k_ele = zeros(n_en,n_en);
    f_ele = zeros(n_en,1);
    x_ele = zeros(n_en,1);
    faih = zeros(n_en,1);

    for aa= 1: n_en
        x_ele(aa) = x_coor( IEN(aa, ee));
    end

    for qua = 1: nqp
        %Geometrical mapping 
        x_qua = 0;
        dx_dxi = 0;
        for aa= 1 : n_en
            x_qua = x_qua + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1: n_en
            Na    = PolyBasis(pp, aa, 0, qp(qua));
            Na_xi = PolyBasis(pp, aa, 1, qp(qua));
            f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi + ...
                        wq(qua) * tau * u * Na_xi * f(x_qua);

            for bb = 1: n_en
                Nb    = PolyBasis(pp, bb, 0, qp(qua));
                Nb_xi = PolyBasis(pp, bb, 1, qp(qua));
                k_ele(aa, bb) =  k_ele(aa, bb) + (-u * wq(qua)*Na_xi*Nb) + ...
                                 wq(qua) * tau * u^2 * Na_xi * Nb_xi * dxi_dx;
            end

        end

    end

    k_ele(n_en, n_en) =  k_ele(n_en, n_en) + u*PolyBasis(pp,n_en,0,1);

    if ee==1
        f_ele(1) = f_ele(1) + u*PolyBasis(pp,1,0,-1)*g0;
    else
        f_ele(1) = f_ele(1) + u*PolyBasis(pp,1,0,-1)*Fai_h(ee-1,n_en);
    end

    faih = k_ele\f_ele;

    for ii = 1 :n_en
        Fai_h(ee,ii) = faih(ii);
    end
    pl = plot(x_ele, faih,'color', 'r','linewidth',1);
    hold on;
 end
xx = linspace(omega_l,omega_r,100);
pl_exact = plot(xx, exact(xx),'b--','LineWidth',2);
legend([pl,pl_exact],'DG','Exact');
hold off;
exportgraphics(gca,['file_DG'  num2str(pn) num2str(ee) '.jpg']);

% error analysis
error_l2  = 0.0;
error_h1  = 0.0;
error_l2s = 0.0;
% new quadrature rule
nqp = 10;
[qp, wq] = Gauss( nqp, -1, 1 );

%Exact = zeros(ee,1);

for ee = 1: nElem
    x_ele = zeros( n_en,1);
    faih_ele = zeros(n_en,1);
    for aa = 1 : n_en
        x_ele(aa) = x_coor( IEN(aa,ee));
        faih_ele(aa) = Fai_h(ee,aa);
    end

    for qua = 1 : nqp
    %Geometrical mapping 
    dx_dxi = 0;
    x_qua = 0;
    fai_h = 0;
    fai_h_dxi = 0;
        for aa = 1 : n_en
            x_qua     = x_qua + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            dx_dxi    = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            fai_h     = fai_h + faih_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
            fai_h_dxi = fai_h_dxi + faih_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
        end
        dxi_dx = 1.0 / dx_dxi;
        fai = exact(x_qua);
        fai_x = exact_x(x_qua);

%       error_l2 = error_l2 + wq(qua) * (fai_h - fai)^2 * dx_dxi;
        error_h1 = error_h1 + wq(qua) * (fai_h_dxi * dxi_dx- fai_x)^2 * dx_dxi;
        error_l2s = error_l2s + tau * u^2 * wq(qua) * (fai_h_dxi* dxi_dx - fai_x)^2 * dx_dxi;
    end
end

error_l2s =  error_l2s + 0.5 * u *(Fai_h(1,1)-exact(x_coor(1))) * (Fai_h(1,1)-exact(x_coor(1))) + ...
                         0.5 * u *(Fai_h(nElem,n_en)-exact(x_coor(end))) * (Fai_h(nElem,n_en)-exact(x_coor(end)));

for ee = 2 : nElem
    error_l2s = error_l2s + 05 * u * (Fai_h(nElem,1)- Fai_h(nElem-1,n_en))^2;
end

       %error_l2_x(num,pn)  = log(sqrt(error_l2));
        error_h1_x(num,pn)  = log10(sqrt(error_h1));
        error_l2s_x(num,pn) = log10(sqrt(error_l2s));
        %log converge rate as h decrease 
        hh_x(num)= hh;
        hh_lg(num) = log10(hh);
        dx_h = zeros( 2,1);
        dy_e = zeros( 2,1);
    if (num >= 2)
        dx_h = [hh_lg(num);hh_lg(num-1)];
        dy_h1e = [error_h1_x(num,pn);error_h1_x(num-1,pn)];
        dy_l2e = [error_l2s_x(num,pn);error_l2s_x(num-1,pn)];
        slop_h1(num,pn) = (dy_h1e(2) - dy_h1e(1)) / (dx_h(2)-dx_h(1));
        slop_l2s(num,pn) = (dy_l2e(2) - dy_l2e(1)) / (dx_h(2)-dx_h(1));    
    else 
        slop_h1(num,pn)  = 0;
        slop_l2s(num,pn) = 0;
    end

    end

  
    figure
    hold on;
    p_h1 = plot(hh_lg,error_h1_x(:,pn),'--r*');
    legend(p_h1,'H1 convergence estimate');
    xlabel('log(hh)');
    ylabel('log(ERROR_H1)');
    hold off;
    exportgraphics(gca,['file_H1' num2str(pn) '.jpg']);

    figure
    hold on
    p_l2s = plot(hh_lg,error_l2s_x(:,pn),'--bO');
    xlabel('log(hh)');
    ylabel('log(ERROR_l2s)');
    legend(p_l2s,'L2 stab convergence estimate');
    hold off
    exportgraphics(gca,['fileL2_s' num2str(pn) '.jpg']);

end
% error convergence data : table 
T = table(hh_x,error_h1_x,error_l2s_x,slop_h1,slop_l2s,'variableNames',{'hh_mesh','error_H1','error_l2s','H1 convergence rate','l2_s convergence rate'})
writetable(T)


    

       



 


       

 



       






       

















































