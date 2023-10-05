% ========================================================================
% Monte-Carlo Size experiments using DGP1 from 
% "A Novel Approach to Predictive Accuracy Testing in Nested Environments" 
% Econometric Theory, 2023, Jean-Yves Pitarakis. 
% 
% ========================================================================
% 
% Remarks: (i) this code is not meant to replicate exactly the tables in the
% paper, although outcomes will have negligible differences overall. This
% code uses a faster evaluation of recursive forecast errors, allowing for
% a much greater number of replications to be feasible. (ii) The program
% outputs empirical sizes for the S(lam1,lam2) and Sbar(tau0;lam20) family
% of test statistics and places them in excel files. It produces 3 excel
% files corresponding to phi=0.75, phi = 0.95 and phi = 0.98 for the AR
% process used as a predictor. Each excel file contains two sheets for the
% two families of test statistics. 
% 
% DGP: y(t) = beta0 + u(t) (model 1)
% Predictor: x(t) = phi x(t-1) + v(t) 
% i.e. fitted model 2: y(t)=b0+b1 x(t-1)+u(t) (model 2)
% 
% V[u(t)] = 3, V[v(t)] = 0.01, Corr[u(t),v(t)] = -0.8
% phi \in {0.75,0.95,0.98}
% T \in {250, 500, 1000}
% R = 50000 replications
% 
% Test Statistic Parameterisations: 
% S(lam10,lam20) uses lam10 = 1 and lam20 ranged between 0.50 and 0.95
% Sbar(tau0;lam20) uses tau0 = 0.8 and lam20 ranged between 0.50 and 1.00
% 
% ======================================================================

R = 50000; 
pi0 = 0.25;
nvec = [250;500;1000];
[nv,~]=size(nvec);

% Parameterizations for the S(lam10,lam20) statistic

lam10 = 1;
lam20_seq = [0.5:0.05:0.95]';
[ns,~]=size(lam20_seq);

% Parameterizations for the Sbar(tau0;lam20) statistic 

tauvec = [0.80];
[nt,~]=size(tauvec);
lam20_seq_bar = [0.5:0.05:1]';
[nsb,~]=size(lam20_seq_bar);

% Initializations (storage preps) 

S0 = nan(R,ns,nv);
S0_nw = nan(R,ns,nv);
S0_adj = nan(R,ns,nv);
S0_adj_nw = nan(R,ns,nv);
pv_S0 = nan(R,ns,nv);
pv_S0_nw = nan(R,ns,nv);
pv_S0_adj = nan(R,ns,nv);
pv_S0_adj_nw = nan(R,ns,nv);

Sbar = nan(R,nt,nsb,nv);
Sbar_nw = nan(R,nt,nsb,nv);
Sbar_adj = nan(R,nt,nsb,nv);
Sbar_adj_nw = nan(R,nt,nsb,nv);
pv_Sbar = nan(R,nt,nsb,nv);
pv_Sbar_nw = nan(R,nt,nsb,nv);
pv_Sbar_adj = nan(R,nt,nsb,nv);
pv_Sbar_adj_nw = nan(R,nt,nsb,nv);

Size_Sbar  = nan(nt,nsb,nv);
Size_Sbar_nw  = nan(nt,nsb,nv);
Size_Sbar_adj  = nan(nt,nsb,nv);
Size_Sbar_adj_nw  = nan(nt,nsb,nv);

% ========================================

% DGP parameterizations (see dgp.m file) 

p=1;
k=1;
trend = 1:p;

sigsq_u = 3;
Omega_vv = 0.01*eye(p);

correl_uv = -0.8.^trend;
beta0 = 0;
beta1 = zeros(p,1);

phivec = [0.75;0.95;0.98];
[nf,~]=size(phivec);

for d=1:nf

    Phimat = phivec(d)*eye(p);


for g = 1:nv
n = nvec(g);

for J = 1:R
    [y,X] = dgp(n,p,k,sigsq_u,correl_uv,Omega_vv,Phimat,beta0,beta1);
    
ehat1 = recursive_hstep_interceptonly(y,pi0,1);
ehat2 = recursive_hstep_fast(y,X(:,1),pi0,1);

for s = 1:ns
[S0(J,s,g), S0_nw(J,s,g), S0_adj(J,s,g), S0_adj_nw(J,s,g), pv_S0(J,s,g), pv_S0_nw(J,s,g), pv_S0_adj(J,s,g), pv_S0_adj_nw(J,s,g)] = Nested_Stats_S0(ehat1,ehat2,lam10,lam20_seq(s));
end

for t = 1:nt
    for m = 1:nsb
[Sbar(J,t,m,g), Sbar_nw(J,t,m,g), Sbar_adj(J,t,m,g), Sbar_adj_nw(J,t,m,g), pv_Sbar(J,t,m,g), pv_Sbar_nw(J,t,m,g), pv_Sbar_adj(J,t,m,g), pv_Sbar_adj_nw(J,t,m,g)] = Nested_Stats_Sbar(ehat1,ehat2,lam20_seq_bar(m),tauvec(t));
    end
end

end

end

size_S0 = nan(nv,ns);
size_S0_adj = nan(nv,ns);

for g = 1:nv
size_S0(g,:) = sum(pv_S0(:,:,g)<0.10)/R;
size_S0_adj(g,:) = sum(pv_S0_adj(:,:,g)<0.10)/R;
end

str_0 = ["lambda2"; "T=250 S0"; "T=500 S0"; "T=1000 S0"; "T=250 S0adj"; "T=500 S0adj"; "T=1000 S0adj"];
out_S0_S0adj = [str_0, [lam20_seq';size_S0;size_S0_adj]];

for t = 1:nt
   for m = 1:nsb
Size_Sbar(t,m,:) = sum(pv_Sbar(:,t,m,:)<0.10)/R;
Size_Sbar_adj(t,m,:) = sum(pv_Sbar_adj(:,t,m,:)<0.10)/R;
    end
end

str_bar = ["lambda2 (tau0=0.8)"; "T=250 S_bar"; "T=500 S_bar"; "T=1000 S_bar";"T=250 S_bar_adj"; "T=500 S_bar_adj"; "T=1000 S_bar_adj"];
out_Sbar_Sbaradj = [str_bar,[lam20_seq_bar';Size_Sbar(1,:,1);Size_Sbar(1,:,2);Size_Sbar(1,:,3);Size_Sbar_adj(1,:,1);Size_Sbar_adj(1,:,2);Size_Sbar_adj(1,:,3)]];
    

BaseName='Fsize_dgp1_';
filename = [BaseName,num2str(d)];
xlswrite(filename, out_S0_S0adj,'Slam1lam2','A1') 
xlswrite(filename, out_Sbar_Sbaradj,'Sbar(tau)','A1') 

clear out_S0_S0adj;
clear out_Sbar_Sbaradj;

end


