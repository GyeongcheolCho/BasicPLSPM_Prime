function [INI,TABLE,ETC] = BasicPLSPM(z0, W0, B0, modetype,scheme,ind_sign,N_Boot,Max_iter,Min_limit,Flag_Parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BasicGSCA() - MATLAB function to perform a basic version of Partial     %
%               Least Sqaures Path Modeling  (PLSPM).                     %
% Author: Heungsun Hwang & Gyeongcheol Cho                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Data = an N by J matrix of scores for N individuals on J indicators   %
%   W = a J by P matrix of weight parameters                              %
%   B = a P by P matrix of path coefficients                              %
%   modetype = a vector of length P indicating the mode of each latent    %
%              variable (1 = mode A, 2 = mode B)                          %
%   scheme = an integer indicating the scheme for updating inner weights  % 
%              (1 = centroid, 2 = factorial, 3 = path weighting)          %
%   ind_sign = a vector of length P representing dominant indicator of    %
%              each latent variable                                       %
%   N_Boot = Integer representing the number of bootstrap samples for     %
%            calculating standard errors (SE) and 95% confidence          %
%            intervals (CI)                                               %
%   Max_iter = Maximum number of iterations for the Alternating Least     % 
%              Squares (ALS) algorithm                                    %
%   Min_limit = Tolerance level for ALS algorithm                         %
%   Flag_Parallel = Logical value to determine whether to use parallel    %
%                   computing for bootstrapping                           %
% Output arguments:                                                       %
%   INI: Structure array containing goodness-of-fit values, R-squared     % 
%        values, and matrices parameter estimates                         %
%     .iter = Number of iterations for the ALS algorithm                  %
%     .W: a J by P matrix of weight estimates                             %
%     .C: a P by J matrix of loading estimates                            %
%     .B: a P by P matrix of path coefficient estimates                   %
%     .CVscore: an N by P matrix of component scores                      % 
%  TABLE: Structure array containing tables of parameter estimates, their %
%         SEs, 95% CIs,and other statistics                               %
%     .W: Table for weight estimates                                      %
%     .C: Table for loading estimates                                     %
%     .B: Table for path coefficients estimates                           %
%  ETC: Structure array including bootstrapped parameter estmates         %
%     .W_Boot: Matrix of bootstrapped weight estimates                    %
%     .C_Boot: Matrix of bootstrapped loading estimates                   %
%     .B_Boot: Matrix of bootstrapped path coefficient estimates          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:                                                                   %
% 1. This function employs the cost function unifying Mode A and Mode B   %
%    of the PLS algorithm as follows:                                     %
%    "f = alpha_j*SS(Z - FW') + beta_j*SS(F - ZW), where F = Gamma*B"     %
%      alpha_j + beta_j = 1                                               %
%      alpha_j = 1: Mode A                                                %
%       beta_j = 1: Mode B                                                %
%      Z = N by J matrix of indicators                                    %
%      F = N by P matrix of sign-weighted composites of latents           %
%      B = P by P scheme matrix for connected latents                     %
%      A = P by J matrix of weights (~ loadings)                          %
%      W = J by nvlv matrix of weights                                    %
%      Gamma = N by P matrix of latents                                   %
% 2. The detailed description of the algorithm can be found in the        %
%    following paper:                                                     %
%      Hwang, H., Takane, Y. & Tenenhaus, A. An Alternative Estimation    %
%      Procedure for Partial Least Squares Path Modeling. Behaviormetrika %
%      42, 63–78 (2015). https://doi.org/10.2333/bhmk.42.63               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = size(W0,2);                    % num of LVs
[N, J] = size(z0);              % num of cases
z = zscore(z0)/sqrt(N-1);         % normalized data
%% Initialization 
W0=W0~=0; Nw=sum(sum(W0,1),2);
C0=W0';   Nc=sum(sum(C0,1),2);
B0=B0~=0; Nb=sum(sum(B0,1),2);
W = double(W0);
B = double(B0);
ind_Bdep=sum(B0,1)>0; Py = sum(ind_Bdep,2);
% Gamma = zeros(N,P);                     
% for j = 1:P
%     windex = find(W0(:,j));
%     z_j = z(:,windex);
%     [~, S V] = svd(z_j'*z_j);
%     Gamma(:,j) = z_j*V(:,1)/sqrt(S(1,1));    % first PC scores
% %     W(windex,j) = rand(length(windex),1);      % random starts
% %     Gamma(:,j) = z_j*W(windex,j);
% %     Gamma(:,j) = Gamma(:,j)/norm(Gamma(:,j));
% end

%vr = jcovcor(z0);
%S = vr.cor;
S=z'*z;
%W = W0./99;
W = W./repmat(sqrt(diag(W'*S*W))',J,1);    % same starts as in GSCA
Gamma = zeros(N,P);                     
for j = 1:P
    windex = find(W0(:,j));
    z_j = z(:,windex);
    Gamma(:,j) = z_j*W(windex,j);
    Gamma(:,j) = Gamma(:,j)/norm(Gamma(:,j));
end

 if scheme == 1              %centroid scheme
       corLV = corrcoef(Gamma);
       for j = 1:P
           bindex = find(B0(:,j));   % DV 
          if ~isempty(bindex)
             B(bindex,j) = sign(corLV(bindex,j));  
          end
       end 
       for j = 1:P
          bindex = find(B0(j,:));   % IV 
          if ~isempty(bindex)
             B(bindex,j) = sign(corLV(j,bindex));
          end
       end
 elseif scheme == 2         % factorial scheme
       corLV = corrcoef(Gamma);
       for j = 1:P
           bindex = find(B0(:,j));   % DV 
          if ~isempty(bindex)
             B(bindex,j) = corLV(bindex,j);  
          end
       end 
       for j = 1:P
          bindex = find(B0(j,:));   % IV 
         if ~isempty(bindex)
             B(bindex,j) = corLV(j,bindex);
          end
        end
 elseif scheme == 3       % path weighting scheme
       for j = 1:P
           bindex = find(B0(:,j));   % DV
           if ~isempty(bindex)
              gj = Gamma(:,bindex);
              B(bindex,j) = (gj'*gj)\gj'*Gamma(:,j);  
           end
       end    
       corLV = corrcoef(Gamma);
       for j = 1:P
           bindex = find(B0(j,:));   % IV
           if ~isempty(bindex)
               B(bindex,j) = corLV(j,bindex);
           end
       end
end
%LD = LD';
%WT;
%Path;
%B;
[est_W,est_C,est_B,it,Converge, Gamma] = ALS_BasicPLSPM(z,Gamma,W0,B0,W,B,modetype,scheme,Max_iter,Min_limit,N,J,P);
INI.iter = it;
INI.Converge=Converge;
INI.W = est_W;
INI.C = est_C;
INI.B = est_B;
INI.CVscore = Gamma*sqrt(N);

if N_Boot<100
   TABLE.W=[est_W(W0),NaN(Nw,5)];
   TABLE.C=[est_C(C0),NaN(Nc,5)];
   TABLE.B=[est_B(B0),NaN(Nb,5)];
   ETC.W_Boot=[];
   ETC.C_Boot=[];
   ETC.B_Boot=[];  
else
   W_Boot=zeros(Nw,N_Boot);
   C_Boot=zeros(Nc,N_Boot);
   B_Boot=zeros(Nb,N_Boot);
   if Flag_Parallel
       parfor b=1:N_Boot
           [Z_ib,~]=GC_Boot(z);
           [W_b,C_b,B_b,~,~]=ALS_BasicPLSPM(Z_ib,Gamma,W0,B0,W,B,modetype,scheme,Max_iter,Min_limit,N,J,P);           
           W_Boot(:,b)=W_b(W0);
           C_Boot(:,b)=C_b(C0);
           B_Boot(:,b)=B_b(B0);
       end
   else
       for b=1:N_Boot
           if rem(b,100)==1; fprintf("Bootstrapping %d\n", b); end
           [Z_ib,~]=GC_Boot(z);
           [W_b,C_b,B_b,~,~]=ALS_BasicPLSPM(Z_ib,Gamma,W0,B0,W,B,modetype,scheme,Max_iter,Min_limit,N,J,P);
           W_Boot(:,b)=W_b(W0);
           C_Boot(:,b)=C_b(C0);
           B_Boot(:,b)=B_b(B0);
       end
   end
% (5) Calculation of statistics
   alpha=.05;
   CI=[alpha/2,alpha,1-alpha,1-(alpha/2)];
   loc_CI=round(CI*(N_Boot-1))+1; % .025 .05 .95 .975
% basic statistics for parameter
   TABLE.W=para_stat(est_W(W0),W_Boot,loc_CI);
   TABLE.C=para_stat(est_C(C0),C_Boot,loc_CI); 
   if Py>0; TABLE.B=para_stat(est_B(B0),B_Boot,loc_CI); end
   ETC.W_Boot=W_Boot;
   ETC.C_Boot=C_Boot;
   ETC.B_Boot=B_Boot;    
end
end
function Table=para_stat(est_mt,boot_mt,CI_mp)
   boot_mt=sort(boot_mt,2);
   SE=std(boot_mt,0,2);
   Table=[est_mt,SE,boot_mt(:,CI_mp(1,1)),boot_mt(:,CI_mp(1,4))]; 
end
function [in_sample,out_sample,index,N_oob]=GC_Boot(Data)
   N=size(Data,1); 
   index=ceil(N*rand(N,1));
   in_sample=Data(index,:); 
   index_oob=(1:N)'; index_oob(index)=[];
   out_sample=Data(index_oob,:);
   N_oob=length(index_oob);
end
