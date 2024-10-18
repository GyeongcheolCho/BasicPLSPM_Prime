%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for BasicPLSPM Prime package                               %
%   Author: Heungsun Hwang Gyeongcheol Cho                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%   - This code aims to illustrate how to use BasicPLSPM Prime package.   %
%   - The dataset is a replica of the ACSI data used in Cho & Hwang       %
%     (2024).                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Cho, G., Hwang, H. Generalized Structured Component Analysis      %
%         Accommodating Convex Components: A Knowledge-Based Multivariate %
%         Method with Interpretable Composite Indexes. Psychometrika 89,  %
%         241–266 (2024). https://doi.org/10.1007/s11336-023-09944-3      %  
%     * Hwang, H., Takane, Y. & Tenenhaus, A. An Alternative Estimation   %
%         Procedure for Partial Least Squares Path Modeling.              %
%         Behaviormetrika 42, 63–78 (2015).                               %
%         https://doi.org/10.2333/bhmk.42.63                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

help BasicPLSPM()

Data=readtable('ACSI_774_Replica.csv');
W0=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 ; ...
   0 0 0 1 1 1 0 0 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 1 1 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 0 0 1 1 1 0 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 1 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 0 1 1 ]';
C0=zeros(6,14);%W0';
B0=[0 1 1 1 0 0;...
    0 0 1 1 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 1;...
    0 0 0 0 0 1;...
    0 0 0 0 0 0];
N_Boot=1000;
Max_iter = 1000;
Min_limit = 10^(-6);
Flag_Parallel = false;
modetype=zeros(1,6); % 1 = reflective, 2 = formative
scheme=3; %  1 = centroid, 2 = factorial, 3 = path weighting
ind_sign = [1,4,7,9,12,13];
[INI,TABLE,ETC]=BasicPLSPM(Data{:,:}, W0, B0, modetype,scheme,ind_sign,N_Boot,Max_iter,Min_limit,Flag_Parallel);
INI
INI.Converge
INI.iter
INI.W
INI.C
INI.B
INI.CVscore
TABLE
TABLE.W
TABLE.C
TABLE.B
ETC