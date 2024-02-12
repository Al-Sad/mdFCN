% Copyright (c) 2024 Mohammad Al-Sa'd
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Email: mohammad.al-sad@helsinki.fi, alsad.mohamed@gmail.com
%
% The following reference should be cited whenever this script is used:
% Al‐Sa'd, M., Vanhatalo, S. and Tokariev, A., 2024. Multiplex dynamic
% networks in the newborn brain disclose latent links with neurobehavioral
% phenotypes. Human Brain Mapping, 45(2), https://doi.org/10.1002/hbm.26610
%
% Last Modification: 12-February-2024
%
% Description:
% It decomposes the latent networks by CPD and selects the best model order.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC'; % HC or AED
subset  = 1;    % 1 or 2
Nb      = 5;    % Number of bands
Qmax    = 40;   % Maximum CPD order
N_rep   = 20;   % CPD repetitions
MaxIter = 1e3;  % Maximum number of iterations
Tol     = 1e-8; % CPD tolerance
per     = 0.05; % Order selection percentage

%% Get NMF extracted latent network
T = [];
load(['Results\' grp_idx '\NMF\subset_' num2str(subset) '_order']);
order = out_entropy.order;
for i = 1:Nb
    fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
        num2str(i) '_order_' num2str(order(i))];
    load(fname,'W','H','g');
    x = W*H;
    T = cat(3,T,x);
end
Ent = entropy(T);
T = permute(T,[1 3 2]);

%% CPD-based latent network decomposition
options.TolFun  = Tol;
options.TolX    = Tol;
options.MaxIter = MaxIter;
for q = 1:Qmax
    U0 = cell(N_rep,3);
    for n = 1:N_rep
        U0{n,1} = rand(size(T,1),q);
        U0{n,2} = rand(size(T,2),q);
        U0{n,3} = rand(size(T,3),q);
    end
    trial_error = zeros(1,N_rep);
    trial_entropy = zeros(1,N_rep);
    U_hat = cell(1,N_rep);
    for n = 1:N_rep
        disp(['Order ' num2str(q) ', Trial ' num2str(n)]);
        U_hat{n} = cpd_positive_als(T,U0(n,:),options);
        trial_error(1,n) = frobcpdres(T,U_hat{n})/frob(T);
        trial_entropy(1,n) = entropy(cpdgen(U_hat{n}));
    end
    [~,In] = min(trial_error);
    U_hat = U_hat{In};
    %%% Saving
    save(['Results\' grp_idx '\CPD\subset_' num2str(subset) '_order_' num2str(q)],...
        'U_hat','g','trial_error','trial_entropy','Ent'); 
end

%% Model automatic order selection
Entropy = zeros(1, Qmax, N_rep);
Error   = zeros(1, Qmax, N_rep);
for q = 1:Qmax
    fname = ['Results\' grp_idx '\CPD\subset_' num2str(subset)...
        '_order_' num2str(q)];
    load(fname,'trial_entropy','trial_error','Ent');
    Entropy(1,q,:) = trial_entropy./Ent;
    Error(1,q,:)   = trial_error;
end
out_error   = order_selection(Error, per, 0);
out_entropy = order_selection(Entropy, per, 1);
%%% Saving
fname = ['Results\' grp_idx '\CPD\subset_' num2str(subset) '_order'];
save(fname,'out_error','out_entropy');
