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
% It extracts dynamic FC latent networks by NMF and selects the best model
% order.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC'; % HC or AED
subset  = 1;    % 1 or 2 for short- or long-term scores
N_rep   = 10;   % NMF repetitions
Pmax    = 40;   % Maximum NMF order
per     = 0.05; % Order selection percentage

%% Mask for the measures
load('Head Model\FidelityOperator.mat');
mask = logical(FidelityOperator);
I = unfold_upper_matrix(mask);

%% Form group-level tensor
data_folder = ['Results\' grp_idx '\FC'];
m = 1:length(dir([data_folder '*\*.mat']));
load(['Results\' grp_idx '\scores.mat'],'idx_all');
m = m(idx_all(:,subset));
T = [];
g = [];
for i = 1:length(m)
    load([data_folder '\Subj_' num2str(m(i)) '.mat']);
    if(~isempty(dwPLI_TA))
        T = cat(2,T,dwPLI_TA(:,:,I));
        g = cat(1,g,m(i));
    end
    if(~isempty(dwPLI_AS))
        T = cat(2,T,dwPLI_AS(:,:,I));
        g = cat(1,g,m(i));
    end
end
Ent = zeros(1,size(T,1));
for i = 1:size(T,1)
    Ent(i) = entropy(squeeze(T(i,:,:)));
end

%% NMF-based latent network extraction
for i = 1:size(T,1)
    x = squeeze(T(i,:,:));
    for p = 1:Pmax
        W = cell(1,N_rep);
        H = cell(1,N_rep);
        trial_error = zeros(1,N_rep);
        trial_entropy = zeros(1,N_rep);
        for n = 1:N_rep
            fprintf('Band %d, Order %d, Trial %d\n',i,p,n);
            [W{n},H{n}] = nmfnnls_mod(x,p);
            trial_error(n) = norm(x-(W{n}*H{n}),'fro')/norm(x,'fro');
            trial_entropy(n) = entropy(W{n}*H{n});
        end
        [~,In] = min(trial_error);
        W = W{In}; H = H{In};
        %%% Saving
        fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
            num2str(i) '_order_' num2str(p)];
        save(fname,'W','H','g','Ent','trial_error','trial_entropy');
    end
end

%% Model automatic order selection
Entropy = zeros(size(T,1), Pmax, N_rep);
Error   = zeros(size(T,1), Pmax, N_rep);
for i = 1:size(T,1)
    for p = 1:Pmax
        fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
            num2str(i) '_order_' num2str(p)];
        load(fname,'trial_entropy','trial_error','Ent');
        Entropy(i,p,:) = trial_entropy./Ent(i);
        Error(i,p,:)   = trial_error;
    end
end
out_error   = order_selection(Error, per, 0);
out_entropy = order_selection(Entropy, per, 1);
%%% Saving
fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_order'];
save(fname,'out_error','out_entropy');