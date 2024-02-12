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
% It verifies the mdFCN analysis pipeline. Note that the neonates raw EEG
% and neurocognitive scores are not supplied.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'AED'; % HC or AED
alpha   = 0.05;  % statistical significance level

%% Main
for test_num = 1:5
    for subset = 1:2
        % Get subject scores
        load(['Results\' grp_idx '\scores.mat'],'Scores','Age','idx_all');
        Scores = Scores(idx_all(:,subset),:);
        Age = Age(idx_all(:,subset),:);
        if(subset == 1)
            Scores = Scores(:,1:2);
        else
            Scores = Scores(:,3:end);
        end

        % Mask for the measures
        load('Head Model\FidelityOperator.mat');
        mask = logical(FidelityOperator);
        I = unfold_upper_matrix(mask);

        % Get values and perform correlations
        r = cell(1,size(Scores,2));
        p = cell(1,size(Scores,2));
        switch test_num
            case 1 % Static FC
                data_folder = ['Results\' grp_idx '\FC'];
                m = 1:length(dir([data_folder '*\*.mat']));
                m = m(idx_all(:,subset));
                T = [];
                g = [];
                for i = 1:length(m)
                    load([data_folder '\Subj_' num2str(m(i)) '.mat']);
                    if(~isempty(sdwPLI_TA))
                        T = cat(3,T,sdwPLI_TA(:,I));
                        g = cat(1,g,m(i));
                    end
                    if(~isempty(sdwPLI_AS))
                        T = cat(3,T,sdwPLI_AS(:,I));
                        g = cat(1,g,m(i));
                    end
                end
                x = [];
                K = unique(g);
                for k = 1:length(K)
                    x = cat(3,x,mean(T(:,:,g==K(k)),3));
                end
                x = permute(x,[1 3 2]);
                clear T dwPLI_TA dwPLI_AS sdwPLI_TA sdwPLI_AS g
                %%% Correlation
                Kp = zeros(size(Scores,2),size(x,1));
                Kn = zeros(size(Scores,2),size(x,1));
                x = tensor_correct_median(x);
                for i = 1:size(Scores,2)
                    [r{i}, p{i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
                    [Kp(i,:), Kn(i,:)] = corr_density(r(:,i), p(:,i), alpha);
                end

            case 2 % Dynamic FC
                data_folder = ['Results\' grp_idx '\FC'];
                m = 1:length(dir([data_folder '*\*.mat']));
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
                x = [];
                K = unique(g);
                for k = 1:length(K)
                    x = cat(2,x,mean(T(:,g==K(k),:),2));
                end
                clear T dwPLI_TA dwPLI_AS sdwPLI_TA sdwPLI_AS g
                %%% Correlation
                Kp = zeros(size(Scores,2),size(x,1));
                Kn = zeros(size(Scores,2),size(x,1));
                x = tensor_correct_median(x);
                for i = 1:size(Scores,2)
                    [r{i}, p{i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
                    [Kp(i,:), Kn(i,:)] = corr_density(r(:,i), p(:,i), alpha);
                end

            case 3 % Latent dynamic FC
                T = [];
                load(['Results\' grp_idx '\NMF\subset_' num2str(subset) '_order']);
                order = out_entropy.order;
                for i = 1:5
                    fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
                        num2str(i) '_order_' num2str(order(i))];
                    load(fname,'W','H','g');
                    T = cat(3,T,W*H);
                end
                x = [];
                K = unique(g);
                for k = 1:length(K)
                    x = cat(1,x,mean(T(g==K(k),:,:),1));
                end
                x = permute(x,[3 1 2]);
                clear T W H g
                %%% Correlation
                Kp = zeros(size(Scores,2),size(x,1));
                Kn = zeros(size(Scores,2),size(x,1));
                x = tensor_correct_median(x);
                for i = 1:size(Scores,2)
                    [r{i}, p{i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
                    [Kp(i,:), Kn(i,:)] = corr_density(r(:,i), p(:,i), alpha);
                end

            case 4 % Decomposed dynamic FC
                load(['Results\' grp_idx '\CPD\subset_' num2str(subset) ...
                    '_block_averaged_network']);
                x = sum(cat(4, block_avg_tensor{:}),4);
                %%% Correlation
                Kp = zeros(size(Scores,2),size(x,1));
                Kn = zeros(size(Scores,2),size(x,1));
                x = tensor_correct_median(x);
                for i = 1:size(Scores,2)
                    [r{i}, p{i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
                    [Kp(i,:), Kn(i,:)] = corr_density(r(:,i), p(:,i), alpha);
                end

            case 5 % Reconstructed dynamic FC
                load(['Results\' grp_idx '\Selection\subset_' num2str(subset) ...
                    '_selected'],'K_pos','K_neg');
                Kp = zeros(size(Scores,2),size(K_pos{1},2));
                Kn = zeros(size(Scores,2),size(K_neg{1},2));
                for i = 1:size(Scores,2)
                    Kp(i,:) = K_pos{i};
                    Kn(i,:) = -1.*K_neg{i};
                end
        end

        %%% Saving
        save(['Results\' grp_idx '\Verification\subset_' ...
            num2str(subset) '_test_' num2str(test_num)],'Kp','Kn');
    end
end