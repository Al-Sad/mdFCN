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
% It reconstructs the decomposed networks and generates null
% distributions for all other reconstruction configurations.
% Note that the neonates raw EEG and neurocognitive scores are not
% supplied.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC';  % HC or AED
subset  = 1;     % 1 or 2
alpha   = 0.05;  % Statistical significance level
Kperm   = 10000; % Number of null permutations

%% Get scores
load(['Results\' grp_idx '\scores.mat'],'Scores','Age','idx_all');
Scores = Scores(idx_all(:,subset),:);
Age = Age(idx_all(:,subset),:);
if(subset == 1)
    Scores = Scores(:,1:2);
else
    Scores = Scores(:,3:end);
end

%% Temporal block-averaging
load(['Results\HC\CPD\subset_' num2str(subset) '_order']);
order_hc = out_entropy.order;
load(['Results\AED\CPD\subset_' num2str(subset) '_order']);
order_aed = out_entropy.order;
order = max([order_hc order_aed]);
load(['Results\' grp_idx '\CPD\subset_' num2str(subset) ...
    '_order_' num2str(order)],'U_hat','g')
U_hat = cpd_network_time_avg(U_hat, 179);
U = cell(1);
for i = 1:size(U_hat{1},2)
    K = unique(g);
    for k = 1:length(K)
        U{1}(k,i) = mean(U_hat{1}(g==K(k),i));
    end
end
U_hat = {U{1}, U_hat{2}, U_hat{3}};
block_avg_tensor = cell(1,size(U_hat{1},2));
for i = 1:size(U_hat{1},2)
    block_avg_tensor{i} = cpdgen({U_hat{2}(:,i), U_hat{1}(:,i),U_hat{3}(:,i)});
end
%%% Saving
save(['Results\' grp_idx '\CPD\subset_' num2str(subset) ...
    '_block_averaged_network'],'block_avg_tensor');

%% Compute correlations for every reconstruction configuration
load(['Results\' grp_idx '\CPD\subset_' num2str(subset) '_block_averaged_network']);
N = 2^length(block_avg_tensor)-1;
r = cell(N,size(Scores,2));
p = cell(N,size(Scores,2));
for n = 1:N
    disp(100*n/N);
    I = fliplr(int2bit(n,length(block_avg_tensor))');
    temp = block_avg_tensor(logical(I));
    x = sum(cat(4, temp{:}),4);
    x = tensor_correct_median(x);
    for i = 1:size(Scores,2)
        [r{n,i}, p{n,i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
    end
end
%%% Saving
save(['Results\' grp_idx '\Selection\subset_' num2str(subset) ...
    '_correlations'],'r','p','-v7.3');

%% Compute fitness function for every reconstruction configuration
load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_correlations']);
[Nb, Ns, Nc] = size(block_avg_tensor{1});
Fit_pos = cell(Nb, size(r,2));
Fit_neg = cell(Nb, size(r,2));
% loop through the neurocognitive scores
for i = 1:size(r,2)
    [Kp, Kn] = corr_density(r(:,i), p(:,i), alpha);
    % loop through the frequency bands
    for f = 1:Nb
        IIp = find(Kp(:,f) > 0);
        IIn = find(Kn(:,f) > 0);
        % Compute fitness for positive correlation densities
        if(~isempty(IIp))
            r2pos = zeros(1,length(IIp));
            kpos  = zeros(1,length(IIp));
            rpos  = zeros(1,length(IIp));
            % loop through reconstruction configurations
            for j = 1:length(IIp)
                kpos(j) = Kp(IIp(j),f);
                I = fliplr(int2bit(IIp(j),length(block_avg_tensor))');
                pvals = p{IIp(j),i}(f,:); pvals(isnan(pvals)) = 1;
                It = ~isnan(pvals) & make_FDR(pvals,alpha);
                II = (r{IIp(j),i}(f,:) > 0) & It; % Positive correlations
                temp = r{IIp(j),i}(f,:); temp(1,~II) = nan;
                rpos(j) = mean(temp,'omitnan');
                T = block_avg_tensor(logical(I));
                T = sum(cat(4,T{:}),4);
                T = tensor_correct_median(T);
                T = squeeze(T(f,:,:));
                T(:,~II) = nan;
                T = mean(T,2,'omitnan');
                if(var(T,'omitnan') < 1e-10)
                    r2pos(j) = 0;
                else
                    [P, S] = polyfit(T,Scores(:,i),1);
                    yfit = polyval(P, T, S);
                    SStot = sum((Scores(:,i)-mean(Scores(:,i))).^2);
                    SSres = sum((Scores(:,i)-yfit).^2);
                    r2pos(j) = 1-SSres/SStot;
                end
            end
            Fit_pos{f,i} = (kpos.*r2pos.*rpos).^(1/3);
        else
            Fit_pos{f,i} = 0;
        end
        % Compute fitness for negative correlation densities
        if(~isempty(IIn))
            r2neg = zeros(1,length(IIn));
            kneg  = zeros(1,length(IIn));
            rneg  = zeros(1,length(IIn));
            % loop through reconstruction configurations
            for j = 1:length(IIn)
                kneg(j) = Kn(IIn(j),f);
                I = fliplr(int2bit(IIn(j),length(block_avg_tensor))');
                pvals = p{IIn(j),i}(f,:); pvals(isnan(pvals)) = 1;
                It = ~isnan(pvals) & make_FDR(pvals,alpha);
                II = (r{IIn(j),i}(f,:) < 0) & It; % Negative correlations
                temp = r{IIn(j),i}(f,:); temp(1,~II) = nan;
                rneg(j) = abs(mean(temp,'omitnan'));
                T = block_avg_tensor(logical(I));
                T = sum(cat(4,T{:}),4);
                T = tensor_correct_median(T);
                T = squeeze(T(f,:,:));
                T(:,~II) = nan;
                T = mean(T,2,'omitnan');
                if(var(T,'omitnan')<1e-10)
                    r2neg(j) = 0;
                else
                    [P, S] = polyfit(T,Scores(:,i),1);
                    yfit = polyval(P, T, S);
                    SStot = sum((Scores(:,i)-mean(Scores(:,i))).^2);
                    SSres = sum((Scores(:,i)-yfit).^2);
                    r2neg(j) = 1-SSres/SStot;
                end
            end
            Fit_neg{f,i} = (kneg.*r2neg.*rneg).^(1/3);
        else
            Fit_neg{f,i} = 0;
        end
    end
end

%% Component selection
I_pos = zeros(Nb, size(r,2));
K_pos = cell(1, size(r,2));
r_pos = cell(1, size(r,2));
p_pos = cell(1, size(r,2));
x_pos = cell(1, size(r,2));
I_neg = zeros(Nb, size(r,2));
K_neg = cell(1, size(r,2));
r_neg = cell(1, size(r,2));
p_neg = cell(1, size(r,2));
x_neg = cell(1, size(r,2));
x_pos_org = cell(1, size(r,2));
x_neg_org = cell(1, size(r,2));
for i = 1:size(r,2)
    [Kp, Kn] = corr_density(r(:,i), p(:,i), alpha);
    % loop through the frequency bands
    for f = 1:Nb
        [Vp, Ip] = max(Fit_pos{f,i});
        [Vn, In] = max(Fit_neg{f,i});
        IIp = find(Kp(:,f) > 0);
        IIn = find(Kn(:,f) > 0);
        if(Vp > 0)
            I_pos(f,i) = IIp(Ip);
            I = fliplr(int2bit(I_pos(f,i),length(block_avg_tensor))');
            T = block_avg_tensor(logical(I));
            T = sum(cat(4,T{:}),4);
            x_pos_org{i}(f,:,:) = T(f,:,:);
            T = tensor_correct_median(T);
            x_pos{i}(f,:,:) = T(f,:,:);
            K_pos{i}(1,f) = Kp(I_pos(f,i),f);
            pvals = p{I_pos(f,i),i}(f,:); pvals(isnan(pvals)) = 1;
            It = ~isnan(pvals) & make_FDR(pvals,alpha);
            II = (r{I_pos(f,i),i}(f,:) > 0) & It;
            r_pos{i}(f,:) = r{I_pos(f,i),i}(f,:);
            p_pos{i}(f,:) = p{I_pos(f,i),i}(f,:);
            r_pos{i}(f,~II) = nan;
            p_pos{i}(f,~II) = nan;
            x_pos{i}(f,:,~II) = nan;
            x_pos_org{i}(f,:,~II) = nan;
        else
            I_pos(f,i) = nan;
            x_pos{i}(f,:,:) = nan(1,Ns,Nc);
            x_pos_org{i}(f,:,:) = nan(1,Ns,Nc);
            K_pos{i}(1,f) = 0;
            r_pos{i}(f,:) = nan(1,Nc);
            p_pos{i}(f,:) = nan(1,Nc);
        end
        if(Vn > 0)
            I_neg(f,i) = IIn(In);
            I = fliplr(int2bit(I_neg(f,i),length(block_avg_tensor))');
            T = block_avg_tensor(logical(I));
            T = sum(cat(4,T{:}),4);
            x_neg_org{i}(f,:,:) = T(f,:,:);
            T = tensor_correct_median(T);
            x_neg{i}(f,:,:) = T(f,:,:);
            K_neg{i}(1,f) = -1.*Kn(I_neg(f,i),f);
            pvals = p{I_neg(f,i),i}(f,:); pvals(isnan(pvals)) = 1;
            It = ~isnan(pvals) & make_FDR(pvals,alpha);
            II = (r{I_neg(f,i),i}(f,:) < 0) & It;
            r_neg{i}(f,:) = r{I_neg(f,i),i}(f,:);
            p_neg{i}(f,:) = p{I_neg(f,i),i}(f,:);
            r_neg{i}(f,~II) = nan;
            p_neg{i}(f,~II) = nan;
            x_neg{i}(f,:,~II) = nan;
            x_neg_org{i}(f,:,~II) = nan;
        else
            I_neg(f,i) = nan;
            x_neg{i}(f,:,:) = nan(1,Ns,Nc);
            x_neg_org{i}(f,:,:) = nan(1,Ns,Nc);
            K_neg{i}(1,f) = 0;
            r_neg{i}(f,:) = nan(1,Nc);
            p_neg{i}(f,:) = nan(1,Nc);
        end
    end
end
save(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'], ...
    'K_pos','K_neg','I_pos','I_neg','r_pos','p_pos', ...
    'r_neg','p_neg','x_pos','x_neg','x_pos_org','x_neg_org');

%% Null distribution
load(['Results\' grp_idx '\Selection\subset_' num2str(subset) ...
    '_selected'],'x_pos','x_neg');
Nc = size(x_pos{1},3);
r_pos = cell(Kperm,size(Scores,2));
p_pos = cell(Kperm,size(Scores,2));
r_neg = cell(Kperm,size(Scores,2));
p_neg = cell(Kperm,size(Scores,2));
for k = 1:Kperm
    disp(['Trial : ' num2str(k)]);
    for i = 1:size(Scores,2)
        x = x_pos{i};
        y = x_neg{i};
        for a = 1:size(x,2)
            x(:,a,:) = x(:,a,randperm(Nc));
        end
        for a = 1:size(y,2)
            y(:,a,:) = y(:,a,randperm(Nc));
        end
        [r_pos{k,i}, p_pos{k,i}] = partialcorr_score_tensor(x, Scores(:,i), Age);
        [r_neg{k,i}, p_neg{k,i}] = partialcorr_score_tensor(y, Scores(:,i), Age);
    end
end
%%% Saving
save(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_permutations'], ...
    'r_pos','p_pos','r_neg','p_neg','-v7.3');
