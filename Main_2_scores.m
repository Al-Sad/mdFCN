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
% It generate masks for the missing scores in each group and each sleep state.
% Note that the neonates raw EEG and neurocognitive scores are not
% supplied.

%% Initialization
clear; close all; clc;

%% Parameters
grp_idx = 'HC'; % HC or AED

%% Get the TA and AS masks
data_folder = ['Results\' grp_idx '\FC'];
mask_TA = false(length(dir([data_folder '*\*.mat'])),1);
mask_AS = false(length(dir([data_folder '*\*.mat'])),1);
for m = 1:length(dir([data_folder '*\*.mat']))
    load([data_folder '\Subj_' num2str(m) '.mat']);
    if(~isempty(dwPLI_TA))
        mask_TA(m) = true;
    else
        mask_TA(m) = false;
    end
    if(~isempty(dwPLI_AS))
        mask_AS(m) = true;
    else
        mask_AS(m) = false;
    end
end

%% Adjust scores according to the masks
x = readmatrix(['Dataset\Scores\Scores_' grp_idx '.xlsx']);
x(:,1) = [];
Age = x(:,3);
Scores = x(:,[1 2 4:8]);
mask_scores = ~isnan(Scores);
mask_score1 = all(mask_scores(:,1:2),2);
mask_score2 = all(mask_scores(:,3:7),2);
idx_TA  = [logical(mask_TA.*mask_score1) logical(mask_TA.*mask_score2)];
idx_AS  = [logical(mask_AS.*mask_score1) logical(mask_AS.*mask_score2)];
idx_all = [mask_score1 mask_score2];

%% Saving
save(['Results\' grp_idx '\scores.mat'],'Age','Scores','idx_TA','idx_AS','idx_all');
