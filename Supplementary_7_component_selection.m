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
% This script generates Fig. 13 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'AED';    % HC or AED

%% Main
I = cell(1,2);
for j = 1:2 % correlation directions
    II = cell(1,2);
    for subset = 1:2
        load(['Results\' grp_idx '\CPD\subset_' num2str(subset) '_block_averaged_network']);
        load(['Results\' grp_idx '\Selection\subset_' num2str(subset) ...
            '_selected'],'I_pos','I_neg');
        if(j == 1)
            Idir = I_pos;
        else
            Idir = I_neg;
        end
        Q = length(block_avg_tensor);
        II{subset} = zeros(Q,size(Idir,1)*size(Idir,2)+size(Idir,2));
        cnt = 1;
        for i = 1:size(Idir,2)
            for f = 1:size(Idir,1)
                if(~isnan(Idir(f,i)))
                    II{subset}(:,cnt) = fliplr(int2bit(Idir(f,i),Q)')';
                else
                    II{subset}(:,cnt) = zeros(Q,1);
                end
                cnt = cnt + 1;
            end
            II{subset}(:,cnt) = nan(Q,1);
            cnt = cnt + 1;
        end
    end
    II{2} = II{2}(:,1:end-1);
    if(size(II{1},1) > size(II{2},1))
        I{j} = [II{1} [II{2}; nan(size(II{1},1)-size(II{2},1),size(II{2},2))]];
    elseif(size(II{2},1) > size(II{1},1))
        I{j} = [[II{1}; nan(size(II{2},1)-size(II{1},1),size(II{1},2))] II{2}];
    else
        I{j} = [II{1} II{2}];
    end
end

%% Plotting
for j = 1:2 % correlation directions
    figure('Color',[1,1,1],'Position',[10 10 1400 500]);
    heatmap(I{j},'MissingDataColor',[1 1 1],'fontsize',14,'XDisplayLabels', ...
        repelem({''},1,size(I{j},2)),'GridVisible','off'); colorbar off;
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,['Results\Figures\Supp_Fig_13_selection_pos_' grp_idx],'-dpdf','-r400');
    print(2,['Results\Figures\Supp_Fig_13_selection_neg_' grp_idx],'-dpdf','-r400');
end