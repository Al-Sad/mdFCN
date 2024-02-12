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
% This script generates plots of Fig. 1 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Load dFC tensor for subject 2
load('Results\HC\FC\Subj_2.mat','dwPLI_TA');
T  = dwPLI_TA;
T  = permute(T,[2 3 1]);
T1 = fold_upper_matrix(squeeze(T(1,:,1)));
T2 = fold_upper_matrix(squeeze(T(2,:,1)));
T3 = fold_upper_matrix(squeeze(T(179,:,1)));

%% Plotting
[N,M,P] = size(T);
figure('Color',[1 1 1],'Position',[100 100 550 450]);
imagesc(T1); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 550 450]);
imagesc(T2); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 550 450]);
imagesc(T3); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 450 600]);
imagesc(squeeze(T(:,:,1))'); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 450 500]); hold on
for i = 1:P
    h = surface(i*ones(1,M),linspace(1-20*i,N-20*i,N),...
        repmat(linspace(M+100*i,1+100*i,M),N,1),repmat(M:-1:1,N,1));
    h.CData = T(:,:,i); caxis([0 1]);
    h.EdgeColor = 'none';
    h.FaceColor = 'interp';
end
axis off; view([-90 0]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Supp_Fig_1_instantaneous_connectivity_1','-dpdf','-r400');
    print(2,'Results\Figures\Supp_Fig_1_instantaneous_connectivity_2','-dpdf','-r400');
    print(3,'Results\Figures\Supp_Fig_1_instantaneous_connectivity_3','-dpdf','-r400');
    print(4,'Results\Figures\Supp_Fig_1_unfolded_connectivity','-dpdf','-r400');
    print(5,'Results\Figures\Supp_Fig_1_subject_mdFCN','-dpdf','-r400');
end