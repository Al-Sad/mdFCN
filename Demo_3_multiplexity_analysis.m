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
% This script generates the multiplexity analysis results in Fig. 5.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
fdr_thresh = 0.025;

%% Load all requirements
load('Head Model\FidelityOperator.mat');
mask = logical(FidelityOperator);
mask = unfold_upper_matrix(mask);

%% Load results
D_pos_hc  = []; X_pos_hc  = []; D_neg_hc  = []; X_neg_hc  = [];
D_pos_aed = []; X_pos_aed = []; D_neg_aed = []; X_neg_aed = [];
for subset = 1:2
    load(['Results\HC\Selection\subset_' num2str(subset) ...
        '_selected'],'K_pos','K_neg','x_pos_org','x_neg_org');
    D_pos_hc = cat(2,D_pos_hc,K_pos); D_neg_hc = cat(2,D_neg_hc,K_neg);
    X_pos_hc = cat(2,X_pos_hc,x_pos_org); X_neg_hc = cat(2,X_neg_hc,x_neg_org);
    load(['Results\AED\Selection\subset_' num2str(subset) ...
        '_selected'],'K_pos','K_neg','x_pos_org','x_neg_org');
    D_pos_aed = cat(2,D_pos_aed,K_pos); D_neg_aed = cat(2,D_neg_aed,K_neg);
    X_pos_aed = cat(2,X_pos_aed,x_pos_org); X_neg_aed = cat(2,X_neg_aed,x_neg_org);
end
mask_pos_hc = false(length(D_pos_hc), length(D_pos_hc{1}));
mask_neg_hc = false(length(D_pos_hc), length(D_pos_hc{1}));
mask_pos_aed = false(length(D_pos_aed), length(D_pos_aed{1}));
mask_neg_aed = false(length(D_pos_aed), length(D_pos_aed{1}));
for i = 1:length(D_pos_hc)
    mask_pos_hc(i,:) = abs(D_pos_hc{i}) >= fdr_thresh;
    mask_neg_hc(i,:) = abs(D_neg_hc{i}) >= fdr_thresh;
    mask_pos_aed(i,:) = abs(D_pos_aed{i}) >= fdr_thresh;
    mask_neg_aed(i,:) = abs(D_neg_aed{i}) >= fdr_thresh;
end

%% Multiplex network features
mx_pos_hc  = multiplex_measures(X_pos_hc,  mask_pos_hc,  mask);
mx_neg_hc  = multiplex_measures(X_neg_hc,  mask_neg_hc,  mask);
mx_pos_aed = multiplex_measures(X_pos_aed, mask_pos_aed, mask);
mx_neg_aed = multiplex_measures(X_neg_aed, mask_neg_aed, mask);
Oi_hc  = mean([mx_pos_hc.Oi; mx_neg_hc.Oi],'omitnan');
Pi_hc  = mean([mx_pos_hc.Pi; mx_neg_hc.Pi],'omitnan');
S_hc   = squeeze(mean([mx_pos_hc.S; mx_neg_hc.S],'omitnan'));
H_hc   = squeeze(mean([mx_pos_hc.H; mx_neg_hc.H],'omitnan'));
Oe_hc  = squeeze(mean([mx_pos_hc.Oe; mx_neg_hc.Oe],'omitnan'));
X_hc   = squeeze(mean([mx_pos_hc.X; mx_neg_hc.X],'omitnan'));
Oi_aed = mean([mx_pos_aed.Oi; mx_neg_aed.Oi],'omitnan');
Pi_aed = mean([mx_pos_aed.Pi; mx_neg_aed.Pi],'omitnan');
S_aed  = squeeze(mean([mx_pos_aed.S; mx_neg_aed.S],'omitnan'));
H_aed  = squeeze(mean([mx_pos_aed.H; mx_neg_aed.H],'omitnan'));
Oe_aed = squeeze(mean([mx_pos_aed.Oe; mx_neg_aed.Oe],'omitnan'));
X_aed  = squeeze(mean([mx_pos_aed.X; mx_neg_aed.X],'omitnan'));
range = [-0.1 1.1];
x = linspace(range(1),range(2),256);
y = linspace(range(1),range(2),256);
[X,Y] = meshgrid(x,y);
data = [Oi_hc; Pi_hc];
pdf_hc = mvnpdf([X(:) Y(:)], mean(data,2)',cov(data'));
pdf_hc = reshape(pdf_hc,256,256);
data = [Oi_aed; Pi_aed];
pdf_aed = mvnpdf([X(:) Y(:)], mean(data,2)',cov(data'));
pdf_aed = reshape(pdf_aed,256,256);

%% Plotting
M = 5;
figure('Color',[1,1,1],'Position',[100 100 750 650]);
subplot(M,M,[(1:M-1)*M+1 ((1:M-1)+1)*M-1]); hold on;
h1 = plot(nan,nan,'Marker','^','MarkerSize',20,'MarkerFaceColor', ...
    [0 0.447 0.741], 'MarkerEdgeColor','none','LineStyle','none');
h2 = plot(nan,nan,'Marker','o','MarkerSize',20,'MarkerFaceColor', ...
    [0.85 0.325 0.098], 'MarkerEdgeColor','none','LineStyle','none');
scatter(Oi_hc,Pi_hc,256,'Marker','^','MarkerFaceColor',[0 0.447 0.741], ...
    'MarkerEdgeColor','none','MarkerFacealpha',0.5);
scatter(Oi_aed,Pi_aed,256,'Marker','o','MarkerFaceColor',[0.85 0.325 0.098], ...
    'MarkerEdgeColor','none','MarkerFacealpha',0.5);
contour(x,y,pdf_hc,8,'color',[0 0.447 0.741],'linewidth',1);
contour(x,y,pdf_aed,8,'color',[0.85 0.325 0.098],'linewidth',1);
legend([h1(1) h2(1)],'HC','AED','Location','southeast', ...
    'fontsize',22,'Interpreter','latex');
set(gca,'FontSize',18,'FontWeight','bold');
xlabel('Overlapping Degree','Interpreter','latex','FontSize',24);
ylabel('Participation Coefficient','Interpreter','latex','FontSize',24);
axis([0 0.5 0 1]); grid on; box on;
subplot(M,M,(2:M)*M);
histogram(Pi_hc,16,'Orientation','horizontal'); hold on;
histogram(Pi_aed,16,'Orientation','horizontal'); ylim([0 1]); box off;
set(gca,'XColor','w','TickDir','out','YTickLabel','');
subplot(M,M,1:M-1); hold on;
histogram(Oi_hc,16); histogram(Oi_aed,16); xlim([0 0.5]); box off;
set(gca,'YColor','w','TickDir','out','XTickLabel','');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

band_str = {'low-\delta','high-\delta','\theta','\alpha','\beta'};
figure('Color',[1,1,1],'Position',[10 10 450 400]);
h = heatmap(band_str,band_str,round(S_hc,2), ...
    'MissingDataColor',[1 1 1],'fontsize',20,'MissingDataLabel','n.s.');
colormap parula; caxis([0.2 0.99]);
h.YDisplayData = flipud(h.YDisplayData); title('HC');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[10 10 450 400]);
h = heatmap(band_str,band_str,round(S_aed,2), ...
    'MissingDataColor',[1 1 1],'fontsize',20,'MissingDataLabel','n.s.');
colormap parula; caxis([0.2 0.99]);
h.YDisplayData = flipud(h.YDisplayData); title('AED');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[10 10 450 400]);
h = heatmap(band_str,band_str,round(H_hc,2), ...
    'MissingDataColor',[1 1 1],'fontsize',20,'MissingDataLabel','n.s.');
colormap parula; caxis([0 0.69]);
h.YDisplayData = flipud(h.YDisplayData); 
title('HC');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[10 10 450 400]);
h = heatmap(band_str,band_str,round(H_aed,2), ...
    'MissingDataColor',[1 1 1],'fontsize',20,'MissingDataLabel','n.s.');
colormap parula; caxis([0 0.69]);
h.YDisplayData = flipud(h.YDisplayData); title('AED');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

st = 0.14;
figure('Color',[1,1,1],'Position',[10 10 350 550]);
bar(1,[Oe_hc Oe_aed],0.5); grid on; hold on;
errorbar(1-st,Oe_hc,0.1,2.*std([mx_pos_hc.Oe; mx_neg_hc.Oe],'omitnan'), ...
    'vertical','color',[0 0.447 0.741],'CapSize',22,'LineWidth',4);
errorbar(1+st,Oe_aed,0.1,2.*std([mx_pos_aed.Oe; mx_neg_aed.Oe],'omitnan'), ...
    'vertical','color',[0.85 0.325 0.098],'CapSize',22,'LineWidth',4);
ylim([0 1.05]);
legend('HC','AED','fontsize',28,'Interpreter','latex', ...
    'orientation','horizontal','location','northwest');
set(gca,'XTickLabel','','fontsize',24,'fontweight','bold');
title({'Total Edge','Overlap Ratio'},'Interpreter','latex');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[10 10 350 550]);
bar(1,[X_hc X_aed],0.5); grid on; hold on;
errorbar(1-st,X_hc,0.1,2.*std([mx_pos_hc.X; mx_neg_hc.X],'omitnan'), ...
    'vertical','color',[0 0.447 0.741],'CapSize',22,'LineWidth',4);
errorbar(1+st,X_aed,0.1,2.*std([mx_pos_aed.X; mx_neg_aed.X],'omitnan'), ...
    'vertical','color',[0.85 0.325 0.098],'CapSize',22,'LineWidth',4);
ylim([0 0.38]);
legend('HC','AED','fontsize',28,'Interpreter','latex', ...
    'orientation','horizontal','location','northwest');
set(gca,'XTickLabel','','fontsize',24,'fontweight','bold');
title({'Overall Edge','Intersection Index'},'Interpreter','latex');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Fig_5_spatial_features','-dpdf','-r400');
    print(2,'Results\Figures\Fig_5_pairwise_multiplexity_HC','-dpdf','-r400');
    print(3,'Results\Figures\Fig_5_pairwise_multiplexity_AED','-dpdf','-r400');
    print(4,'Results\Figures\Fig_5_hamming_distance_HC','-dpdf','-r400');
    print(5,'Results\Figures\Fig_5_hamming_distance_AED','-dpdf','-r400');
    print(6,'Results\Figures\Fig_5_total_edge_overlap_ratio','-dpdf','-r400');
    print(7,'Results\Figures\Fig_5_overall_edge_intersection_index','-dpdf','-r400');
end