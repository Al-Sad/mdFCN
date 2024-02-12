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
% This script generates Fig. 11 in the supplementary material. 

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC';    % HC or AED
fdr_thr = 0.025;    % 5% false discovery rate +/- threshold

%% Load results
X_pos = []; D_pos = []; R_pos = []; P_pos = [];
X_neg = []; D_neg = []; R_neg = []; P_neg = [];
for subset = 1:2
    load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'], ...
        'K_pos','K_neg','r_pos','r_neg','p_pos', 'p_neg','x_pos','x_neg');
    X_pos = cat(2,X_pos,x_pos); X_neg = cat(2,X_neg,x_neg);
    D_pos = cat(2,D_pos,K_pos); D_neg = cat(2,D_neg,K_neg);
    R_pos = cat(2,R_pos,r_pos); R_neg = cat(2,R_neg,r_neg);
    P_pos = cat(2,P_pos,p_pos); P_neg = cat(2,P_neg,p_neg);
end

%% Get averaged correlation
corr_m_pos = zeros(length(D_pos),length(D_pos{1}));
corr_m_neg = zeros(length(D_neg),length(D_neg{1}));
corr_s_pos = zeros(length(D_pos),length(D_pos{1}));
corr_s_neg = zeros(length(D_neg),length(D_neg{1}));
corr_p_pos = zeros(length(D_pos),length(D_pos{1}));
corr_p_neg = zeros(length(D_neg),length(D_neg{1}));
mask_pos = false(length(D_pos),length(D_pos{1}));
mask_neg = false(length(D_neg),length(D_neg{1}));
for i = 1:length(D_pos)
    mask_pos(i,:) = abs(D_pos{i}) >= fdr_thr;
    mask_neg(i,:) = abs(D_neg{i}) >= fdr_thr;
end
for i = 1:length(D_pos)
    corr_m_pos(i,:) = mean(R_pos{i},2,'omitnan')';  corr_m_pos(i,~mask_pos(i,:)) = nan;
    corr_m_neg(i,:) = mean(R_neg{i},2,'omitnan')';  corr_m_neg(i,~mask_neg(i,:)) = nan;
    corr_s_pos(i,:) = std(R_pos{i},0,2,'omitnan')';  corr_s_pos(i,~mask_pos(i,:)) = nan;
    corr_s_neg(i,:) = std(R_neg{i},0,2,'omitnan')';  corr_s_neg(i,~mask_neg(i,:)) = nan;
    corr_p_pos(i,:) = mean(P_pos{i},2,'omitnan')';  corr_p_pos(i,~mask_pos(i,:)) = nan;
    corr_p_neg(i,:) = mean(P_neg{i},2,'omitnan')';  corr_p_neg(i,~mask_neg(i,:)) = nan;
end

%% Plotting
ccc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; ...
    0.494 0.184 0.556; 0.929 0.694 0.125]; p1 = [];
legned_str = {'low-$\delta$','high-$\delta$','$\theta$','$\alpha$','$\beta$'};
XTickString{1} = '\begin{tabular}{c}Neurological\\Score (C1)\end{tabular}';
XTickString{2} = '\begin{tabular}{c}Neurological\\Score (C2)\end{tabular}';
XTickString{3} = '\begin{tabular}{c}Cognition\end{tabular}';
XTickString{4} = '\begin{tabular}{c}Language\\Comprehension\end{tabular}';
XTickString{5} = '\begin{tabular}{c}Language\\Production\end{tabular}';
XTickString{6} = '\begin{tabular}{c}Fine\\Motor Skills\end{tabular}';
XTickString{7} = '\begin{tabular}{c}Gross\\Motor Skills\end{tabular}';

figure('Color',[1,1,1],'Position',[25 25 1400 450]);
colororder({'k','k'}); yyaxis left; p1 = [];
hold on; grid on; box on; t_shift = 0.15; cnt = 1;
R_ylim = [-0.8 0.8]; P_ylim = [-0.08 0.08]; PN = 4;
t = zeros(1,numel(corr_p_pos));
for i = 1:size(corr_m_pos,1)
    for j = 1:size(corr_m_pos,2)
        t(cnt) = i + (j-1)*t_shift-2*t_shift;
        h = stem(t(cnt),corr_m_pos(i,j),'filled','Color',ccc(j,:), ...
            'LineStyle','-','LineWidth',12,'Marker','none');
        p1(j) = h(1);
        stem(t(cnt),corr_m_neg(i,j),'filled','Color',ccc(j,:), ...
            'LineStyle','-','LineWidth',12,'Marker','none');
        errorbar(t(cnt),corr_m_pos(i,j),0,2*corr_s_pos(i,j),'LineStyle','none', ...
            'Color',ccc(j,:),'LineWidth',3,'CapSize',12,'Marker','none');
        errorbar(t(cnt),corr_m_neg(i,j),2*corr_s_neg(i,j),0,'LineStyle','none', ...
            'Color',ccc(j,:),'LineWidth',3,'CapSize',12,'Marker','none');
        plot(repmat(0.5+(i-1),1,10),linspace(-1,1,10),'k-');
        cnt = cnt + 1;
    end
end
set(gca,'fontweight','bold','FontSize',22,'XGrid','off');
axis([0.5 size(corr_m_pos,1)+0.5 R_ylim]);
a = gca; a.XAxis.FontSize = 18;
set(gca,'xtick',1:size(corr_m_pos,1),'xticklabel',XTickString,'TickLabelInterpreter','latex', ...
    'XTickLabelRotationMode','manual','YTick',R_ylim(1):R_ylim(2)/PN:R_ylim(2));
ylabel('Averaged Correlation','Interpreter','latex','FontSize',28);
yyaxis right; hold on;
y = reshape(corr_p_pos',1,numel(corr_p_pos));
stem(t,y,'Color','k','LineStyle','-','LineWidth',6,'Marker','none');
y = -1.*reshape(corr_p_neg',1,numel(corr_p_neg));
h = stem(t,y,'Color','k','LineStyle','-','LineWidth',6,'Marker','none');
p1(end+1) = h(1);
ylabel('Averaged p-Value','Interpreter','latex','FontSize',30);
axis([0.5 size(corr_p_pos,1)+0.5 P_ylim(1) P_ylim(2)]);
set(gca,'xtick',1:size(corr_p_pos,1), ...
    'TickLabelInterpreter','latex','YTick',P_ylim(1):P_ylim(2)/PN:P_ylim(2),...
    'YTickLabel',[P_ylim(2):P_ylim(1)/PN:0 P_ylim(2)/PN:P_ylim(2)/PN:P_ylim(2)]);
legend(p1,[legned_str, 'p-value'],'Orientation','horizontal', ...
    'Location','southeast','interpreter','latex','fontsize',24);
title(grp_idx,'interpreter','latex','fontsize',32);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,['Results\Figures\Supp_Fig_11_correlation_' grp_idx],'-dpdf','-r400');
end