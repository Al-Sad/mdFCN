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
% This script generates the mdFCN verification results in Fig. 6.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx   = 'AED';    % HC or AED
freq_band = 1;        % 1-5
focus_idx = 3;        % Expand results for this score index (1-7)
alpha     = 0.05;     % Statistical significance level
fdr_thr   = 0.025;    % 5% false discovery rate +/- threshold

%% Get raw group-level tensor
if(focus_idx > 2), subset = 2; else, subset = 1; end
load('Head Model\FidelityOperator.mat');
mask = logical(FidelityOperator);
mask = unfold_upper_matrix(mask);
data_folder = ['Results\' grp_idx '\FC'];
m = 1:length(dir([data_folder '*\*.mat']));
load(['Results\' grp_idx '\scores.mat'],'idx_all');
m = m(idx_all(:,subset));
T_group = [];
for i = 1:length(m)
    load([data_folder '\Subj_' num2str(m(i)) '.mat']);
    if(~isempty(dwPLI_TA))
        T_group = cat(1,T_group,squeeze(dwPLI_TA(freq_band,:,mask)));
    end
    if(~isempty(dwPLI_AS))
        T_group = cat(1,T_group,squeeze(dwPLI_AS(freq_band,:,mask)));
    end
end

%% Get NMF-latent network
load(['Results\' grp_idx '\NMF\subset_' num2str(subset) '_order']);
nmf_order = out_entropy.order;
fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
    num2str(freq_band) '_order_' num2str(nmf_order(freq_band))];
load(fname,'W','H');
T_nmf = W*H;

%% Get CPD network
load(['Results\HC\CPD\subset_' num2str(subset) '_order']);
order_hc = out_entropy.order;
load(['Results\AED\CPD\subset_' num2str(subset) '_order']);
order_aed = out_entropy.order;
cpd_order = max([order_hc order_aed]);
load(['Results\' grp_idx '\CPD\subset_' num2str(subset) ...
    '_order_' num2str(cpd_order)],'U_hat');
T_cpd = cpdgen(U_hat);
T_cpd = squeeze(T_cpd(:,freq_band,:));

%% Get Reconstructed network
load(['Results\' grp_idx '\Selection\subset_1_selected'],'I_pos'); temp1 = I_pos;
load(['Results\' grp_idx '\Selection\subset_2_selected'],'I_pos'); temp2 = I_pos;
I_pos = [temp1 temp2];
I = logical(fliplr(int2bit(I_pos(freq_band,focus_idx),size(U_hat{1},2))'));
U_hat_select = {U_hat{1}(:,I) U_hat{2}(:,I) U_hat{3}(:,I)};
T_select = cpdgen(U_hat_select);
T_select = squeeze(T_select(:,freq_band,:));

%% Verification results
Km_pos = zeros(1,5);
Km_neg = zeros(1,5);
Ks_pos = zeros(1,5);
Ks_neg = zeros(1,5);
for test_num = 1:5
    load(['Results\' grp_idx '\Verification\subset_1_test_' num2str(test_num)]);
    temp1 = Kp;  temp2 = Kn;
    load(['Results\' grp_idx '\Verification\subset_2_test_' num2str(test_num)]);
    KKp = 100.*[temp1; Kp]; KKn = -100.*[temp2; Kn];
    Km_pos(1,test_num) = mean(KKp(focus_idx,:));
    Km_neg(1,test_num) = mean(KKn(focus_idx,:));
    Ks_pos(1,test_num) = std(KKp(focus_idx,:));
    Ks_neg(1,test_num) = std(KKn(focus_idx,:));
end

%% Null distribution
load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_permutations']);
if(subset > 1)
    p_pos = p_pos(:,focus_idx-2);
    p_neg = p_neg(:,focus_idx-2);
    r_pos = r_pos(:,focus_idx-2);
    r_neg = r_neg(:,focus_idx-2);
else
    p_pos = p_pos(:,focus_idx);
    p_neg = p_neg(:,focus_idx);
    r_pos = r_pos(:,focus_idx);
    r_neg = r_neg(:,focus_idx);
end
K_pos_null = zeros(size(r_pos,1),size(r_pos{1},1));
K_neg_null = zeros(size(r_pos,1),size(r_pos{1},1));
for i = 1:size(r_pos,1)
    [K_pos_null(i,:), ~] = corr_density(r_pos(i), p_pos(i), alpha);
    [~, K_neg_null(i,:)] = corr_density(r_neg(i), p_neg(i), alpha);
end

%% Plotting
x1 = 2000; x2 = 3000;
y1 = 150; y2 = 450;
figure('Color',[1 1 1],'Position',[100 100 650 450]);
imagesc(T_group'); axis xy; caxis([0 1]); hold on; axis off;
plot(repelem(x1,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(repelem(x2,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y1,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y2,100),'k-','linewidth',4);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 650 450]);
imagesc(T_nmf'); axis xy; caxis([0.3 0.75]); hold on; axis off;
plot(repelem(x1,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(repelem(x2,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y1,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y2,100),'k-','linewidth',4);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 650 450]);
imagesc(T_cpd'); axis xy; caxis([0.3 0.75]); hold on; axis off;
plot(repelem(x1,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(repelem(x2,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y1,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y2,100),'k-','linewidth',4);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 650 450]);
imagesc(T_select'); axis xy; caxis([0 0.4]); hold on; axis off;
plot(repelem(x1,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(repelem(x2,100),linspace(y1,y2,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y1,100),'k-','linewidth',4);
plot(linspace(x1,x2,100),repelem(y2,100),'k-','linewidth',4);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 350 450]);
imagesc(T_group(x1:x2,y1:y2)'); axis xy; caxis([0 1]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 350 450]);
imagesc(T_nmf(x1:x2,y1:y2)'); axis xy; caxis([0.3 0.75]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 350 450]);
imagesc(T_cpd(x1:x2,y1:y2)'); axis xy; caxis([0.3 0.75]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 350 450]);
imagesc(T_select(x1:x2,y1:y2)'); axis xy; caxis([0 0.4]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));


figure('Color',[1,1,1],'Position',[25 25 750 350]);
stem(1:length(Km_pos),Km_pos,'filled','Color',[0.85,0.33,0.1], ...
    'LineStyle','-','LineWidth',12,'Marker','none'); hold on
stem(1:length(Km_neg),Km_neg,'filled','Color',[0,0.45,0.74], ...
    'LineStyle','-','LineWidth',12,'Marker','none'); grid on;
errorbar(1:length(Km_pos),Km_pos,zeros(1,length(Km_pos)),2.*Ks_pos, ...
    'LineStyle','none','Color',[0.85,0.33,0.1],'LineWidth',3,'CapSize',12,'Marker','none');
errorbar(1:length(Km_neg),Km_neg,2.*Ks_neg,zeros(1,length(Km_neg)), ...
    'LineStyle','none','Color',[0,0.45,0.74],'LineWidth',3,'CapSize',12,'Marker','none');
for j = 1:length(Km_pos)
    plot(repmat(0.5+(j-1),1,10),linspace(-100,100,10),'k-');
end
ymax = max([-1.*min(Km_neg-2.*Ks_neg) max(Km_pos+2.*Ks_pos)]);
ymin = min([min(Km_neg-2.*Ks_neg) -1.*max(Km_pos+2.*Ks_pos)]);
if(ymax < 5 && ymin > -5)
    set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
        'xtick',1:5,'xticklabel','','YTick',-100:1:100);
    axis([0.5 length(Km_pos)+0.5 -5.5 5.5]);
elseif(ymax < 15 || ymin > -15)
    set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
        'xtick',1:5,'xticklabel','','YTick',-100:5:100);
    axis([0.5 length(Km_pos)+0.5 -17 17]);
else
    set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
        'xtick',1:5,'xticklabel','','YTick',-100:20:100);
    axis([0.5 length(Km_pos)+0.5 ymin-5 ymax+5]);
end
ylabel('Network Density (\%)','Interpreter','latex','FontSize',20);
a = gca; a.XAxis.FontSize = 14;
XTickString{1} = '\begin{tabular}{c}Static FCN\end{tabular}';
XTickString{2} = '\begin{tabular}{c}Raw\\mdFCN\end{tabular}';
XTickString{3} = '\begin{tabular}{c}Latent\\mdFCN\end{tabular}';
XTickString{4} = '\begin{tabular}{c}Decomposed\\mdFCN\end{tabular}';
XTickString{5} = '\begin{tabular}{c}Reconstructed\\mdFCN\end{tabular}';
set(gca,'xticklabel',XTickString,'TickLabelInterpreter','latex','XTickLabelRotationMode','manual');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[50 50 1400 700]); ax1 = axes();
histogram(1128.*K_pos_null(:),'Normalization','probability','FaceColor',[1,0.41,0.16]); hold on;
histogram(-1128.*K_neg_null(:),'Normalization','probability','FaceColor',[0.07,0.62,1]); grid on;
set(gca,'fontsize',26,'fontweight','bold','YScale','log','Xtick',-5:5,'Xticklabels',[5:-1:0 1:5]);
axis([-5 5 1e-6 1.5]);
xlabel('Number of Significant Connections','Interpreter','latex','FontSize',32);
ylabel({'Probability of Finding','Significant Connections'},'Interpreter','latex','FontSize',32);
annotation(gcf,'textbox',[0.120142857142857 0.65 0.3 0.188923076923077],...
    'VerticalAlignment','middle',...
    'String',{'Null Distribution','(Permutation)'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
title({'Network Density (\%)',''},'Interpreter','latex','FontSize',32);
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
set(ax2,'fontsize',26,'fontweight','bold');
ax2.YAxis.Visible = 'off';
ax2.XLim = [-5 5];
ax2.XTick = -5:5;
ax2.XTickLabel = [round(100.*(5:-1:1)./1128,2) 0 round(100.*(1:5)./1128,2)];
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Fig_6_cohort_mdFCN','-dpdf','-r400');
    print(2,'Results\Figures\Fig_6_latent_mdFCN','-dpdf','-r400');
    print(3,'Results\Figures\Fig_6_decomposed_mdFCN','-dpdf','-r400');
    print(4,'Results\Figures\Fig_6_reconstructed_mdFCN','-dpdf','-r400');
    print(5,'Results\Figures\Fig_6_cohort_mdFCN_zoomed','-dpdf','-r400');
    print(6,'Results\Figures\Fig_6_latent_mdFCN_zoomed','-dpdf','-r400');
    print(7,'Results\Figures\Fig_6_decomposed_mdFCN_zoomed','-dpdf','-r400');
    print(8,'Results\Figures\Fig_6_reconstructed_mdFCN_zoomed','-dpdf','-r400');
    print(9,'Results\Figures\Fig_6_verification','-dpdf','-r400');
    print(10,'Results\Figures\Fig_6_null_distribution','-dpdf','-r400');
end
