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
% This script generates Figs. 8-10 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'AED';    % HC or AED
fdr_thr = 0.025;    % 5% false discovery rate +/- threshold
alpha   = 0.05;     % Statistical significance level

%% Load density results
D_pos = []; D_neg = [];
for subset = 1:2
    load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'],'K_pos','K_neg');
    D_pos = cat(2,D_pos,K_pos);
    D_neg = cat(2,D_neg,K_neg);
end

%% Get correlation density
y_pos = []; y_neg = [];
mask_pos = false(length(D_pos),length(D_pos{1}));
mask_neg = false(length(D_neg),length(D_neg{1}));
for i = 1:length(D_pos)
    y_pos = cat(1,y_pos,D_pos{i}(:)');
    y_neg = cat(1,y_neg,D_neg{i}(:)');
    mask_pos(i,:) = abs(D_pos{i}) >= fdr_thr;
    mask_neg(i,:) = abs(D_neg{i}) >= fdr_thr;
end

%% Null distribution
R_pos = []; R_neg = [];
P_pos = []; P_neg = [];
for subset = 1:2
    load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_permutations']);
    R_pos = cat(2,R_pos,r_pos); P_pos = cat(2,P_pos,p_pos);
    R_neg = cat(2,R_neg,r_neg); P_neg = cat(2,P_neg,p_neg);
end
clear p_pos p_neg r_pos r_neg;
K_pos = []; K_neg = [];
for i = 1:size(R_pos,2)
    for k = 1:size(R_pos,1)
        [k_pos, ~] = corr_density(R_pos(k,i), P_pos(k,i), alpha);
        [~, k_neg] = corr_density(R_neg(k,i), P_neg(k,i), alpha);
        K_pos = [K_pos k_pos]; K_neg = [K_neg k_neg];
    end
end
clear R_pos R_neg P_pos P_neg k_pos k_neg;

%% Density verification results
Km_pos = zeros(7,5); Km_neg = zeros(7,5);
Ks_pos = zeros(7,5); Ks_neg = zeros(7,5);
for test_num = 1:5
    load(['Results\' grp_idx '\Verification\subset_1_test_' num2str(test_num)]);
    temp1 = Kp;  temp2 = Kn;
    load(['Results\' grp_idx '\Verification\subset_2_test_' num2str(test_num)]);
    KKp = 100.*[temp1; Kp]; KKn = -100.*[temp2; Kn];
    Km_pos(:,test_num) = mean(KKp,2);
    Km_neg(:,test_num) = mean(KKn,2);
    Ks_pos(:,test_num) = std(KKp,[],2);
    Ks_neg(:,test_num) = std(KKn,[],2);
end

%% Plotting
XTickString_a{1} = '\begin{tabular}{c}Neurological\\Score (C1)\end{tabular}';
XTickString_a{2} = '\begin{tabular}{c}Neurological\\Score (C2)\end{tabular}';
XTickString_a{3} = '\begin{tabular}{c}Cognition\end{tabular}';
XTickString_a{4} = '\begin{tabular}{c}Language\\Comprehension\end{tabular}';
XTickString_a{5} = '\begin{tabular}{c}Language\\Production\end{tabular}';
XTickString_a{6} = '\begin{tabular}{c}Fine\\Motor Skills\end{tabular}';
XTickString_a{7} = '\begin{tabular}{c}Gross\\Motor Skills\end{tabular}';
XTickString_b{1} = '\begin{tabular}{c}Static FCN\end{tabular}';
XTickString_b{2} = '\begin{tabular}{c}Raw\\mdFCN\end{tabular}';
XTickString_b{3} = '\begin{tabular}{c}Latent\\mdFCN\end{tabular}';
XTickString_b{4} = '\begin{tabular}{c}Decomposed\\mdFCN\end{tabular}';
XTickString_b{5} = '\begin{tabular}{c}Reconstructed\\mdFCN\end{tabular}';
ccc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; ...
    0.494 0.184 0.556; 0.929 0.694 0.125]; p1 = [];
legned_str = {'low-$\delta$','high-$\delta$','$\theta$','$\alpha$','$\beta$'};


figure('Color',[1,1,1],'Position',[25 25 1400 450]);
h1 = bar(100.*y_pos,1,'hist'); hold on;
h2 = bar(100.*y_neg,1,'hist'); grid on;
for i = 1:5, h1(i).FaceColor = ccc(i,:); h2(i).FaceColor = ccc(i,:); p1(i) = h1(i); end
h3 = fill([-10 -10 10 10], 100.*[-1*fdr_thr fdr_thr fdr_thr -1*fdr_thr], [0.75 0.75 0.75]);
set(h3,'FaceAlpha',0.5,'Linestyle','none');
for i = 1:size(y_pos,1)
    for j = 1:size(y_pos,2)
        plot(repmat(0.5+(i-1),1,10),linspace(-100,100,10),'k-');
    end
end
set(gca,'fontweight','bold','FontSize',22,'XGrid','off');
ylabel('Network Density (\%)','Interpreter','latex','FontSize',26);
legend(p1,legned_str,'Orientation','horizontal', ...
    'Location','southeast','interpreter','latex','fontsize',26);
a = gca; a.XAxis.FontSize = 18;
set(gca,'xtick',1:size(y_pos,1),'xticklabel',XTickString_a,'TickLabelInterpreter','latex', ...
    'XTickLabelRotationMode','manual');
axis([0.5 size(y_pos,1)+0.5 -89 89]);
set(gca,'TickLabelInterpreter','latex','YTick',-80:20:80);
title(grp_idx,'interpreter','latex','fontsize',32);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[50 50 1400 700]); ax1 = axes();
histogram(1128.*K_pos,'Normalization','probability','FaceColor',[1,0.41,0.16]); hold on;
histogram(-1128.*K_neg,'Normalization','probability','FaceColor',[0.07,0.62,1]); grid on;
set(gca,'fontsize',26,'fontweight','bold','YScale','log','Xtick',-5:5,'Xticklabels',[5:-1:0 1:5]);
legend('Positively Correlated','Negatively Correlated','Interpreter','latex','FontSize',32);
xlabel('Number of Significant Connections','Interpreter','latex','FontSize',32);
ylabel({'Probability of Finding','Significant Connections'},'Interpreter','latex','FontSize',32);
axis([-5 5 1e-6 1.5]);
annotation(gcf,'textbox',[0.120142857142857 0.55 0.3 0.188923076923077],...
    'VerticalAlignment','middle',...
    'String',grp_idx,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',50,...
    'FitBoxToText','off',...
    'EdgeColor','none','Interpreter','latex');
title({'Network Density (\%)',''},'Interpreter','latex','FontSize',32);
ax2 = axes('Position', get(ax1,'Position'),'XAxisLocation','top','Color','none','XColor','k');
set(ax2,'fontsize',26,'fontweight','bold');
ax2.YAxis.Visible = 'off'; ax2.XLim = [-5 5]; ax2.XTick = -5:5;
ax2.XTickLabel = [round(100.*(5:-1:1)./1128,2) 0 round(100.*(1:5)./1128,2)];
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

for i = 1:size(Km_pos,1)
    mean_pos = Km_pos(i,:); mean_neg = Km_neg(i,:);
    std_pos  = Ks_pos(i,:); std_neg  = Ks_neg(i,:);
    figure('Color',[1,1,1],'Position',[25 25 750 350]);
    stem(1:length(mean_pos),mean_pos,'filled','Color',[0.85,0.33,0.1], ...
        'LineStyle','-','LineWidth',12,'Marker','none'); hold on
    stem(1:length(mean_neg),mean_neg,'filled','Color',[0,0.45,0.74], ...
        'LineStyle','-','LineWidth',12,'Marker','none'); grid on;
    errorbar(1:length(mean_pos),mean_pos,zeros(1,length(mean_pos)),2.*std_pos, ...
        'LineStyle','none','Color',[0.85,0.33,0.1],'LineWidth',3,'CapSize',12,'Marker','none');
    errorbar(1:length(mean_neg),mean_neg,2.*std_neg,zeros(1,length(mean_neg)), ...
        'LineStyle','none','Color',[0,0.45,0.74],'LineWidth',3,'CapSize',12,'Marker','none');
    for j = 1:length(mean_pos)
        plot(repmat(0.5+(j-1),1,10),linspace(-100,100,10),'k-');
    end
    ymax = max([-1.*min(mean_neg-2.*std_neg) max(mean_pos+2.*std_pos)]);
    ymin = min([min(mean_neg-2.*std_neg) -1.*max(mean_pos+2.*std_pos)]);
    if(ymax < 5 && ymin > -5)
        set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
            'xtick',1:5,'xticklabel','','YTick',-100:1:100);
        axis([0.5 length(mean_pos)+0.5 -5.5 5.5]);
    elseif(ymax < 15 || ymin > -15)
        set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
            'xtick',1:5,'xticklabel','','YTick',-100:5:100);
        axis([0.5 length(mean_pos)+0.5 -17 17]);
    else
        set(gca,'fontweight','bold','FontSize',18,'XGrid','off', ...
            'xtick',1:5,'xticklabel','','YTick',-100:20:100);
        axis([0.5 length(mean_pos)+0.5 ymin-5 ymax+5]);
    end
    ylabel('Network Density (\%)','Interpreter','latex','FontSize',20);
    a = gca; a.XAxis.FontSize = 14;
    set(gca,'xticklabel',XTickString_b,'TickLabelInterpreter','latex','XTickLabelRotationMode','manual');
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,['Results\Figures\Supp_Fig_8_correlation_density_' grp_idx],'-dpdf','-r400');
    print(2,['Results\Figures\Supp_Fig_8_null_distribution_' grp_idx],'-dpdf','-r400');
    print(3,['Results\Figures\Supp_Fig_9_verification_1_' grp_idx],'-dpdf','-r400');
    print(4,['Results\Figures\Supp_Fig_9_verification_2_' grp_idx],'-dpdf','-r400');
    print(5,['Results\Figures\Supp_Fig_10_verification_1_' grp_idx],'-dpdf','-r400');
    print(6,['Results\Figures\Supp_Fig_10_verification_2_' grp_idx],'-dpdf','-r400');
    print(7,['Results\Figures\Supp_Fig_10_verification_3_' grp_idx],'-dpdf','-r400');
    print(8,['Results\Figures\Supp_Fig_10_verification_4_' grp_idx],'-dpdf','-r400');
    print(9,['Results\Figures\Supp_Fig_10_verification_5_' grp_idx],'-dpdf','-r400');
end