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
% This script generates the correlation analysis results in Fig. 4. Note that the
% neonates neurocognitive scores are not supplied and the code uses random
% numbers as example.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx   = 'HC';    % HC or AED
corr_dir  = '+';      % Correlation direction +/-
plot_idx  = 2:6;      % Plot results for these score indices
focus_idx = 3;        % Expand results for this score index (1-7)
alpha     = 0.05;     % Statistical significance level
fdr_thr   = 0.025;    % 5% false discovery rate +/- threshold

%% Load all requirements
load('Head Model\FidelityOperator.mat');
load('Head Model\Atlas.mat');
load('Head Model\surfaces\cortex_smoothed.mat');
load('Head Model\MyAtlas_n58.mat');
mask = logical(FidelityOperator);
mask = unfold_upper_matrix(mask);

%% Load results
D_pos = []; D_neg = [];
R_pos = []; R_neg = [];
P_pos = []; P_neg = [];
X_pos = []; X_neg = [];
for subset = 1:2
    load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'], ...
        'K_pos','K_neg','r_pos','r_neg','p_pos', 'p_neg','x_pos','x_neg');
    D_pos = cat(2,D_pos,K_pos); D_neg = cat(2,D_neg,K_neg);
    R_pos = cat(2,R_pos,r_pos); R_neg = cat(2,R_neg,r_neg);
    P_pos = cat(2,P_pos,p_pos); P_neg = cat(2,P_neg,p_neg);
    X_pos = cat(2,X_pos,x_pos); X_neg = cat(2,X_neg,x_neg);
end

%% Get correlation density
y_pos = []; y_neg = [];
for i = 1:length(D_pos)
    y_pos = cat(1,y_pos,D_pos{i}(:)');
    y_neg = cat(1,y_neg,D_neg{i}(:)');
end
y_pos = y_pos(plot_idx,:);
y_neg = y_neg(plot_idx,:);

%% Get averaged correlation
if(strcmp(corr_dir,'+'))
    mask_band = abs(D_pos{focus_idx}) >= fdr_thr;
    R = R_pos{focus_idx};
    P = P_pos{focus_idx};
else
    mask_band = abs(D_neg{focus_idx}) >= fdr_thr;
    R = R_neg{focus_idx};
    P = P_neg{focus_idx};
end
corr_m = mean(R,2,'omitnan');  corr_m(~mask_band) = nan;
corr_s = std(R,0,2,'omitnan'); corr_s(~mask_band) = nan;
corr_q = mean(P,2,'omitnan');  corr_q(~mask_band) = nan;
corr_d = std(P,0,2,'omitnan'); corr_d(~mask_band) = nan;

%% Get brain network
if(strcmp(corr_dir,'+'))
    X = X_pos{focus_idx};
else
    X = X_neg{focus_idx};
end
XX = squeeze(mean(X,2,'omitnan'));               % average the subjects connections
XX = (XX - min(XX(:)))./(max(XX(:))-min(XX(:))); % Standardize to be from 0 to 1
XX(isnan(XX)) = 0;                               % remove connections with insignificant correlations
XX = XX.^2;                                      % Improve visualization by squaring

%% Get linear predictive models
xx = mean(X(mask_band,:,:),[1 3],'omitnan');     % average connectivity across all bands and connections


%%%%%%%%%%%%%%%%%%%%%%% Original Code %%%%%%%%%%%%%%%%%%%%%%%%%
% load(['Results\' grp_idx '\scores.mat'],'Scores','idx_all'); %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Alternative Code %%%%%%%%%%%%%%%%%%%%%%
load(['Results\' grp_idx '\scores.mat'],'idx_all');          %%
Scores = [4*rand(size(idx_all,1),2)-2, ...                   %%              
    randi([5 20],size(idx_all,1),5)];                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(focus_idx > 2)
    score = Scores(idx_all(:,2),focus_idx);
else
    score = Scores(idx_all(:,1),focus_idx);
end
if(any(mask_band))
    [PP, SS] = polyfit(xx,score,1);
    [yfit, delta] = polyval(PP,xx,SS);
end

%% Plotting
str = {'Neurological Score (C1)','Neurological Score (C2)','Cognition',...
    'Language Comprehension','Language Production','Fine Motor Skills','Gross Motor Skills'};
ccc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; ...
    0.494 0.184 0.556; 0.929 0.694 0.125]; p1 = [];
legned_str = {'low-$\delta$','high-$\delta$','$\theta$','$\alpha$','$\beta$'};

figure('Color',[1,1,1],'Position',[25 25 750 450]); 
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
    'Location','southwest','interpreter','latex','fontsize',26);

a = gca; a.XAxis.FontSize = 14;
XTickString{1} = '\begin{tabular}{c}Neurological\\Score (C2)\end{tabular}';
XTickString{2} = '\begin{tabular}{c}Cognition\end{tabular}';
XTickString{3} = '\begin{tabular}{c}Language\\Comprehension\end{tabular}';
XTickString{4} = '\begin{tabular}{c}Language\\Production\end{tabular}';
XTickString{5} = '\begin{tabular}{c}Fine\\Motor Skills\end{tabular}';
set(gca,'xtick',1:size(y_pos,1),'xticklabel',XTickString,'TickLabelInterpreter','latex', ...
    'XTickLabelRotationMode','manual');
if(strcmp(grp_idx,'HC'))
    axis([0.5 size(y_pos,1)+0.5 -47 47]);
    set(gca,'TickLabelInterpreter','latex','YTick',-50:10:50);
else
    axis([0.5 size(y_pos,1)+0.5 -89 89]);
    set(gca,'TickLabelInterpreter','latex','YTick',-80:20:80);
end
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

N = 101; dx = 1.05; t = linspace(0, 1, N);
figure('Color',[1,1,1],'Position',[20 100 1400 300]); hold on;
for f = 1:size(XX,1)
    axis([-0.5 0.5+(size(XX,1)-1)*dx -0.5 0.75]);
    if(mask_band(f))
        [~, y] = fold_connectivity(XX(f,:), mask, 0);
        for i = 1:length(MyAtlas.Parcels)
            for j = i+1:length(MyAtlas.Parcels)
                if y(i, j) ~= 0 % skip zeros (no connections)
                    pts = kron((1-t).^2, MyAtlas.Circular_xy(i, :)') + ...
                        kron(2*(1-t).^t, [0; 0]) + kron(t.^2, MyAtlas.Circular_xy(j, :)');
                    plot3(pts(1,:)+(f-1)*dx, pts(2,:), repelem(0,1,N), 'Color', ...
                        [ccc(f,:) y(i,j)], 'Linewidth', 1);
                end
            end
        end

        % Nodes
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2), 0, 40, ...
            'MarkerFaceColor', [0.00 0.00 0.00], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2), 0, 40, ...
            'MarkerFaceColor', [1.00 0.00 1.00], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2), 0, 40, ...
            'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2), 0, 40, ...
            'MarkerFaceColor', [0.20 0.60 0.00], 'MarkerEdgeColor', 'none');
        % titles
        pos = [0.14+(0.156)*(f-1) 0.798 0.128714285714286 0.106666666666667];
        annotation(gcf,'textbox',pos,'String',['$r=' num2str(round(corr_m(f),3)) '$\,\,\,$p=' ...
            num2str(round(corr_q(f),3)) '$'],'Interpreter','latex','HorizontalAlignment', ...
            'center','FontWeight','bold','FontSize',16,'FitBoxToText','off','EdgeColor','none', ...
            'VerticalAlignment','middle','Color',ccc(f,:));
    else
        annotation(gcf,'textbox',...
            [0.12+(0.156)*(f-1) 0.225 0.45*(2/5) 0.425],...
            'VerticalAlignment','middle',...
            'String','n.s.',...
            'Interpreter','latex',...
            'HorizontalAlignment','center',...
            'FontSize',40,...
            'FitBoxToText','off',...
            'EdgeColor','none');
        % Nodes
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2), 0, 40, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2), 0, 40, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2), 0, 40, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
        scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1)+(f-1)*dx, ...
            MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2), 0, 40, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
    end
end
axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[25 25 850 500]); hold on;
if(any(mask_band))
    [~,II] = sort(yfit + delta);
    scatter(xx,score,256,'Marker','o','MarkerFaceColor',[0.5,0.5,0.5], ...
        'MarkerEdgeColor','none','MarkerFacealpha',0.5);
    plot(xx,yfit,'-','LineWidth',6,'Color','k');
    plot(xx(II),yfit(II)+2.*delta(II),'-.','Color','k','linewidth',4);
    plot(xx(II),yfit(II)-2.*delta(II),'-.','Color','k','linewidth',4);
    axis([min(xx)-0.1*(max(xx)-min(xx)) max(xx)+0.1*(max(xx)-min(xx)) min(yfit(II)-2.*delta(II))-2 max(yfit(II)+2.*delta(II))+2]);

    if(strcmp(corr_dir,'+'))
        pos = [0.662588235294117 0.2 0.22041176470588 0.15];
    else
        pos = [0.658294117647057 0.775200000000001 0.22041176470588 0.1376];
    end    
    str1 = ['\begin{tabular}{c}$r=' num2str(round(mean(corr_m,'omitnan'),3))...
        '$\\$p=' num2str(round(mean(corr_q,'omitnan'),3)) '$\end{tabular}'];
    annotation(gcf,'textbox',pos,'String',str1,'Interpreter','latex','HorizontalAlignment', ...
        'center','FontWeight','bold','FontSize',24,'FitBoxToText','off','EdgeColor','k', ...
        'BackgroundColor','w','VerticalAlignment','middle');
else
    axis off;
    annotation(gcf,'textbox',...
        [0.0274 0.0492 0.9726 0.9321],...
        'VerticalAlignment','middle',...
        'String','n.s.',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',50,...
        'FitBoxToText','off',...
        'EdgeColor','k');
end
set(gca,'FontSize',24,'FontWeight','bold','YAxisLocation','right'); box on; grid on;
xlabel('Standardized Connectivity','FontSize',30,'FontWeight','bold','Interpreter','latex');
ylabel([str{focus_idx} ' Score'],'FontSize',30,'FontWeight','bold','Interpreter','latex');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,['Results\Figures\Fig_4_Correlation_density_' grp_idx],'-dpdf','-r400');
    print(2,['Results\Figures\Fig_4_brain_network_' grp_idx],'-dpdf','-r400');
    print(3,['Results\Figures\Fig_4_linear_regression_' grp_idx],'-dpdf','-r400');
end