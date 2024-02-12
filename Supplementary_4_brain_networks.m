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
% This script generates Figs. 4-7 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
corr_dir  = '+';      % Correlation direction +/-
fdr_thr   = 0.025;    % 5% false discovery rate +/- threshold

%% Load all requirements
load('Head Model\FidelityOperator.mat');
load('Head Model\Atlas.mat');
load('Head Model\surfaces\cortex_smoothed.mat');
load('Head Model\MyAtlas_n58.mat');
mask = logical(FidelityOperator);
mask = unfold_upper_matrix(mask);

%% Load results
X_hc = []; D_hc = []; X_ad = []; D_ad = [];
for subset = 1:2
    if(strcmp(corr_dir,'+'))
        load(['Results\HC\Selection\subset_' num2str(subset) '_selected'],'x_pos','K_pos');
        X_hc = cat(2,X_hc,x_pos); D_hc = cat(2,D_hc,K_pos);
        load(['Results\AED\Selection\subset_' num2str(subset) '_selected'],'x_pos','K_pos');
        X_ad = cat(2,X_ad,x_pos); D_ad = cat(2,D_ad,K_pos);
    else
        load(['Results\HC\Selection\subset_' num2str(subset) '_selected'],'x_neg','K_neg');
        X_hc = cat(2,X_hc,x_neg); D_hc = cat(2,D_hc,K_neg);
        load(['Results\AED\Selection\subset_' num2str(subset) '_selected'],'x_neg','K_neg');
        X_ad = cat(2,X_ad,x_neg); D_ad = cat(2,D_ad,K_neg);
    end
end

%% Get brain network
mask_hc = cell(1,length(X_hc));
mask_ad = cell(1,length(X_hc));
for i = 1:length(X_hc)
    mask_hc{i} = abs(D_hc{i}) >= fdr_thr;
    mask_ad{i} = abs(D_ad{i}) >= fdr_thr;

    XX = squeeze(mean(X_hc{i},2,'omitnan'));         % average the subjects connections
    XX = (XX - min(XX(:)))./(max(XX(:))-min(XX(:))); % Standardize to be from 0 to 1
    XX(isnan(XX)) = 0;                               % remove connections with insignificant correlations
    X_hc{i} = XX.^2;                                 % Improve visualization by squaring

    XX = squeeze(mean(X_ad{i},2,'omitnan'));         % average the subjects connections
    XX = (XX - min(XX(:)))./(max(XX(:))-min(XX(:))); % Standardize to be from 0 to 1
    XX(isnan(XX)) = 0;                               % remove connections with insignificant correlations
    X_ad{i} = XX.^2;                                 % Improve visualization by squaring
end

%% Plotting
ccc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.929 0.694 0.125];
N = 101; dx = 1.05; t = linspace(0, 1, N);
Band_str = {'low-$\delta$','high-$\delta$','$\theta$','$\alpha$','$\beta$'};
Score_str{1} = 'Neurological Score (C1)';
Score_str{2} = 'Neurological Score (C2)';
Score_str{3} = 'Cognition';
Score_str{4} = 'Language Comprehension';
Score_str{5} = 'Language Production';
Score_str{6} = 'Fine Motor Skills';
Score_str{7} = 'Gross Motor Skills';
for n = 1:length(X_hc)
    figure('Color',[1,1,1],'Position',[20 100 1400 500]); hold on;
    for f = 1:size(X_hc{n},1)
        if(mask_hc{n}(f))
            [~, y] = fold_connectivity(X_hc{n}(f,:), mask, 0);
            for i = 1:length(MyAtlas.Parcels)
                for j = i+1:length(MyAtlas.Parcels)
                    if y(i, j) ~= 0 % skip zeros (no connections)
                        pts = kron((1-t).^2, MyAtlas.Circular_xy(i, :)') + ...
                            kron(2*(1-t).^t, [0; 0]) + kron(t.^2, MyAtlas.Circular_xy(j, :)');
                        plot3(pts(1,:)+(f-1)*dx, pts(2,:)+1.1, repelem(0,1,N), ...
                            'Color', [ccc(f,:) y(i,j)], 'Linewidth', 1);
                    end
                end
            end
            % Nodes
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.00 0.00 0.00], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [1.00 0.00 1.00], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.20 0.60 0.00], 'MarkerEdgeColor', 'none');
        else
            annotation(gcf,'textbox',...
                [0.12+(0.156)*(f-1) 0.525 0.45*(2/5) 0.425],...
                'VerticalAlignment','middle',...
                'String','n.s.',...
                'Interpreter','latex',...
                'HorizontalAlignment','center',...
                'FontSize',40,...
                'FitBoxToText','off',...
                'EdgeColor','none');
            % Nodes
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
            scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1)+(f-1)*dx, ...
                MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2)+1.075, 0, 40, ...
                'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
        end

        if(mask_ad{n}(f))
            [~, y] = fold_connectivity(X_ad{n}(f,:), mask, 0);
            for i = 1:length(MyAtlas.Parcels)
                for j = i+1:length(MyAtlas.Parcels)
                    if y(i, j) ~= 0 % skip zeros (no connections)
                        pts = kron((1-t).^2, MyAtlas.Circular_xy(i, :)') + ...
                            kron(2*(1-t).^t, [0; 0]) + kron(t.^2, MyAtlas.Circular_xy(j, :)');
                        plot3(pts(1,:)+(f-1)*dx, pts(2,:), repelem(0,1,N), ...
                            'Color', [ccc(f,:) y(i,j)], 'Linewidth', 1);
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
        else
            annotation(gcf,'textbox',...
                [0.12+(0.156)*(f-1) 0.1 0.45*(2/5) 0.425],...
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
        if(n == 1 || n == 3)
            annotation(gcf,'textbox',[0.12+(0.156)*(f-1) 0.75 0.45*(2/5) 0.425],'VerticalAlignment', ...
                'middle','String',Band_str{f},'Interpreter','latex','HorizontalAlignment', ...
                'center','FontSize',24,'FitBoxToText','off','EdgeColor','none');
        elseif(n == 2 || n == 7)
            annotation(gcf,'textbox',[0.12+(0.156)*(f-1) 0 0.45*(2/5) 0.1],'VerticalAlignment', ...
                'middle','String',Band_str{f},'Interpreter','latex','HorizontalAlignment', ...
                'center','FontSize',24,'FitBoxToText','off','EdgeColor','none');
        end
    end
    annotation(gcf,'textbox',[0.1 0.155 0.25 0.2],'VerticalAlignment', ...
        'middle','String',Score_str{n},'Interpreter','latex','HorizontalAlignment', ...
        'center','FontSize',38,'FitBoxToText','off','EdgeColor','none','rotation',90);
    annotation(gcf,'textbox',[0.77 0.625 0.35 0.2],'VerticalAlignment', ...
        'middle','String','HC','Interpreter','latex','HorizontalAlignment', ...
        'center','FontSize',24,'FitBoxToText','off','EdgeColor','none');
    annotation(gcf,'textbox',[0.77 0.2 0.35 0.2],'VerticalAlignment', ...
        'middle','String','AED','Interpreter','latex','HorizontalAlignment', ...
        'center','FontSize',24,'FitBoxToText','off','EdgeColor','none');
    xlim([-0.5 0.5+(size(X_hc{1},1)-1)*dx]); ylim([-0.5 1.6]); axis off;
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    for n = 1:length(X_hc)
        if(strcmp(corr_dir,'+'))
            if(n <= 2)
                print(n,['Results\Figures\Supp_Fig_4_brain_network_' num2str(n)],'-dpdf','-r400');
            else
                print(n,['Results\Figures\Supp_Fig_6_brain_network_' num2str(n-2)],'-dpdf','-r400');
            end
        else
            if(n <= 2)
                print(n,['Results\Figures\Supp_Fig_5_brain_network_' num2str(n)],'-dpdf','-r400');
            else
                print(n,['Results\Figures\Supp_Fig_7_brain_network_' num2str(n-2)],'-dpdf','-r400');
            end
        end
    end
end