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
% This script generates plots of Fig. 3 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));
addpath(genpath('Head Model'));

%% Load data
load('MyAtlas_n58.mat');
load('cortex_smoothed.mat');

%% Main
K = 8;
Mk = cell(1,5);
Mk{1} = zeros(58,58);
Mk{2} = zeros(58,58);
Mk{3} = zeros(58,58);
Mk{4} = zeros(58,58);
Mk{5} = zeros(58,58);
Mk{1}(1:K,1*K+2:2*K+1) = rand(K).^2;
Mk{2}(1:K,3*K+2:4*K+1) = rand(K).^2;
Mk{3}(1:K,5*K+2:6*K+1) = rand(K).^2;
Mk{4}(2*K+1:3*K,6*K+1:7*K) = rand(K).^2;
Mk{5}(4*K+1:5*K,6*K+1:7*K) = rand(K).^2;
for i = 1:5
    Mk{i} = Mk{i} + Mk{i}';
    Mk{i} = Mk{i}./max(Mk{i}(:));
end

%% Plotting
N = 101; dx = 1.05; ax = cell(1,5);
clr = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; ...
    0.494 0.184 0.556; 0.929 0.694 0.125];
for k = 1:5
    figure('Color',[1 1 1],'Position',[10 10 550 450]);
    imagesc(Mk{k}); colormap(gray);
    set(gca,'Xticklabels','','Yticklabels','');
    xlabel('Connections','fontsize',22,'FontWeight','bold');
    ylabel('Connections','fontsize',22,'FontWeight','bold');
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

figure('Color',[1 1 1],'Position',[10 10 550 450]); ax0 = axes;
plot3d_atlas(MyAtlas, flat_cx,60); hold(ax0,'on'); axis off; view([-90 90]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
for k = 1:5
    figure('Color',[1 1 1],'Position',[10 10 550 450]); ax{k} = axes;
    plot3d_atlas(MyAtlas, flat_cx,60); hold(ax{k},'on'); axis off; view([-90 90]);
    for i = 1:58
        for j = 1:58
            plot3(ax0,[MyAtlas.Centroids(i,1) MyAtlas.Centroids(j,1)], ...
                [MyAtlas.Centroids(i,2) MyAtlas.Centroids(j,2)], ...
                [MyAtlas.Centroids(i,3) MyAtlas.Centroids(j,3)], ...
                'Color', [clr(k,:) Mk{k}(i,j)], 'Linewidth',1.25);
            plot3(ax{k},[MyAtlas.Centroids(i,1) MyAtlas.Centroids(j,1)], ...
                [MyAtlas.Centroids(i,2) MyAtlas.Centroids(j,2)], ...
                [MyAtlas.Centroids(i,3) MyAtlas.Centroids(j,3)], ...
                'Color', [clr(k,:) Mk{k}(i,j)], 'Linewidth',1.25);
        end
    end
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

figure('Color',[1,1,1],'Position',[20 20 450 450]); hold on;
for k = 1:5
    t = linspace(0, 1, N);
    for i = 1:58
        for j = i+1:58
            if Mk{k}(i,j) ~= 0 % skip zeros (no connections)
                pts = kron((1-t).^2, MyAtlas.Circular_xy(i, :)') + ...
                    kron(2*(1-t).^t, [0; 0]) + kron(t.^2, MyAtlas.Circular_xy(j, :)');
                plot3(pts(1,:), pts(2,:), repelem(0,1,N), ...
                    'Color', [clr(k,:) Mk{k}(i,j)], 'Linewidth', 1);
            end
        end
    end
    % Nodes
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1), ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2), 0, 60, ...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1), ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2), 0, 60, ...
        'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1), ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2), 0, 60, ...
        'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1), ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2), 0, 60, ...
        'MarkerFaceColor', [0.2 0.6 0], 'MarkerEdgeColor', 'none');
end
xlim([-0.5 0.5]); ylim([-0.5 0.5]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[20 20 1400 250]); hold on;
for k = 1:5
    t = linspace(0, 1, N);
    for i = 1:58
        for j = i+1:58
            if Mk{k}(i,j) ~= 0 % skip zeros (no connections)
                pts = kron((1-t).^2, MyAtlas.Circular_xy(i, :)') + ...
                    kron(2*(1-t).^t, [0; 0]) + kron(t.^2, MyAtlas.Circular_xy(j, :)');
                plot3(pts(1,:)+(k-1)*dx, pts(2,:), repelem(0,1,N), ...
                    'Color', [clr(k,:) Mk{k}(i,j)], 'Linewidth', 1);
            end
        end
    end
    % Nodes
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 1)+(k-1)*dx, ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'O', 2), 0, 40, ...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 1)+(k-1)*dx, ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'C', 2), 0, 40, ...
        'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 1)+(k-1)*dx, ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'F', 2), 0, 40, ...
        'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', 'none');
    scatter3(MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 1)+(k-1)*dx, ...
        MyAtlas.Circular_xy(cell2mat(MyAtlas.Areas) == 'T', 2), 0, 40, ...
        'MarkerFaceColor', [0.2 0.6 0], 'MarkerEdgeColor', 'none');
end
xlim([-0.5 0.5+4*dx]); ylim([-0.5 0.5]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Supp_Fig_3_tensor_1','-dpdf','-r400');
    print(2,'Results\Figures\Supp_Fig_3_tensor_2','-dpdf','-r400');
    print(3,'Results\Figures\Supp_Fig_3_tensor_3','-dpdf','-r400');
    print(4,'Results\Figures\Supp_Fig_3_tensor_4','-dpdf','-r400');
    print(5,'Results\Figures\Supp_Fig_3_tensor_5','-dpdf','-r400');
    print(6,'Results\Figures\Supp_Fig_3_brain_network','-dpdf','-r400');
    print(7,'Results\Figures\Supp_Fig_3_brain_network_1','-dpdf','-r400');
    print(8,'Results\Figures\Supp_Fig_3_brain_network_2','-dpdf','-r400');
    print(9,'Results\Figures\Supp_Fig_3_brain_network_3','-dpdf','-r400');
    print(10,'Results\Figures\Supp_Fig_3_brain_network_4','-dpdf','-r400');
    print(11,'Results\Figures\Supp_Fig_3_brain_network_5','-dpdf','-r400');
    print(12,'Results\Figures\Supp_Fig_3_multiplex','-dpdf','-r400');
    print(13,'Results\Figures\Supp_Fig_3_multiplex_detailed','-dpdf','-r400');
end