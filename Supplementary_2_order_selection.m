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
% This script generates plots of Fig. 2 in the supplementary material.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC'; % HC or AED
subset  = 1;    % 1 or 2
N_rep   = 10;   % NMF repetitions
Pmax    = 40;   % Maximum NMF order
Nb      = 5;    % Number of bands
per     = 0.05; % Order selection percentage

%% Model automatic order selection
Entropy = zeros(Nb, Pmax, N_rep);
Error   = zeros(Nb, Pmax, N_rep);
for i = 1:Nb
    for p = 1:Pmax
        fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
            num2str(i) '_order_' num2str(p)];
        load(fname,'trial_entropy','trial_error','Ent');
        Entropy(i,p,:) = trial_entropy./Ent(i);
        Error(i,p,:)   = trial_error;
    end
end
out_error   = order_selection(1-Error, per, 1);
out_entropy = order_selection(Entropy, per, 1);

%% Plotting NMF order selection
ccc = [0 0.447 0.741; 0.85 0.325 0.098];
figure('Color',[1,1,1],'Position',[25 50 750 500]);
colororder({'k','k'});
yyaxis left;
ax1 = plot(out_entropy.mfit(1,:),'color',[0 0.447 0.741],'LineStyle','-',...
        'LineWidth',2,'Marker','none'); hold on;
ax2 = plot(out_entropy.order(1),out_entropy.mfit(1,out_entropy.order(1)), ...
    'o','MarkerFaceColor',[0 0.447 0.741],'MarkerEdgeColor','k','MarkerSize',15);
grid on; axis([0.75 Pmax 0.7 1]); set(gca,'fontweight','bold','FontSize',18);
ylabel('Normalized Estimation Entropy','Interpreter','latex','FontSize',22);
xlabel('NMF Model Order $P$','Interpreter','latex','FontSize',22);
yyaxis right
ax3 = plot(2:Pmax,out_entropy.mdif(1,:),'color',[0 0.447 0.741], ...
        'LineStyle','-.','LineWidth',2,'Marker','none');
ax4 = plot(linspace(-1,41,100),repelem(out_entropy.thrs(1),100), ...
    '--k','LineWidth',2);
ylabel('Entropy Rate of Change','Interpreter','latex','FontSize',22);
axis([0.75 Pmax 0 0.04]);
legend([ax1(1) ax3(1) ax4(1) ax2(1)],'Averaged Normalized Entropy', ...
    'Averaged Rate of Change','Threshold $\epsilon$','Selected NMF order', ...
    'Location','northwest','Interpreter','latex','Orientation','vertical', ...
    'NumColumns',2,'fontsize',14);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[25 50 750 500]);
colororder({'k','k'});
yyaxis left;
ax1 = plot(out_error.mfit(1,:),'color',[0.85 0.325 0.098],'LineStyle','-',...
        'LineWidth',2,'Marker','none'); hold on;
grid on; axis([0.75 Pmax 0.42 0.51]); set(gca,'fontweight','bold','FontSize',18);
ylabel('1 - Estimation Error','Interpreter','latex','FontSize',22);
xlabel('NMF Model Order $P$','Interpreter','latex','FontSize',22);
yyaxis right
ax2 = plot(2:Pmax,out_error.mdif(1,:),'color',[0.85 0.325 0.098], ...
        'LineStyle','-.','LineWidth',2,'Marker','none');
ylabel('(1 - Error) Rate of Change','Interpreter','latex','FontSize',22);
axis([0.75 Pmax 0 0.04]);
legend([ax1(1) ax2(1)],'Averaged 1-Error', ...
    'Averaged Rate of Change','Location','northwest', ...
    'Interpreter','latex','Orientation','vertical','NumColumns',2,'fontsize',14)
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));


%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Supp_Fig_2_entropy','-dpdf','-r400');
    print(1,'Results\Figures\Supp_Fig_2_error','-dpdf','-r400');
end