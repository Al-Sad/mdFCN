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
% This script generates the study overview in Fig. 1. Note that the
% neonates raw EEG and neurocognitive scores are not supplied and the code
% uses random numbers as example.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));
addpath(genpath('Head Model'));

%% Parameters
grp_idx   = 'AED'; % HC or AED
subset    = 2;    % 1 or 2 for short- or long-term scores
sub_idx   = 2;
focus_idx = 7;
fs        = 100;
T         = 3*60;
fdr_thr   = 0.025;
ccc = [0 0.447 0.741
    0.85 0.325 0.098
    0.466 0.674 0.188
    0.494 0.184 0.556
    0.929 0.694 0.125];

%% EEG recording
%%%%%%%%%%%%%%%%%%%%%%% Original Code %%%%%%%%%%%%%%%%%%%%%%%%%
% data_folder = ['Dataset\' grp_idx '\'];                      %%                          
% load([data_folder '\Subj_' sprintf('%03d',sub_idx) ...       %% 
%     '\TA_eeg.mat'],'data');                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Alternative Code %%%%%%%%%%%%%%%%%%%%%%%%%
data.Fs = 100;                                                  %%
data.eeg = movmean(randn(19,18000),32,2);                       %%
data.Labels = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3',...    %%
    'C3','Cz','C4','T4','T5','P3','Pz','P4','T6','O1','O2'}';   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = data.eeg';
xx = [x(:,[1 2]) nan(size(x,1),2) x(:,end)];
figure('Color',[1 1 1],'Position',[100 100 750 450]);
plot_signals(xx, data.Labels, fs, 'eeg'); box on;
set(gca,'FontWeight','bold','yticklabels','','xticklabels','', ...
    'FontName','times','fontsize',26);
xlabel('Time (s)','fontsize',28);
ylabel('EEG channels','fontsize',28); xlim([0 T]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% EEG pre-processing
flt_eeg = iir_filters(fs);
L = size(x, 1);
% Flipped signal
signal_ud = flipud(x);
% Add pieces to signal (to fix edge effects)
signal_big = [signal_ud; x; signal_ud];
% Number of bandpass filters (= fr.bands); columns = bandpass filters
B = size(flt_eeg, 2);
% Init cell array for bandpass filtered signals
y{1,B} = [];
yy{1,B} = [];
% Filter signals (bandpass filter = HPF + LPF)
for fr = 1:B
    buf = filtfilt(flt_eeg{1, fr}, signal_big); % HPF/cutoff = 0.85xFc
    buf = filtfilt(flt_eeg{2, fr}, buf);        % LPF/cutoff = 1.15xFc
    y{1,fr} = buf(L+1:L+L,:);    % cut signal_big >> orig.sig.
    yy{1,fr} = [y{1,fr}(:,[1 2]) nan(size(x,1),3) y{1,fr}(:,end)];
end
figure('Color',[1 1 1],'Position',[100 100 750 550]);
for i = 1:B
    [h1,f1] = freqz(flt_eeg{1,i},2^10,fs);
    [h2,f2] = freqz(flt_eeg{2,i},2^10,fs);
    hold on; grid on;
    plot(f1,abs(h1).*abs(h2),'LineWidth',5,'Color',ccc(i,:));
end
set(gca,'FontWeight','bold','FontName','times','fontsize',26,'yticklabels','');
axis([0 25 0 1.2]); xlabel('Frequency (Hz)','fontsize',28);
legend('low-$\delta$','high-$\delta$','$\theta$','$\alpha$','$\beta$', ...
    'interpreter','latex','fontsize',48,'NumColumns',5,...
    'orientation','horizontal'); box on;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100 100 950 650]);
plot_mult_signals(yy, fs, data.Labels, 'eeg');
view([22 25]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Cortical signals reconstruction
load('Atlas.mat','Atlas');
load('InverseOperator.mat','InverseOperator');
load('CollapseOperator.mat','CollapseOperator');
Np = length(Atlas.Parcels); % number of parcels
z{1,B} = [];
zz{1,B} = [];
CollapseOperator = repmat(CollapseOperator, 1, L);
for k = 1:B
    % source signals
    src_buf = (InverseOperator * y{1,k}').* CollapseOperator;
    % parcel/cortical signals
    parcel_buf = zeros(Np, L);
    for j = 1:Np
        parcel_buf(j, :) = mean(src_buf(Atlas.Parcels{j, 1}, :));
    end
    z{1,k} = parcel_buf';
    zz{1,k} = [z{1,k}(:,[1 2]) nan(size(x,1),3) z{1,k}(:,end)];
end
figure('Color',[1 1 1],'Position',[100 100 950 650]);
plot_mult_signals(zz, fs, Atlas.Areas, 'eeg');
view([22 25]); axis off;
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Group level mdFCN tensor
load('Head Model\FidelityOperator.mat');
load(['Results\' grp_idx '\scores.mat'],'idx_all');
mask = logical(FidelityOperator);
mask = unfold_upper_matrix(mask);
load(['Results\' grp_idx '\FC\Subj_' num2str(sub_idx)],'dwPLI_TA');
T = dwPLI_TA(:,:,mask);
T = permute(T,[2 3 1]);
[N,M,P] = size(T);
figure('Color',[1 1 1],'Position',[100 100 600 500]); hold on;
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

data_folder = ['Results\' grp_idx '\FC'];
m = 1:length(dir([data_folder '*\*.mat']));
m = m(idx_all(:,subset));
T = [];
for i = 1:length(m)
    load([data_folder '\Subj_' num2str(m(i)) '.mat']);
    if(~isempty(dwPLI_TA))
        T = cat(2,T,dwPLI_TA(:,:,mask));
    end
    if(~isempty(dwPLI_AS))
        T = cat(2,T,dwPLI_AS(:,:,mask));
    end
end
T = permute(T,[2 3 1]);
[N,M,P] = size(T);
figure('Color',[1 1 1],'Position',[100 100 1000 250]); hold on;
for i = 1:P
    h = surface(i*ones(1,M),linspace(1-200*i,N-200*i,N),...
        repmat(linspace(M+100*i,1+100*i,M),N,1),repmat(M:-1:1,N,1));
    h.CData = T(:,:,i); caxis([0 1]);
    h.EdgeColor = 'none';
    h.FaceColor = 'interp';
end
axis off; view([-90 0]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Latent network extraction
T = [];
load(['Results\' grp_idx '\NMF\subset_' num2str(subset) '_order']);
order = out_entropy.order;
for i = 1:5
    fname = ['Results\' grp_idx '\NMF\subset_' num2str(subset) '_f'...
        num2str(i) '_order_' num2str(order(i))];
    load(fname,'W','H','g');
    x = W*H;
    T = cat(3,T,x);
end
[N,M,P] = size(T);
figure('Color',[1 1 1],'Position',[100 100 1000 250]); hold on;
for i = 1:P
    h = surface(i*ones(1,M),linspace(1-200*i,N-200*i,N),...
        repmat(linspace(M+100*i,1+100*i,M),N,1),repmat(M:-1:1,N,1));
    h.CData = T(:,:,i); caxis([0 1]);
    h.EdgeColor = 'none';
    h.FaceColor = 'interp';
end
axis off; view([-90 0]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Latent network decomposition
load(['Results\' grp_idx '\CPD\subset_' num2str(subset) '_block_averaged_network'])
for j = [1 5 6 7 9 13]
    T = permute(block_avg_tensor{j},[2 3 1]);
    [N,M,P] = size(T);
    figure('Color',[1 1 1],'Position',[100 100 250 400]); hold on;
    for i = 1:P
        h = surface(i*ones(1,M),linspace(1-5*i,N-5*i,N),...
            repmat(linspace(M+50*i,1+50*i,M),N,1),repmat(M:-1:1,N,1));
        h.CData = T(:,:,i); caxis([0 0.2]);
        h.EdgeColor = 'none';
        h.FaceColor = 'interp';
    end
    axis off; view([-90 0]);
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Component selection
if(focus_idx > 2)
    focus = focus_idx - 2;
else
    focus = focus_idx;
end
load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'],'I_pos');
load(['Results\' grp_idx '\CPD\subset_' num2str(subset) '_block_averaged_network']);
I = zeros(size(block_avg_tensor{1},1),length(block_avg_tensor));
for i = 1:size(block_avg_tensor{1},1)
    if(~isnan(I_pos(i,focus)))
        I(i,:) = fliplr(int2bit(I_pos(i,focus),length(block_avg_tensor))');
    end
end
for j = [1 5 6 7 9 13]
    T = permute(block_avg_tensor{j},[2 3 1]);
    [N,M,P] = size(T);
    figure('Color',[1 1 1],'Position',[100 100 250 400]); hold on;
    for i = 1:P
        h = surface(i*ones(1,M),linspace(1-5*i,N-5*i,N),...
            repmat(linspace(M+50*i,1+50*i,M),N,1),repmat(M:-1:1,N,1));
        h.CData = T(:,:,i); caxis([0 0.2]);
        h.EdgeColor = 'none';
        h.FaceColor = 'interp';
        h.FaceAlpha = I(i,j).*(0.9) + 0.1;
    end
    axis off; view([-90 0]);
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Correlation to neurobehavioral clinical scores
if(focus_idx > 2)
    focus = focus_idx - 2;
else
    focus = focus_idx;
end
%%%%%%%%%%%%%%%%%%%%%%% Original Code %%%%%%%%%%%%%%%%%%%%%%%%%
% load(['Results\' grp_idx '\scores.mat'],'Scores','idx_all'); %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Alternative Code %%%%%%%%%%%%%%%%%%%%%%
load(['Results\' grp_idx '\scores.mat'],'idx_all');          %%
Scores = [4*rand(size(idx_all,1),2)-2, ...                   %%              
    randi([5 20],size(idx_all,1),5)];                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
score = Scores(idx_all(:,subset),focus_idx);
load(['Results\' grp_idx '\Selection\subset_' num2str(subset) '_selected'], ...
    'K_pos','K_neg','r_pos','r_neg','p_pos','p_neg','x_pos','x_neg');
m_pos  = mean(r_pos{focus},2,'omitnan')';
m_neg  = mean(r_neg{focus},2,'omitnan')';
pp_pos = mean(p_pos{focus},2,'omitnan')';
pp_neg = mean(p_neg{focus},2,'omitnan')';
s_pos  = std(r_pos{focus},0,2,'omitnan')';
s_neg  = std(r_neg{focus},0,2,'omitnan')';
x_pos  = x_pos(focus);
x_neg  = x_neg(focus);

figure('Color',[1,1,1],'Position',[50 50 550 350]);
colororder({'k','k'}); yyaxis left;
hold on; grid on; box on; t_shift = 0.05; cnt = 1;
for j = 1:size(m_pos,2)
    t(cnt) = 1 + (j-1)*t_shift-2*t_shift;
    h = stem(t(cnt),100.*m_pos(j),'filled','Color',ccc(j,:), ...
        'LineStyle','-','LineWidth',33,'Marker','none');
    p1(j) = h(1);
    stem(t(cnt),100.*m_neg(j),'filled','Color',ccc(j,:), ...
        'LineStyle','-','LineWidth',33,'Marker','none');
    errorbar(t(cnt),100.*m_pos(j),0,100.*2*s_pos(j),'LineStyle','none', ...
        'Color',ccc(j,:),'LineWidth',3,'CapSize',20,'Marker','none');
    errorbar(t(cnt),100.*m_neg(j),100.*2*s_neg(j),0,'LineStyle','none', ...
        'Color',ccc(j,:),'LineWidth',3,'CapSize',20,'Marker','none');
    cnt = cnt + 1;
end
set(gca,'fontweight','bold','FontSize',22,'XGrid','off');
axis([0.85 1.15 -80 80]);
ylabel('Correlation (\%)','fontsize',26,'interpreter','latex');
set(gca,'xtick',1,'xticklabel',[],'TickLabelInterpreter','latex','YTick',-80:20:80);
a = gca; a.XAxis.FontSize = 22;
yyaxis right; hold on;
y = pp_pos;
stem(t,y,'Color','k','LineStyle','-','LineWidth',24,'Marker','none');
y = -1.*pp_neg;
h = stem(t,y,'Color','k','LineStyle','-','LineWidth',24,'Marker','none');
p1(end+1) = h(1); yLIM = 0.02;
axis([0.85 1.15 -1*yLIM yLIM]);
ylabel('Averaged p-Value','Interpreter','latex','FontSize',26);
set(gca,'xtick',1,'xticklabel',[], ...
    'TickLabelInterpreter','latex','YTick',-yLIM:yLIM/4:yLIM,...
    'YTickLabel',[yLIM:-yLIM/4:0 yLIM/4:yLIM/4:yLIM]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[50 50 950 550]);
[x, I] = sort(score);
y_pos = x_pos{1}; y_neg = x_neg{1};
% Averaging
y_pos = mean(y_pos,[1 3],'omitnan');
y_neg = mean(y_neg,[1 3],'omitnan');
% Extraction
y_pos = y_pos(:,I);
y_neg = y_neg(:,I);
% Linear regression
[P_pos, S_pos] = polyfit(x,y_pos,1);
[yfit_pos, delta_pos] = polyval(P_pos,x,S_pos);
[P_neg, S_neg] = polyfit(x,y_neg,1);
[yfit_neg, delta_neg] = polyval(P_neg,x,S_neg);
% Plotting
colororder([1 0.41 0.16; 0.07 0.62 1]);
hold on; grid on; grid minor; box on;
yyaxis left;
plot(x,y_pos,'Marker','o','MarkerSize',12,'linestyle','none', ...
    'MarkerFaceColor','none', 'MarkerEdgeColor',[1 0.41 0.16],'linewidth',1);
plot(x,yfit_pos,'-','LineWidth',3,'Color',[1 0.41 0.16]);
plot(x,yfit_pos+2.*delta_pos,'--','Color',[1 0.41 0.16],'linewidth',2);
plot(x,yfit_pos-2.*delta_pos,'--','Color',[1 0.41 0.16],'linewidth',2);
axis([min(x)-0.2, max(x)+0.2 2*min(y_pos)-max(y_pos)-0.01 max(y_pos)+0.01]);
yyaxis right;
plot(x,y_neg,'Marker','^','MarkerSize',12,'linestyle','none', ...
    'MarkerFaceColor','none', 'MarkerEdgeColor',[0.07 0.62 1],'linewidth',1);
plot(x,yfit_neg,'-','LineWidth',3,'Color',[0.07 0.62 1]);
plot(x,yfit_neg+2.*delta_neg,'--','Color',[0.07 0.62 1],'linewidth',2);
plot(x,yfit_neg-2.*delta_neg,'--','Color',[0.07 0.62 1],'linewidth',2);
axis([min(x)-0.2, max(x)+0.2 min(y_neg)-0.01 2*max(y_neg)-min(y_neg)+0.01]);
set(gca,'FontSize',22);
xlabel('Clinical Score','interpreter','latex','fontsize',26);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Results\Figures\Fig_1_EEG','-dpdf','-r400');
    print(2,'Results\Figures\Fig_1_Filter','-dpdf','-r400');
    print(3,'Results\Figures\Fig_1_multi_freq_EEG','-dpdf','-r400');
    print(4,'Results\Figures\Fig_1_multi_freq_parcel','-dpdf','-r400');
    print(5,'Results\Figures\Fig_1_subject_mdFCN','-dpdf','-r400');
    print(6,'Results\Figures\Fig_1_cohort_mdFCN','-dpdf','-r400');
    print(7,'Results\Figures\Fig_1_latent_mdFCN','-dpdf','-r400');
    print(8,'Results\Figures\Fig_1_decomposed_mdFCN_1','-dpdf','-r400');
    print(9,'Results\Figures\Fig_1_decomposed_mdFCN_2','-dpdf','-r400');
    print(10,'Results\Figures\Fig_1_decomposed_mdFCN_3','-dpdf','-r400');
    print(11,'Results\Figures\Fig_1_decomposed_mdFCN_4','-dpdf','-r400');
    print(12,'Results\Figures\Fig_1_decomposed_mdFCN_5','-dpdf','-r400');
    print(13,'Results\Figures\Fig_1_decomposed_mdFCN_6','-dpdf','-r400');
    print(14,'Results\Figures\Fig_1_selected_mdFCN_1','-dpdf','-r400');
    print(15,'Results\Figures\Fig_1_selected_mdFCN_2','-dpdf','-r400');
    print(16,'Results\Figures\Fig_1_selected_mdFCN_3','-dpdf','-r400');
    print(17,'Results\Figures\Fig_1_selected_mdFCN_4','-dpdf','-r400');
    print(18,'Results\Figures\Fig_1_selected_mdFCN_5','-dpdf','-r400');
    print(19,'Results\Figures\Fig_1_selected_mdFCN_6','-dpdf','-r400');
    print(20,'Results\Figures\Fig_1_correlation','-dpdf','-r400');
    print(21,'Results\Figures\Fig_1_linear_regression','-dpdf','-r400');
end