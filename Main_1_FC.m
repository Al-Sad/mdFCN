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
% It calculates the static and dynamic FC of each neonate's EEG.
% Note that the neonates raw EEG and neurocognitive scores are not
% supplied.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
grp_idx = 'HC'; % HC or AED
T       = 3*60; % Parcel signal duration in seconds
fs      = 100;  % Sampling frequency in hertz
L       = 2;    % Segment window duration in seconds
P       = 0.5;  % Segment window overlap in % of L

%% Main
data_folder = ['Dataset\' grp_idx '\'];
load('Head Model\FidelityOperator.mat');
mask = logical(FidelityOperator);
for m = 1:(length(dir(data_folder))-2)
    disp(['Subject: ' num2str(m)]);
    [parcels_AS, parcels_TA] = compute_parcel_signals(m,fs,T,data_folder);
    % Compute static dwpli
    [~, ~, sdwPLI_AS] = static_pli_measures(parcels_AS,mask);
    [~, ~, sdwPLI_TA] = static_pli_measures(parcels_TA,mask);
    % Segment the parcel signals using L and P
    parcels_seg_AS = parcel_seg(parcels_AS,L*fs,P);
    parcels_seg_TA = parcel_seg(parcels_TA,L*fs,P);
    % Compute dynamic dwpli
    [~, ~, dwPLI_AS] = dynamic_pli_measures(parcels_seg_AS,mask);
    [~, ~, dwPLI_TA] = dynamic_pli_measures(parcels_seg_TA,mask);
    % Saving
    save(['Results\' grp_idx '\FC\Subj_' num2str(m) ...
        '.mat'],'sdwPLI_AS','sdwPLI_TA','dwPLI_AS','dwPLI_TA');
end