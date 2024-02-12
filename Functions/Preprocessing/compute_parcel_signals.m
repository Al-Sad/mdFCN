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
% It runs the preprocessing stage in Fig. 1.

function [parcels_AS, parcels_TA] = compute_parcel_signals(m,fs,T,data_folder)
%% Loading all requirements
load('Head Model\Atlas.mat','Atlas');
load('Head Model\InverseOperator.mat','InverseOperator');
load('Head Model\CollapseOperator.mat','CollapseOperator');
flt_eeg = iir_filters(fs);

%% Compute TA parcel signals
if(isfile([data_folder '\Subj_' sprintf('%03d',m) '\TA_eeg.mat']))
    % Loading and collecting TA epochs
    load([data_folder '\Subj_' sprintf('%03d',m) '\TA_eeg.mat'],'data');
    eeg_org = data.eeg';
    fs_org = data.Fs;
    clear data;
    % Preprocess
    eeg = pre_process(eeg_org,fs,fs_org,T,flt_eeg);
    % Compute parcel signals
    parcels_TA = get_parcel_signals(eeg,InverseOperator,CollapseOperator,Atlas);
else
    parcels_TA = [];
end
%% Compute AS parcel signals
if(isfile([data_folder '\Subj_' sprintf('%03d',m) '\AS_eeg.mat']))
    % Loading and collecting AS epochs
    load([data_folder '\Subj_' sprintf('%03d',m) '\AS_eeg.mat'],'data');
    eeg_org = data.eeg';
    fs_org = data.Fs;
    clear data;
    % Preprocess
    eeg = pre_process(eeg_org,fs,fs_org,T,flt_eeg);
    % Compute parcel signals
    parcels_AS = get_parcel_signals(eeg,InverseOperator,CollapseOperator,Atlas);
else
    parcels_AS = [];
end
end

%% Pre-process function
function y = pre_process(x,fs,fs_org,T,flt)
% Resampling
if(fs ~= fs_org)
    [p,q] = rat(fs/fs_org, 0.0001);
    y = resample(x,p,q);
    y(1:10,:) = []; y(end-10+1:end,:) = [];
else
    y = x;
end
% Take T seconds duration
y = y(1:(T*fs),:);
% Get the 5 frequency bands
y = filter_eeg(y, flt);
end

%% Anton's functions
function signal_src = get_parcel_signals(signal, InverseOperator, CollapseOperator, MyAtlas)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% see also Tokariev et al., 2019, Cerebral Cortex
%
% Input:
% signal: cell array {1 x N freq.} of filtered EEG data [samples x ch]
% InverseOperator: Inverse solution for infant head with N = 19 trodes
% CollapseOperator: weights for src signals within 'host' parcels
% MyAtlas: assignment of src/verticies to parcels (in MyAtlas.Parcels)
%
% Output:
% signal_src: cell array {1 x N freq.} of filtered parcel signals [samples x ch]

N_fr = size(signal, 2);       % number of frequencies

Np = length(MyAtlas.Parcels); % number of parcels

L = size(signal{1, 1}, 1);    % EEG length

signal_src{1, N_fr} = [];     % init output array

CollapseOperator = repmat(CollapseOperator, 1, L);

for k = 1:N_fr

    % source signals
    src_buf = (InverseOperator * signal{1, k}').* CollapseOperator;     % [sources x samples]

    % parcel/cortical signals
    parcel_buf = zeros(Np, L);

    for j = 1:Np
        parcel_buf(j, :) = mean(src_buf(MyAtlas.Parcels{j, 1}, :));     % [parcels x samples]
    end

    signal_src{1, k} = parcel_buf';                                     % [samples x parcels]

end

end
function signal_flt = filter_eeg(signal, flt)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Input:
% signal: [samples x ch]
%    flt: cell array with filter objects (top - HPF, bottom - LPF)
%    Fs: sampling rate
%
% Output:
% signal_flt: cell array{1, N of fr.bands} of filtered data [samples x ch]

% Signal Length
L = size(signal, 1);

% Flipped signal
signal_ud = flipud(signal);

% Add pieces to signal (to fix edge effects)
signal_big = [signal_ud; signal; signal_ud];

% Number of bandpass filters (= fr.bands); columns = bandpass filters
N_flt = size(flt, 2);

% Init cell array for bandpass filtered signals
signal_flt{1, N_flt} = [];

% Filter signals (bandpass filter = HPF + LPF)
for fr = 1:N_flt

    buf = [];                                                            %#ok<NASGU>
    buf = filtfilt(flt{1, fr}, signal_big); % HPF/cutoff = 0.85xFc
    buf = filtfilt(flt{2, fr}, buf);        % LPF/cutoff = 1.15xFc

    signal_flt{1, fr} = buf(L+1:L+L, :);    % cut signal_big >> orig.sig.

end
end