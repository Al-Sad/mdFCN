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
% It plots multichannel signals.

function plot_signals(s, lables, Fs, type)
hold on;

L = size(s, 1); % length
ch = size(s, 2); % number of channels
shift = 5*(max(std(s))); % b/w different channels for visualization
% shift = 75;
baseline = (1:1:ch)*shift;
baseline = fliplr(repmat(baseline, L, 1));

% time
t = (0:1:L-1)/Fs;
t = repmat(t, ch, 1);
t = t';

% add baseline to each channel for visualization
sig_plot = s + baseline;

% Plot EEG/Filtered EEG (all black)
if strcmp(type, 'eeg') == 1
    for i = 1:ch
        plot(t(:,i), sig_plot(:,i), 'Color', 'k');
    end
end

% Plot parcel signals (colorful)
if strcmp(type, 'src') == 1
    % assign parcel colors
    clr = zeros(ch, 3); % Occipital = black (default)
    clr(cell2mat(lables) == 'T', 1:3) = repmat([0.47 0.67 0.19],sum(cell2mat(lables) == 'T'), 1); % Temporal
    clr(cell2mat(lables) == 'C', 1:3) = repmat([0.75 0.00 0.75],sum(cell2mat(lables) == 'C'), 1); % Central
    clr(cell2mat(lables) == 'F', 1:3) = repmat([0.93 0.69 0.13],sum(cell2mat(lables) == 'F'), 1); % Frontal
    for i = 1:ch
        plot(t(:,i), sig_plot(:,i), 'Color', clr(i, :));
    end
end
set(gca, 'YTick', fliplr(baseline(1,:)));
ylim([min(min(sig_plot)) max(max(sig_plot))]); hold off
end