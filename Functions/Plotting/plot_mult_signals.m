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
% It plots multichannel signals at different frequencies.

function plot_mult_signals(z, fs, lables, type)
hold on;
L = size(z{1}, 1); % length
ch = size(z{1}, 2); % number of channels
shift = 5*(max(std(z{1}))); % b/w different channels for visualization
t = (0:1:L-1)/fs;
t = repmat(t, ch, 1);
t = t';
if strcmp(type, 'eeg')
    ccc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; ...
    0.494 0.184 0.556; 0.929 0.694 0.125];
else
    ccc = zeros(ch, 3); % Occipital = black (default)
    ccc(cell2mat(lables) == 'T', 1:3) = repmat([0.47 0.67 0.19],sum(cell2mat(lables) == 'T'), 1); % Temporal
    ccc(cell2mat(lables) == 'C', 1:3) = repmat([0.75 0.00 0.75],sum(cell2mat(lables) == 'C'), 1); % Central
    ccc(cell2mat(lables) == 'F', 1:3) = repmat([0.93 0.69 0.13],sum(cell2mat(lables) == 'F'), 1); % Frontal
end

for i = 1:length(z)
    s = z{i};
    baseline = (1:1:ch)*shift;
    baseline = fliplr(repmat(baseline, L, 1));
    % add baseline to each channel for visualization
    sig_plot = s + baseline;
    % Plotting
    if strcmp(type, 'eeg')
        for j = 1:ch
            plot3(t(:,j), repelem(i,length(t)), sig_plot(:,j), 'Color', ccc(i,:));
        end
    else
        for j = 1:ch
            plot3(t(:,j), repelem(i,length(t)), sig_plot(:,j), 'Color', ccc(j,:));
        end
    end
    set(gca, 'ZTick', fliplr(baseline(1,:)));
end
% ylim([min(min(sig_plot)) max(max(sig_plot))]); hold off
end