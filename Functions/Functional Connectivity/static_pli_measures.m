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
% It computes the static PLI measures of connectivity for different frequency bands.

function [PLI, wPLI, dwPLI] = static_pli_measures(x,mask)
if(~isempty(x))
    M     = (size(x{1},2)^2 - size(x{1},2))/2;
    PLI   = zeros(size(x,2), M);
    wPLI  = zeros(size(x,2), M);
    dwPLI = zeros(size(x,2), M);
    for i = 1:size(x,2)
        [pli, wpli, dwpli] = phase_lag_index(x{i},mask);
        PLI(i,:)   = unfold_upper_matrix(pli);
        wPLI(i,:)  = unfold_upper_matrix(wpli);
        dwPLI(i,:) = unfold_upper_matrix(dwpli);
    end
else
    PLI   = [];
    wPLI  = [];
    dwPLI = [];
end
end