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
% Calculates the positive/negative correlation density.

function [Kp, Kn] = corr_density(r,p,alpha)
N = size(r,1);
[Nb, Nc] = size(r{1});
Kp = zeros(N,Nb);
Kn = zeros(N,Nb);
for i = 1:N
    for j = 1:Nb
        pvals = p{i}(j,:); pvals(isnan(pvals)) = 1;
        It = ~isnan(pvals) & make_FDR(pvals,alpha);
        I_pos = (r{i}(j,:) > 0) & It;
        I_neg = (r{i}(j,:) < 0) & It;
        Kp(i,j) = nnz(I_pos)/Nc;
        Kn(i,j) = nnz(I_neg)/Nc;
        Rp(i,j) = mean(r{i}(j,I_pos),'omitnan');
        Rn(i,j) = mean(r{i}(j,I_neg),'omitnan');
    end
end
end