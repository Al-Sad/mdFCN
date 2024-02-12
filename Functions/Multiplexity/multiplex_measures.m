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
% It calculates various measures of structural multiplexity.

function measures = multiplex_measures(z,band_mask,fideltymask)
%% The multiplex network
N_feat = length(z);
A = cell(1,N_feat);
K = cell(1,N_feat);
for nf = 1:N_feat
    ML = sum(band_mask(nf,:));
    z_bar = squeeze(mean(z{nf},2,'omitnan')); % average across infants
    z_bar(isnan(z_bar)) = 0;                  % set insignificant connections to zero    
    z_bar = z_bar(band_mask(nf,:),:);         % get significant bands to zero
    for m = 1:ML
        [~, y] = fold_connectivity(z_bar(m,:),fideltymask,0);
        y = logical(y);
        A{nf}(m,:,:) = y;
        K{nf}(m,:) = sum(y,1);
    end
end

%% Node Structure
MD = length(y);
Oi = nan(N_feat, MD); % normalized overlapping degree
Pi = nan(N_feat, MD); % multiplex participation coefficient
for nf = 1:N_feat
    ML = sum(band_mask(nf,:));
    if(any(band_mask(nf,:)))
        Oi(nf,:) = sum(K{nf},1)./(ML*(MD-1));
        if(ML > 1)
            Pi(nf,:) = (ML/(ML-1))*(1 - sum((K{nf}./sum(K{nf},1)).^2,1));
        end
    end
end

%% Layer Structure
ML = size(band_mask,2);
S = nan(N_feat, ML, ML); % pairwise multiplexity
H = nan(N_feat, ML, ML); % Hamming similarity
for nf = 1:N_feat
    ML = sum(band_mask(nf,:));
    idx = find(band_mask(nf,:));
    bi = logical(K{nf});
    for m1 = 1:ML
        for m2 = 1:ML
            S(nf,idx(m1),idx(m2)) = sum(bi(m1,:).*bi(m2,:),2)./MD;
            num = sum(bi(m1,:).*(1-bi(m2,:)) + bi(m2,:).*(1-bi(m1,:)),2);
            den = min([sum(bi(m1,:)+bi(m2,:),2) MD]);
            H(nf,idx(m1),idx(m2)) = num./den;
        end
    end
end

%% Edge Structure
Oe = nan(N_feat,1);
X  = nan(N_feat,1);
ME = (MD^2-MD)/2;
for nf = 1:N_feat
    ML = sum(band_mask(nf,:));
    Oe(nf) = sum(K{nf},'all')./(ML.*ME);
    if(ML > 1)
        X(nf) = ML.*sum(min(A{nf},[],1),'all')./sum(K{nf},'all');
    end
end

%% Output
measures.Oi = Oi;
measures.Pi = Pi;
measures.S  = S;
measures.H  = H;
measures.Oe = Oe;
measures.X  = X;

end