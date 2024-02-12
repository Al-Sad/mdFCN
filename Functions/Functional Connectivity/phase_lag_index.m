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
% It computes the PLI measures of connectivity.

function [pli, wpli, dwpli] = phase_lag_index(z,mask)
if(~ismatrix(z))
    x = z;
else
    x(1,:,:) = z;
end
[Nw, L, Nc] = size(x);
pli   = zeros(Nw,Nc,Nc);
wpli  = zeros(Nw,Nc,Nc);
dwpli = zeros(Nw,Nc,Nc);
for n = 1:Nw
    zt = squeeze(x(n,:,:));
    zf = fft(zt,2*L,1);
    zf = zf(1:end/2,:);
    for p = 1:(Nc-1)
        phi1 = angle(hilbert(zt(:,p)));
        for q = (p+1):Nc
            if(mask(p,q))
                phi2 = angle(hilbert(zt(:,q)));
                pli(n,p,q) = abs(mean(sign(imag( exp(1i*(phi1-phi2)) ))));
                Sxy_img = imag(zf(:,p).*conj(zf(:,q)));
                wpli(n,p,q) = abs(mean(Sxy_img)./mean(abs(Sxy_img)));
                num = sum(Sxy_img).^2 - sum(Sxy_img.^2);
                den = sum(abs(Sxy_img)).^2 - sum(Sxy_img.^2);
                dwpli(n,p,q) = abs(num./den);
                pli(n,q,p)   = pli(n,p,q);
                wpli(n,q,p)  = wpli(n,p,q);
                dwpli(n,q,p) = dwpli(n,p,q);
            end
        end
    end
end
pli   = squeeze(pli);
wpli  = squeeze(wpli);
dwpli = squeeze(dwpli);
end