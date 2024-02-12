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
% Modified version of nmfnnls.m in \nmfv1_4.

function [A,Y]=nmfnnls_mod(X,k,option)
% NMF based on NNLS: X=AY, s.t. X,A,Y>=0.
% Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=nmfnnls(X,k)
%     [A,Y,numIter,tElapsed,finalResidual]=nmfnnls(X,k,option)
% X: non-negative matrix, dataset to factorize, each column is a sample, and each row is a feature.
% k: scalar, number of clusters.
% option: struct:
% option.iter: max number of interations. The default is 1000.
% option.dis: boolen scalar, It could be 
%     false: not display information,
%     true: display (default).
% option.residual: the threshold of the fitting residual to terminate. 
%     If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-4.
% option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 1e-4.
% A: matrix, the basis matrix.
% Y: matrix, the coefficient matrix.
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
%  [1]\bibitem{NMF_ANLS_Kim2008}
%     H. Kim and H. Park,
%     ``Nonnegative matrix factorization based on alternating nonnegativity constrained least squares and active set method,''
%     {\it SIAM J. on Matrix Analysis and Applications},
%     vol. 30, no. 2, pp. 713-730, 2008.
%  [2]\bibitem{NMF_Sparse_Kim2007}
%     H. Kim and H. Park,
%     ``Sparse non-negatice matrix factorization via alternating non-negative-constrained least squares for microarray data analysis,''
%     {\it Bioinformatics},
%     vol. 23, no. 12, pp. 1495-1502, 2007.
%%%%
% Copyright (C) <2012>  <Yifeng Li>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 01, 2011
%%%%

optionDefault.iter=200;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<3
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

% iter: number of iterations
[~,c] = size(X); % c is # of samples, r is # of features
Y=rand(k,c);
XfitPrevious=Inf;
for i=1:option.iter
    A=kfcnnls(Y',X');
    A=A';
    A=normc(A);
    Y=kfcnnls(A,X);
    if mod(i,20)==0 || i==option.iter
        XfitThis=A*Y;
        fitRes=matrixNorm(XfitPrevious-XfitThis);
        XfitPrevious=XfitThis;
        curRes=norm(X-XfitThis,'fro');
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            break;
        end
    end
end
end

