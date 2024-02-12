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
% It plots the brain atlas and the cortical regions in 3D.

function plot3d_atlas(Atlas, flat_cx, sz)
hold on;
% plot colour coded parcel centroids
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 1), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 2), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 3), sz, 'ok', 'filled');
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 1), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 2), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 3), ...
    sz, 'MarkerFaceColor', [0.47 0.67 0.19], ...
    'MarkerEdgeColor', [0.47 0.67 0.19]);
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 1), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 2), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 3), ...
    sz, 'MarkerFaceColor', [0.75 0.00 0.75], ...
    'MarkerEdgeColor', [0.75 0.00 0.75]);
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 1), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 2), ...
    Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 3), ...
    sz, 'MarkerFaceColor', [0.93 0.69 0.13], ...
    'MarkerEdgeColor', [0.93 0.69 0.13]);

% plot smoothed cx
patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, ...
    'FaceColor', [0.80 0.80 0.80], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlim([-0.06 0.06]);
ylim([-0.05 0.05]);
zlim([0 0.08]);
axis off; hold off;
end