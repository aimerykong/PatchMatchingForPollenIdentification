% This code demonstrates the proposed objective function for discovering
% the most informative receptive fields in images.
% Specifically, a synthetic dataset is generated, including three cluters
% to stand for three images, and the most correlated image regions or the
% common foreground objects are represented by the points lying in the
% intersection of the three clusters.
%
%
% Readers can directly run main_syntheticData.m to see the results, and the
% intermediate results/figures are automatically stored in 'figures' folder.
%
%
% For details, readers are suggested to refer to the following report:
%       Shu Kong, "Collaborative Receptive Field Learning", arXiv, 2014
%
% 
% The code is writen by
%           Shu Kong (Aimery)
%           aimerykong@gmail.com
%           www.aimerykong.me
%           Dec. 2013, release version available on Jan. 13, 2014


