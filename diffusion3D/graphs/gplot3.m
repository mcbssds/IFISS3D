function h = gplot3(A,xyz,linespec)
%GPLOT3 plots edge-node graph specified by adjacency and coordinate matrices
%  h = gplot3(A,xyz,linespec) 
%  inputs:
%         A          matrix with 1 in (i,j)th position if nodes i and j 
%                    are connected and 0 otherwise
%         xyz        matrix containing (x,y,z) coordinates of nodes
%         linespec   specifies line type and/or colour for plotting; see
%                    MATLAB function PLOT for possibilities
%
%   output:
%         h          handle for created plot
%
% IFISS function CP; 11 August 2022.
% Modified from code available at: https://github.com/cmccomb/gplot3
% which is itself a modified version of the in-built MATLAB gplot function

% Returns i and j, lists of connected nodes
[i,j] = find(A);
% Extract
X = [xyz(i,1) xyz(j,1)]';
Y = [xyz(i,2) xyz(j,2)]';
Z = [xyz(i,3) xyz(j,3)]';
% Add NaN values to break between line segments
X = [X; NaN(size(i))'];
Y = [Y; NaN(size(i))'];
Z = [Z; NaN(size(i))'];
% Serialize the x, y and z data
X = X(:);
Y = Y(:);
Z = Z(:);

h = plot3(X,Y,Z,linespec);
