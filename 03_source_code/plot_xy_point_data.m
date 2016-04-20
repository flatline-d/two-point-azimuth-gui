function [n_pts] = plot_xy_point_data(xy_pts_file_name)
%
% Function name = read_xy_point_data.m
% Author = Bradley Thomson (bjt@bu.edu)
% Created on 20 Oct 2012
% Last tinkered with on 28 Nov 2012
%
% The goal of this function is to import the x,y shield field position 
% data and plot the spatial distribution.
%
% 	Function inputs:
%       -xy_pts_file_name : file where each x,y point pair is the latitude
%       and longitude of each shield field. (.txt or .csv)
%
% 	Output:
%       -Scatterplot of shield locations
%       -Count of number of point pairs in input file
%
%   Sub-functions called:
%       -read_xy_point_data : read in designated input file and return
%       arrays of X and Y values



%% Read in data file
%

% Open x,y data file for reading using subfuction "read_xy_point_data.m"
%
[x_vals, y_vals] = read_xy_point_data(xy_pts_file_name);



%% Count number of points in input file
%
n_pts = size(x_vals,1);



%% Create output figure for visual check
%

plot(x_vals,y_vals,'^') % plot symbol is upward-pointing triangle
xlabel('Lon (°E)')
ylabel('Lat (°N)')
axis image  % forces aspect ratio to be one-to-one






