function [x_vals, y_vals] = read_xy_point_data(xy_pts_file_name)
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
%       -2 arrays, one for x-values and one for y-values



%% Read in data file
%

% Open x,y data file for reading
%
fprintf('The x,y point file name is %s\n', xy_pts_file_name) % display file name for debugging
%
fid = fopen(xy_pts_file_name);  % open data file for reading
A = textscan(fid, '%f %f', 'delimiter', ',');  % open ascii data into array using textscan
fclose(fid);


%% Separate cell structure into separate x and y arrays
%
x_vals = [A{:,1}];
y_vals = [A{:,2}];







