function [az_hist_values, bin_centers] = azimuth_method_v3_mt(...
    xy_pts_file_name, axes2_handle, axes3_handle, planet)
%
% Goal: The point of this program is to create 2 plots of the results of
% the two point azimuth method [Lutz, 1986]. Here, the azimuth or angle
% between all points n are tabulated, for a total of N = n(n-1)/2 azimuth
% values.
%
% Author: Brad Thomson (bjt@bu.edu)
% Last tinkered with: 26 Aug 2013
%
% 	Function inputs:
%       -xy_pts_file_name : file where each x,y point pair is the latitude
%       and longitude of each shield field. (.txt or .csv)
%       -axes2_handle : handle of axes to plot histogram
%       -axes3_handle : handle of axes to plot rose plot
%       -planet : string that indentifies the planetary body of interest
%
% 	Output:
%       -Histogram of azimuth values
%       -Rose plot of azimuth values
%       -Observed azimuth histogram values as function output
%       -Bin center values of histogram as function output
%
%   Sub-functions called:
%       -read_xy_point_data : read in designated input file and return
%       arrays of X and Y values
%       -two_pt_azimuth_calcs_v3_mt : calculate 2-point azimuth values
%
% Change list:
%
% 29 Dec 2012: Changes from version1, Subfunction "two_pt_azimuth_calcs_v2"
% called instead of "two_pt_azimuth_calcs"
%
% 26 Aug 2013: Added bin center values as second direct output of this
% function.
%
% 29 Dec 2012: Changes from version3, Subfunction "two_pt_azimuth_calcs_v3"
% called instead of "two_pt_azimuth_calcs_v2"
%
% 4 Aug 2015: Changed to subfunction "two_pt_azimuth_calcs_v3_mt"
% instead of "two_pt_azimuth_calcs_v3."


% Initialize variables, counters, and results matricies 
bin_edges = -90:10:90; % create 18 histogram bins from -90 to 90 (10° bins)



%% Read in data file
%

% Open x,y data file for reading using subfuction "read_xy_point_data.m"
%
[x_vals, y_vals] = read_xy_point_data(xy_pts_file_name);



%% Compute 2-point azimuth values using subfunction "two_pt_azimuth_calcs_v3"
%

[az] = two_pt_azimuth_calcs_v3_mt(x_vals, y_vals, planet);



%% Plot results
%
axes(axes2_handle) % designate active axes
% bin_centers = 5:10:175; % 10 deg bins, n = 18 between 0-180 degrees
bin_centers = -85:10:85; % 10 deg bins, n = 18 between -90 to 90 degrees
% Plot azimuth values in 2D histogram. Recall N is 0 deg increasing
% clockwise (90° = east)
hist(az, bin_centers)
% set(gca, 'XTick',0:30:180, 'XLim',[0 180]) % set X-axis limit to 0-180 with 30° ticks
set(gca, 'XTick',-90:30:90, 'XLim',[-90 90]) % set X-axis limit to -90 to +90 with 30° ticks
xlabel('Azimuth (°E)'); % label x-axis of plot
ylabel('Frequency'); % label y-axis of plot


axes(axes3_handle) % re-designate second active axes
az_rad = az*pi/180; % convert azimuth values from decimal to radians
% Plot azimuth values in rose plot
%   Useage: rose(theta, nbins)
%   theta is array of vectors (in radians)
%   nbins is the number of angle bins used (default is 20)
h2 = rose(az_rad, 36);
x = get(h2,'Xdata');
y = get(h2,'Ydata');
g=patch(x,y,'b');
% Rotate and mirror rose plot so that north (0 deg) is at top; east (+90
% deg) is at right.
set(gca,'View',[-90 90],'YDir','reverse');



%% Assign observed azimuth histogram to function output
%
n = histc(az,bin_edges); % n is array with number of values in each bin
                         % note last bin counts values equal to bin edge
% Assign results to output array, removing last bin
az_hist_values = n(1:(size(n)-1));




