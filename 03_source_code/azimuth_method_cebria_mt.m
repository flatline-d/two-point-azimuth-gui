function [az_hist_values, bin_centers] = azimuth_method_cebria_mt(xy_pts_file_name, axes1_handle, axes2_handle, axes3_handle, planet)
%
% Goal: The point of this program is to create 2 plots of the results of
% a modification of the  two point azimuth method [Lutz, 1986] developed by
% Cebria. Here, the azimuth or angle between all points n that are less
% than a threshold value of (mean - standard_deviation)/3 are tabulated.
%
% The exact number of azimuths thus calculated is not known in advance, but
% it is typically <5% of the number in the standard two-point azimuth
% method.
%
% Author: Brad Thomson (bjt@bu.edu)
% Last tinkered with: 18 Apr 2015
%
% 	Function inputs:
%       -xy_pts_file_name : file where each x,y point pair is the latitude
%       and longitude of each shield field. (.txt or .csv)
%       -axes2_handle : handle of axes to plot histogram
%       -axes3_handle : handle of axes to plot rose plot'
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
%       -two_pt_azimuth_cebria_v3_mt : calculate 2-point azimuth values
%
% Change list:
% 14 Aug 2013: Changed subfunction "two_pt_azimuth_cebria"
% called instead of "two_pt_azimuth_calcs_v2"
%
% 23 Aug 2013: Added bin center values as second direct output of this
% function.
%
% 2015-04-18: Subfunction 'two_pt_azimuth_cebria_v3' called instead of
% 'two_pt_azimuth_cebria_v2'; new version outputs paris of x,y coordinates
% from points that are closer than the Cebria cutoff distance.



% Initialize variables, counters, and results matricies 
bin_edges = -90:10:90; % create 18 histogram bins from -90 to 90 (10° bins)



%% Read in data file
%

% Open x,y data file for reading using subfuction "read_xy_point_data.m"
%
[x_vals, y_vals] = read_xy_point_data(xy_pts_file_name);



%% Compute 2-point azimuth values using subfunction "two_pt_azimuth_cebria_v3_mt"
%

tic % start timer for debugging
[az, xy_pts] = two_pt_azimuth_cebria_v3_mt(x_vals, y_vals, planet);
toc % stop timer for debugging

% Process xy_pts into format suitable for plotting. Note that we need to
% append a 2-col by N_lines number of rows matrix of NaNs to the end of
% xy_pts. This creates a series of NaN in the 2-column matrix of x-y points
% to send to plot (every 3rd row is NaN NaN). This only plots a line to the
% segment in question, not connecting different segments.
%
% Find size of xy_pts matrix
N_lines = size(xy_pts,1);
% Create array N-by-2 of all ones, multiple by NaN
null_array = ones(N_lines,2);
null_array = NaN * null_array;
xy_pts_w_nulls = cat(2, xy_pts, null_array); % concatenate null array to
                                             % end of xy_pts
% Use 'reshape' to resort m-by-6 matrix into 3m-by-2 matrix
xy_pts_2col = reshape(xy_pts_w_nulls', [2 3*N_lines])';

                                            

%% Plot results
%
axes(axes1_handle) % designate active axes
% Linespec options: k=black, b=blue, '-'=solid line (default)
plot(xy_pts_2col(:,1), xy_pts_2col(:,2), 'k-');
axis image 


axes(axes2_handle) % designate next active axes
bin_centers = -85:10:85; % 10 deg bins, n = 18 between -90 to 90 degrees
% bin_centers = -82.5:15:82.5; % 15 deg bins, n = 12 between -90 to 90 degrees
% Plot azimuth values in 2D histogram. Recall N is 0 deg increasing
% clockwise (90° = east)
hist(az, bin_centers);
% set(gca, 'XTick',0:30:180, 'XLim',[0 180]) % set X-axis limit to 0-180 with 30° ticks
set(gca, 'XTick',-90:30:90, 'XLim',[-90 90]) % set X-axis limit to -90 to +90 with 30° ticks
xlabel('Azimuth (°E)'); % label x-axis of plot
ylabel('Frequency'); % label y-axis of plot


axes(axes3_handle) % re-designate third active axes
az_rad = az*pi/180; % convert azimuth values from decimal to radians
% Plot azimuth values in rose plot
%   Useage: rose(theta, nbins)
%   theta is array of vectors (in radians)
%   nbins is the number of angle bins used (default in 'rose' is 20)
%   nbins = 24 for 15° bins, 36 for 10° bins
% h2 = rose(az_rad, 24);
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
% Note: This seems inefficient as the histogram is computed twice...


