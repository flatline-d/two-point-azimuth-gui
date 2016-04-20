function [az2, xy_endpts] = two_pt_azimuth_cebria_v3_mt(x_vals, y_vals, planet)
%
% Goal: The point of this program is to implement the J. M. Cebria's 2011 
% alternative version of 2-point azimuth method, which is a method of
% looking for preferred orientations in a field of point features. In the
% original method, the azimuth or angle between all N points are tabulated,
% for a total of N*(N-1)/2 azimuth values.
%
%   In this method, only azimuth values between pairs of points that are
% less than or equal to one standard deviation of the mean separation
% value divided by 3, i.e., d <= (X - 1_sigma)/3. Filtering the azimuth
% results by this distance criterion or threshold effectively examines only
% short-range order within a given field of points.
%
% Author: Brad Thomson
% Last tinkered with: 7 Apr 2016
%
% 	Function inputs:
%       -2 arrays, one for x-values and one for y-values
%       -planetary ID string (e.g., Venus, Mars, Earth)
%
% 	Output:
%       -az2 : for N numbers of x,y point pairs, return all azimuth values 
%       that are less than or equal to (X - 1_sigma)/3, where X is the mean
%       separation distance and sigma is the standard deviation. Note the
%       number of azimuth values that meet this criterion is not known in
%       advance.
%       -xy_pts : a m-by-4 matrix of pairs of x,y coordinates that
%       represent the start and end points of line segmenets that are less
%       than the Cebria criterion.
%
%   Sub-functions called:
%       -referenceEllipsoid : From mapping toolbox, used to define a  
%        planetary ellipsoid with specific radii.
%       -distance : From mapping toolbox, calculates distance and azimuth
%       between pair of lat, lon points.
%
% Change list:
% 29 Dec 2012: swapped function "azimuth" for "maptool_azimuth" stored in
% local directory (mapping toolbox needed for former function).
%
% 2013-08-25: combined calculation of distance and azimuth for
% computational efficieny using just function 'maptool_distance'
%
% 2015-04-18: added code to output x,y coordinates of pairs of features
% that meet Cebria maximum distance criterion. 
%
% 2015-08-04: Updated to use Mapping Toolbox functions.
%
% 2016-04-07: Replaced referenceSphere with referenceEllipsoid.


%% Initialize variables, counters, and results matricies 
%
N_lines = size(x_vals,1); % # of lines of input text = number of xy pairs
N_dists = (N_lines*(N_lines -1))/2; % # of distances that are calculated
N_angles = (N_lines*(N_lines -1))/2; % # of azimuth angles that are calculated
az=zeros(N_angles,1); % pre-allocate and initialize azimuth results array
dist=zeros(N_dists,1); % pre-allocate and initialize distance results array
dist2=zeros(N_dists,1); % pre-allocate and initialize distance holder array
xy_ind_start=zeros(N_lines,1); % initialize start indicies holder array
xy_ind_end=zeros(N_lines,1); % initialize end point indicies holder array
i = 1; % initialize counter for 1st xy pair
j = 2; % initialize counter for 2nd xy pair (not necessary; here for completeness)
k = 1; % initialize counter for number of distances or angles calculated
dist_mean = 0; % initialize mean distance value
dist_stdev = 0; % initialize standard deviation of mean distance value
dist_thresh = 0; % intialize distance threshold value
%
% Set up reference ellipsoid S based on planet selected, units in kilometers
fprintf('Selected planet is %s \n', planet); 
S = referenceEllipsoid(planet,'km');



%% Compute distance and azimuth values
%
for i = 1:(N_lines-1)

	for j = i:(N_lines-1)

        % Use built-in MATLAB function distance, which measures the
        % distance between points on a sphere or ellipsoid and also
        % computes the azimuth value. Azimuth values are measured clockwise
        % from the north in default units of degrees (0 to 360°).
        %
        % Useage: [dist,az]=distance(lat1,lon1,lat2,lon2,planet)
        %             =distance(x-y_vals(row_i), x-y_vals(row_j+1), planet)
%         [dist(k,1),az(k,1)]=maptool_distance(y_vals(i,1),x_vals(i,1),y_vals(j+1,1),x_vals(j+1,1),S);
        [dist(k,1),az(k,1)]=distance(y_vals(i,1),x_vals(i,1),y_vals(j+1,1),x_vals(j+1,1),S);
%         fprintf('x_vals(j+1,1) is %.2f\n',x_vals(j+1,1)) % display dist value for debugging
%         fprintf('azimuth is %.2f\n',az(k,1)) % display azimuth value for debugging
        %
        xy_ind_start(k) = i; % assign starting point index to holder array
        xy_ind_end(k) = (j + 1); % assign end point index to holder array
        %
		j = j + 1; % increment 2nd xy pair index
		k = k + 1; % increment distance/angle counter index

	end % end inner FOR loop

	i = i + 1; % increment 1st xy pair index

end % end outer FOR loop



% Calculate mean distance value and standard deviation of mean
dist_mean = mean(dist,1);
dist_stdev = std(dist,1);
%
% Calculate distance threshold value = |(mean-std_dev)|/3
dist_thresh = abs(dist_mean-dist_stdev)/3;
%
% Display values for debugging
fprintf('The mean separation distance is %.2f km\n',dist_mean)
fprintf('The standard deviation of the mean is %.2f km\n',dist_stdev)
fprintf('The threshold value (mean-stdev)/3 is %.2f km\n',dist_thresh)



%% Check which distance values are < threshold, tabulate azimuth values
% for just these point pairs.
%
k = 1; % re-initialize counter for number of distance calculated
m = 1; % initialize 2nd counter for number of angles that pass threshold
%
for k = 1:(size(dist,1))

    % Step though each distance value (total of N*(N-1)/2)
	if dist(k,1)<=dist_thresh       
        %
        % If condition met, assign this azimuth value to a second array
        az2(m,1)=az(k,1);
        xy_endpts(m,1) = x_vals(xy_ind_start(k,1)); % x-value of start pt
        xy_endpts(m,2) = y_vals(xy_ind_start(k,1)); % y-value of start pt
        xy_endpts(m,3) = x_vals(xy_ind_end(k,1)); % x-value of end pt
        xy_endpts(m,4) = y_vals(xy_ind_end(k,1)); % y-value of end pt
        % Note: matricies 'az2' and 'xy_endpts' are dynamically grown.
        % Re-write this code to avoid this and improve speed 
        %
        % (1) convert 180°-270° to 0°-90°
        if az2(m,1) >= 180 && az2(m,1) < 270
            az2(m,1) = az2(m,1) - 180;
        end % end IF statement
        %
        % (2) convert 270-360° to (0° to -90°)
        if az2(m,1) >= 270 && az2(m,1) < 360
            az2(m,1) = az2(m,1) - 360;
        end % end IF statement
        %
        % (3) convert 90-180° to (0° to -90°)
        if az2(m,1) >= 90 && az2(m,1) < 180
            az2(m,1) = az2(m,1) - 180;
        end % end IF statement
        %
        %
        m = m + 1; % increment counter of angles that meet threshold

    end % end IF statement for distance check threshold

    k = k + 1; % increment distance counter index
    
end % end FOR loop

% Debugging: Display number of azimuths calculated for debugging
num_az = size(az2,1); % debugging: Sum up the number of values in histogram
fprintf('The number of azimuths calculated in this run is %.0f\n\n', num_az);


