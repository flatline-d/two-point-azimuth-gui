function [az] = two_pt_azimuth_calcs_v3_mt(x_vals, y_vals, planet)
%
% Goal: The point of this program is to test the algorithm for the 2-point
% azimuth method, which is a method of looking for preferred orientations
% in a field of point features. Here, the azimuth or angle between all
% points N are tabulated, for a total of (N-1)! azimuth values.
%
% Author: Brad Thomson  (bjt@bu.edu)
% Last tinkered with: 7 Apr 2016
%
% 	Function inputs:
%       -2 arrays, one for x-values and one for y-values
%       -planetary ID string (e.g., Venus, Mars, Earth)
%
% 	Output:
%       -For n x,y point pairs, N = n*(n-1)/2 azimuth values between all n
%       points
%
%   Sub-functions called:
%       -referenceEllipsoid : From mapping toolbox, used to define a  
%        planetary ellipsoid with specific radii.
%       -azimuth : From mapping toolbox, calculated azimuth between pairs
%       of lat, lon points.
%
% Change list:
% 29 Dec 2012: swapped function "azimuth" for "maptool_azimuth" stored in
% local directory (mapping toolbox needed for former function).
%
% 2013-08-25: combined calculation of distance and azimuth for
% computational efficieny using just function 'maptool_distance'
%
% 2015-08-04: Updated to use Mapping Toolbox functions.
%
% 2016-04-07: Replaced referenceSphere with referenceEllipsoid.


%% Initialize variables, counters, and results matricies 
%
N_lines = size(x_vals,1); % number of lines of input text = number of xy pairs
N_angles = (N_lines*(N_lines -1))/2; % number of azimuth angles that will be calculated
az=zeros(N_angles,1); % pre-allocate and initialize results array
i = 1; % initialize counter for 1st xy pair
j = 2; % initialize counter for 2nd xy pair (not necessary; here for completeness)
k = 1; % initialize counter for number of angles calculated
%
% Set up reference ellipsoid S based on planet selected, units in kilometers
fprintf('Selected planet is %s \n', planet); 
S = referenceEllipsoid(planet,'km');



%% Compute azimuth values
%
for i = 1:(N_lines-1)

	for j = i:(N_lines-1)

        % Use built-in MATLAB function azimuth, which measures the azimuth
        % between points on a sphere or ellipsoid. Azimuth values are
        % measured clockwise from the north in default units of degrees (0
        % to 360°).
        % Useage: az=azimuth(lat1,lon1,lat2,lon2, planet)
        %            azimuth(x-y_vals(row_i), x-y_vals(row_j+1))
%         az(k,1)=maptool_azimuth(y_vals(i,1),x_vals(i,1),y_vals(j+1,1),x_vals(j+1,1),S);
        az(k,1)=azimuth(y_vals(i,1),x_vals(i,1),y_vals(j+1,1),x_vals(j+1,1),S);
%         fprintf('x_vals(j+1,1) is %.2f\n',x_vals(j+1,1)) % display az value for debugging
%         fprintf('azimuth is %.2f\n',az(k,1)) % display az value for debugging
        %
%         % Conversion if out values should be limited from 0-360 to 0-180
%         if az(k,1) >= 180 % check for values between 180-360
%             az(k,1) = az(k,1) - 180; % subtract 180 to covert to 0-180
%         end % end IF statement
        %
        % (1) convert 180°-270° to 0°-90°
        if az(k,1) >= 180 && az(k,1) < 270
            az(k,1) = az(k,1) - 180;
        end % end IF statement
        %
        % (2) convert 270-360° to (0° to -90°)
        if az(k,1) >= 270 && az(k,1) < 360
            az(k,1) = az(k,1) - 360;
        end % end IF statement
        %
        % (3) convert 90-180° to (0° to -90°)
        if az(k,1) >= 90 && az(k,1) < 180
            az(k,1) = az(k,1) - 180;
        end % end IF statement
        %
		j = j + 1; % increment 2nd xy pair index
		k = k + 1; % increment angle counter index

	end % end inner FOR loop

	i = i + 1; % increment 1st xy pair index

end % end outer FOR loop


% Debugging: Display number of azimuths calculated for debugging
num_az = size(az,1); % debugging: Sum up the number of values in histogram
fprintf('The number of azimuths calculated in this run is %.0f\n\n', num_az);

