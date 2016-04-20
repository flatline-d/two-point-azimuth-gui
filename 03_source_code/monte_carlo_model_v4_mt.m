function [bn,nh,th] = monte_carlo_model_v4_mt(xy_pts_file_name, observed_hist, axes7_handle, MC_runs, planet)
%
% Goal: The point of this function is to run a series of Monte Carlo
% models, tabulate the results, and output a normalized azimuth histogram.
% Note this version confine points to the convex polygon defined by the
% outermost points in the distribution (i.e., the convex hull)
%
% Author: Brad Thomson  (bjt@bu.edu)
% Last tinkered with: 26 Aug 2012
%
% 	Function inputs:
%       -xy_pts_file_name : file where each x,y point pair is the latitude
%       and longitude of each shield field. (.txt or .csv)
%       -observed_hist : histogram of observed data
%       -axes7_handle : handle of axes to plot normalized histogram
%       -MC_runs : number of Monte Carlo runs to perform
%       -planet : selected planetary body
%
% 	Output:
%       -Histogram of normalized azimuth values
%       -bn = bin_centers = center values of histogram, direct function output
%       -nh = norm_freq =  normalized histogram values, direct function output
%       -th = thresh = 95% confidence threshold, direct function output
%
%   Sub-functions called:
%       -read_xy_point_data : read in designated input file and return
%       arrays of X and Y values
%       -two_pt_azimuth_calcs_v3_mt : calculate 2-point azimuth values,
%       uses mapping toolbox
%       -triangle_area_calc : calculate area enclosed by 3 x,y points
%
% Change list:
%
% 2013-08-23: Added b,n,t as direction output of function (bins, normalized
% histogram values, and 95% threshold values).
%
% 2013-08-26: Added planet as input
%
% 4 Aug 2015: Changed to subfunction "two_pt_azimuth_calcs_v3_mt"
% instead of "two_pt_azimuth_calcs_v3."


% Initialize variables, counters,  and results matricies 
v = MC_runs; % number of Monte Carlo runs, taken from function input
i = 1; % initialize counter for number of Monte Carlo runs
k = 18; % k is the number of bins (there are 18 10° bins between +/- 90°
bin_edges = -90:10:90; % create 18 histogram bins from -90 to 90 (10° bins)
bin_centers = -85:10:85; % 10 deg bins, n = 18 between -90 to 90 degrees
MC_out = zeros(k,v); % intialize output: # bins (rows) vs. # runs (cols)



%% Read in data file
%

% Open x,y data file for reading using subfuction "read_xy_point_data.m"
%
[x_vals, y_vals] = read_xy_point_data(xy_pts_file_name);

% Count number of points (i.e., shield volcanoes) in input file
%
n_shields = size(x_vals,1);



%% Calculate bounding polygon using convex hull method (NOT bounding box!)
%

% (1) Find convex hull
%
% K = convhull(X,Y) returns the 2-D convex hull of the points (X,Y), where
% X and Y are column vectors. The convex hull K is a vector of point 
% indices arranged in a counterclockwise cycle around the hull.
% [K,V] = convhull(...) returns the convex hull K and the corresponding
% area/volume V bounded by K.
[K, V] = convhull(x_vals, y_vals);

% (2) Divide convex hull into triangles using Delaunay triangulation
%
% DT = DelaunayTri(x,y) creates a Delaunay triangulation from a set of 
% points. The triangulation divides a polygon into a set of constituent
% triangles. Note that only the indicies identified by convhull are used.
hullpts = [x_vals(K), y_vals(K)];
DT = DelaunayTri(hullpts);

% (3) Compute the area/volume of each triangle using subfunction
% "triangle_area_calc"
%
area = zeros(size(DT,1),1); % initialize area data holder array
for j = 1:size(DT,1) % size of 1st dimension of DT is number of triangles
    area(j) = triangle_area_calc([ hullpts(DT(j,1),:); hullpts(DT(j,2),:); hullpts(DT(j,3),:) ]);
end % end FOR loop

% % (4) Assign each triangle a weighing factor proportional to its area.
% % From http://www.mathworks.com/matlabcentral/newsreader/view_thread/284645
% %
% c = cumsum([0 area']); % create a new array c that is the cumulative sum of each value of area
% c = c/c(end); % divide each value by the last value (i.e., total area)
% c(end) = inf; % replace last value with positive infinity
% % Count number of points randomly assigned to each triangle in proportion
% % to its (fractional) area
% b = histc(rand(n_shields,1),c)
% b(end) = []; % shorten result b by one cell to match number of triangles



%% Monte Carlo model
%
% (5) In each triangle, generate the number of random points determined in
% step (4)
%

for i = 1:v % v = number of Monte Carlo runs

fprintf('Monte Carlo run number %.0f of %.0f total runs.\n', i, v)

% (4) Assign each triangle a weighing factor proportional to its area.
% From http://www.mathworks.com/matlabcentral/newsreader/view_thread/284645
%
c = cumsum([0 area']); % create a new array c that is the cumulative sum of each value of area
c = c/c(end); % divide each value by the last value (i.e., total area)
c(end) = inf; % replace last value with positive infinity
% Count number of points randomly assigned to each triangle in proportion
% to its (fractional) area
b = histc(rand(n_shields,1),c);
b(end) = []; % shorten result b by one cell to match number of triangles

p = 1; % initialize counter for number of random points
rand_pts = zeros(n_shields,2);
for j = 1:size(DT,1) % size of 1st dimension of DT is number of triangles
    if b(j) ~= 0 % check if there are any points in given triangle
    % Note elements of b contain the number of points in each triangle    
        for m = 1:b(j)

            % If we have two random numbers a0 and a1 (between 0-1),
            % we calculate the random point x as:
            %   x = a0 (v1 - v0) + a1 (v2 - v0)
            % where v1, v2, v3 are the vertices of the triangle.
            v0 = hullpts(DT(j,1),:);
            v1 = hullpts(DT(j,2),:);
            v2 = hullpts(DT(j,3),:);
            rand_pts(p,:) = v0 + ( rand(1)*(v1-v0) + rand(1)*(v2-v0) );

            % Check if point is in triangle using MATLAB function inpolygon
            % Syntax: IN = inpolygon(X,Y,xv,yv)
            % X,Y are points to test
            % xv,yv are verticies of enclosed polygonal region (a triangle in this case)
            xv = [v0(1); v1(1); v2(1)];
            yv = [v0(2); v1(2); v2(2)];
            % IN is 1 if point is inside triangle, 0 otherwise
            IN = inpolygon(rand_pts(p,1),rand_pts(p,2),xv,yv);
            if IN == 0
                % rotate point 180 deg around line v1-v2
                % From http://stackoverflow.com/questions/240778/random-points-inside-a-polygon
                % comment "answered Oct 27 '08 at 18:09"
                rand_pts(p,:) = v0 + (v1-v0) + (v2-v0) - (rand_pts(p,:)-v0); 
            end % end IF statement
            
            p = p+1; % increment counter for number of random points

        end % end inner FOR loop
    end % end IF statement
end % end outer FOR loop


r_lons = rand_pts(:,1);
r_lats = rand_pts(:,2);
% figure; % debugging: open new figure for visual check
% triplot(DT); % debugging: plot Delaunay triangulation 
% hold on % debugging
% plot(r_lons,r_lats,'^'); % debugging: plot scatter plot for visual check
% xlabel('Lon (°E)') % debugging
% ylabel('Lat (°N)')% debugging
% axis image  % debugging: forces aspect ratio to be one-to-one
% hold off % debugging

% Tabulate results of Monte Carlo model run
%
% Compute 2-point azimuth values using subfunction "two_pt_azimuth_calcs_v3_mt"
[az] = two_pt_azimuth_calcs_v3_mt(r_lons, r_lats, planet);


% Compile azimuth values into histogram
n = histc(az,bin_edges); % n is array with number of values in each bin
                         % note last bin counts values equal to bin edge
% Assign results to output array, removing last bin
MC_out(:,i) = n(1:(size(n)-1));

%     figure; % debugging: open new figure for visual check
%     bar(n2); % debugging: plot bar plot for visual check

end % end FOR loop



%% Gather statistics of Monte Carlo results
%

% Calculate mean and standard deviation
%
means = mean(MC_out,2); % mean(A,2) is a column vector containing the mean
                        % value of each row
sigma = std(MC_out,0,2); % calculate standard deviation of each row

% Calculate expected frequency per cell
%   expect_freq = n*(n-1)/2k
%   n*(n-1)/2 is the total number of azimuths
%   k is the number of bins
%
exp_freq = (n_shields*(n_shields-1))/(2*k); % expected frequency per cell

% Determine normalized frequency
%   norm = (expected / Monte Carlo) * observed
norm_freq = (exp_freq ./ means) .* observed_hist;
% norm_factor = (exp_freq ./ means) % debugging: report normalization factor

% Determine normaled 95th% threshold value using the Student's t distribution
%	threshold value = (expected / Monte Carlo) * (Monte Carlo + [sigma * tinv(0.95, v)])
%	v equals N-1 number of runs
%   "tinv" is the Student's t inverse cumulative distribution function
thresh = (exp_freq ./ means) .* ( means+(sigma*tinv(0.95,v)) );



%% Plot results
%
axes(axes7_handle) % designate active axes
%
% Plot azimuth values in 2D histogram. Recall N is 0 deg increasing
% clockwise (90° = east)
bar(bin_centers, norm_freq, 'hist')
%
set(gca, 'XTick',-90:30:90, 'XLim',[-90 90]) % set X-axis limit to -90 to +90 with 30° ticks
xlabel('Azimuth (°E)'); % label x-axis of plot
ylabel('Normalized Freq.'); % label y-axis of plot
% Get y-axis limits
y_range = ylim; % ylim with no arguments returns the limits of the current axes

% Overlay normalized threshold values on bar plot
hold on
plot(bin_centers,thresh,'*');
set(gca, 'YLim', y_range); % match y-axis limits of bar plot
hold off



%% Assign data to function output
%
bn = bin_centers; % center values of histogram
nh = norm_freq; % normalized histogram values
th = thresh; % 95 percent confidence threshold