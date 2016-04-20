readme_source_code_file_list.txt

File description: This text files contains a detailed list of the source file components used in the Graphical User Inferface (GUI) gui_shield_field_azimuth_v1p1_mt.m. All files listed are .m files, which are MATLAB program files.

Author : Brad Thomson (bjt@bu.edu), Boston University Center for Remote Sensing
Date : 2016-04-19

Main program, i.e., GUI for two-point azimuth method:
	gui_shield_field_azimuth_v1p2_mt.m
%
% Notes: The “_mt” appended to the file name indicates that the MATLAB mapping toolbox was used by some functions. So if the user decides to modify the code an re-run, the Mapping Toolbox will have to be enabled.

Functions called:
(1) plot_xy_point_data.m
% Notes: called upon execution of “Plot shield locations” button, plots point data in left most panel of GUI.

(2a) azimuth_method_v3_mt.m
% Notes: Called if radio button “Lutz method” is selected upon selection of the “Compute Azimuth” button in middle panel. This function uses the Lutz two-point azimuth method and outputs a histogram and rose diagram in the middle panel of the GUI. 
%   Sub-functions called:
%       (2a.1) read_xy_point_data : read in designated input file and return arrays of X and Y values.
%       (2a.2) two_pt_azimuth_calcs_v3_mt : calculate 2-point azimuth values.
%     Sub-Sub-functions called:
%           (2a.2.1) referenceEllipsoid : From Mapping Toolbox, used to define a planetary ellipsoid
%		with specific radii.
%           (2a.2.2) azimuth : From Mapping Toolbox, used to find angle and distance between two
%           points on a given reference ellipsoid.

(2b) azimuth_method_cebria_mt.m
% Notes: Called if radio button “Cebria method” is selected upon selection of the “Compute Azimuth” button in middle panel. This function uses the Cebria two-point azimuth method and outputs a histogram and rose diagram in the middle panel of the GUI. 
%   Sub-functions called:
%       (2b.1) read_xy_point_data : read in designated input file and return
%       arrays of X and Y values.
%       (2b.2) two_pt_azimuth_cebria_v3_mt : calculate 2-point azimuth values.
%     Sub-sub-functions called:
%           (2b.2.1) referenceEllipsoid : From Mapping Toolbox, used to define a planetary ellipsoid
%		 with specific radii.
%           (2b.2.2) azimuth : From Mapping Toolbox, used to find angle and distance between two 
%           points on a given reference ellipsoid.
%       (2b.3) plot_xy_point_data.m : replots x-y data in left hand panel, but adds line segments
%       between points that meet Cebria’s criteria.

(3a) monte_carlo_model_v4_mt.m
% Notes: Called if radio button “Lutz method” is selected upon selection of the “Run Model” button in right-hand panel. This function runs a user-specific number of Monte Carlo models, computes the Lutz method on each, and then averages the results. The output is a normalized histogram in the right-hand panel of the GUI. 
%   Sub-functions called:
%       (3a.1) read_xy_point_data : read in designated input file and return arrays of X and Y values.
%       (3a.2) two_pt_azimuth_calcs_v3_mt : calculate 2-point azimuth values.
%     Sub-Sub-functions called:
%           (3a.2.1) referenceEllipsoid : From Mapping Toolbox, used to define a planetary ellipsoid
%		with specific radii.
%           (3a.2.2) azimuth : From Mapping Toolbox, used to find angle and distance between two 
%           points on a given reference ellipsoid.
%       (3a.3) triangle_area_calc : calculate area enclosed by 3 x,y points.

(3b) monte_carlo_model_v4_cebria_mt.m
% Notes: Called if radio button “Cebria method” is selected upon selection of the “Run Model” button in right-hand panel. This function runs a user-specific number of Monte Carlo models, computes the Cebria method on each, and then averages the results. The output is a normalized histogram in the right-hand panel of the GUI. 
%     Sub-functions called:
%         (3b.1) read_xy_point_data.m : read in desigated input file and return
%         arrays of X and Y values
%         (3b.2) two_pt_azimuth_cebria_MC_only_mt.m : calculate 2-point azimuth values using
%		Cebria et al. (2011) method. Similar to two_pt_azimuth_calcs_v3_mt, but here the number
%		of azimuth values in each empirical Monte Carlo model matches that in the observed
%		distribution.
%     Sub-sub-functions called:
%           (3b.2.1) referenceEllipsoid : From Mapping Toolbox, used to define a planetary ellipsoid
%		with specific radii.
%           (3b.2.2) azimuth : From Mapping Toolbox, used to find angle and distance between two 
%           points on a given reference ellipsoid.
%       (3b.3) triangle_area_calc.m : calculate area enclosed by 3 x,y points.
