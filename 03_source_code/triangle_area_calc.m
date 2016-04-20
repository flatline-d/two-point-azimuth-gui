function [area] = triangle_area_calc(P)
%
% Function name = triangle_area_calc.m
% Author = Bradley Thomson (bjt@bu.edu)
% Created on 3 Dec 2012
% Last tinkered with on 3 Dec 2012
%
% The goal of this function is to import the a set of three x,y points and
% determine the volume of the triangular area contained within them.
%
% 	Function inputs:
%       -P : an 3 row x 2 col array with three x,y point sets
%
% 	Output:
%       -area : cartesian area between traid of input points
%
%
% Notes: Algorithm used is stable alternative to Heron's formula:
% http://en.wikipedia.org/wiki/Heron%27s_formula
% http://http.cs.berkeley.edu/~wkahan/Triangle.pdf


% Step 1: Determine a,b,c = lengths of three sides of a triangle.
% Sides a,b,c = L(1), L(2), L(3), order subject to sorting in next step.
% P is a 3x2 matrix of 3 x,y points with one x,y pair per line.
%
L=[sqrt(sum((P(1,:)-P(2,:)).^2)) sqrt(sum((P(2,:)-P(3,:)).^2)) sqrt(sum((P(3,:)-P(1,:)).^2))];


% Step 2: Sort distances so that a >= b >= c
%
L=sort(L); % Replace L with sorted version of L


% Step 3: Calculate area using stabilized Heron's formula 
area = sqrt( (L(3)+(L(2)+L(1)))*(L(1)-(L(3)-L(2)))*(L(1)+(L(3)-L(2)))*(L(3)+(L(2)-L(1))) )/4; 







