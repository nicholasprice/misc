function [W,S] = ellipseIntermediateAxis(A,B,th,beta)
% for an ellipse with properties:
% - A = MajorAxisLength
% - B = MinorAxisLength
% - th= orientation in radians
% This function returns 'W', the radial distance along angle beta
%
% We do this by solving for the intersection between the equations for:
% 1) a rotated ellipse with orientation th relative to the x-axis
% 2) a straight line with angle beta relative to the x-axis
% the equation for a rotated ellipse is here
% https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
% This uses the symbolic toolbox
% 
% OUTPUTS
% W = radial distance from the ellipses centre to the intersection along
% line with angle beta
% S = position of the intersection
%
% Note that there can be problems evaluating things if the angles are 0 or
% pi. It may help to add a small error, e.g. use (pi/2 + 0.01) as an angle
%
% Example usage:
% - obtain A,B and th from regionprops
% - obtain beta based on angle between two ellipse centers
%
% NP 13 Jul 2023

if nargin==0
%     % redundant case - a circle
%     A = 1; B = 1;
%     th = 0; beta = 0;

% an ellipse aligned with the x-axis
A = 2; B = 1;
th = 0;
beta = pi/4;
% beta = pi/2 + 0.01; 

% a complicated example
%     th = 45*pi/180;
%     beta = 60*pi/180;
%     A = 2;
%     B = 1;
end

syms x y
A = sym(A);
B = sym(B);
th = sym(th);
beta = sym(beta);

% our two equations. The second one is the ellipse formula
eqns = [tan(beta) == y/x, ((x*cos(th)+y*sin(th))^2/A^2) + ((x*sin(th)-y*cos(th))^2/B^2) == 1];

S = solve(eqns,[x y],'PrincipalValue',true); % return the first solution to the two equations with the specified inputs

W = eval(sqrt(S.x^2 + S.y^2));