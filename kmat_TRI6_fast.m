function [k,Bmat] = kmat_TRI6_fast(coords,t,Emat)
% This function computes the stiffness matrix for a 6-node isoparametric
% triangle element (v2)
% Inputs:
% coords = a 6x2 array of x,y coordinates for the 6 nodes used to form
%         this element, rows are in the same order as the node ID's, units
%         are (mm)
% t = uniform thickness of the element in (mm)
% Emat = the 3x3 elasticity matrix
% Ouputs:
% k = the 12x12 stiffness matrix for this element, units are (N/mm)
% Bmat = a 3x12 x 3 3D array of B matrices evaluated at the integration
%        points, the third index is the page... i.e., page 1 is a 3x12 B array
%        evaluated at integration point 1, pages 2 and 3 are similarly
%        defined

% integration points for a 6-node TRI, coords in r, s
ip = [1/6 1/6;2/3 1/6;1/6 2/3];
% shape functions for 6-node TRI
% N = [(1-r-s)*(1-2*r-2*s), r*(2*r-1), s*(2*s-1), 4*r*(1-r-s), 4*r*s, 4*s*(1-r-s)];
% derivatives for later finding B
% Bx = diff(N,'r'); By = diff(N,'s');
% Bx = [ 4*r + 4*s - 3, 4*r - 1, 0, 4 - 4*s - 8*r, 4*s, -4*s];
% By = [ 4*r + 4*s - 3, 0, 4*s - 1, -4*r, 4*r, 4 - 8*s - 4*r];
% Bx = { @(r,s) 4*r + 4*s - 3, @(r,s) 4*r - 1, @(r,s) 0, @(r,s) 4 - 4*s - 8*r, @(r,s) 4*s, @(r,s) -4*s };
% By = { @(r,s) 4*r + 4*s - 3, @(r,s) 0, @(r,s) 4*s - 1, @(r,s) -4*r, @(r,s) 4*r, @(r,s) 4 - 8*s - 4*r };
% Bfun = @(r,s) [Bx{1}(r,s)   0   Bx{2}(r,s)   0   Bx{3}(r,s)   0   Bx{4}(r,s)   0   Bx{5}(r,s)   0   Bx{6}(r,s)   0  ;...
%                By{1}(r,s)   0   By{2}(r,s)   0   By{3}(r,s)   0   By{4}(r,s)   0   By{5}(r,s)   0   By{6}(r,s)   0  ;...
%                 0   Bx{1}(r,s)   0   Bx{2}(r,s)   0   Bx{3}(r,s)   0   Bx{4}(r,s)   0   Bx{5}(r,s)   0   Bx{6}(r,s) ;...
%                 0   By{1}(r,s)   0   By{2}(r,s)   0   By{3}(r,s)   0   By{4}(r,s)   0   By{5}(r,s)   0   By{6}(r,s)];
Bfun1 = [...
   -1.666666666666667         0   -0.333333333333333         0         0         0    2.0000         0    0.666666666666667         0   -0.666666666666667         0;...
   -1.666666666666667         0         0         0   -0.333333333333333         0   -0.666666666666667         0    0.666666666666667         0    2.0000         0;...
         0   -1.666666666666667         0   -0.333333333333333         0         0         0    2.0000         0    0.666666666666667         0   -0.666666666666667;...
         0   -1.666666666666667         0         0         0   -0.333333333333333         0   -0.666666666666667         0    0.666666666666667         0    2.0000];
Bfun2 = [...
    0.333333333333333         0    1.666666666666667         0         0         0   -2.0000         0    0.666666666666667         0   -0.666666666666667         0;...
    0.333333333333333         0         0         0   -0.333333333333333         0   -2.666666666666667         0    2.666666666666667         0    0.0000         0;...
         0    0.333333333333333         0    1.666666666666667         0         0         0   -2.0000         0    0.666666666666667         0   -0.666666666666667;...
         0    0.333333333333333         0         0         0   -0.333333333333333         0   -2.666666666666667         0    2.666666666666667         0    0.0000];
Bfun3 = [...
    0.333333333333333         0   -0.333333333333333         0         0         0    0.0000         0    2.666666666666667         0   -2.666666666666667         0;...
    0.333333333333333         0         0         0    1.666666666666667         0   -0.666666666666667         0    0.666666666666667         0   -2.0000         0;...
         0    0.333333333333333         0   -0.333333333333333         0         0         0    0.0000         0    2.666666666666667         0   -2.666666666666667;...
         0    0.333333333333333         0         0         0    1.666666666666667         0   -0.666666666666667         0    0.666666666666667         0   -2.0000];
% Jacobian calculations
[J,Jinv] = jac_fast(coords,ip(1,:));
% B matrix for element 1
B = [1 0 0 0;0 0 0 1;0 1 1 0]*[Jinv zeros(2,2);zeros(2,2) Jinv];
% B11 = B*Bfun(ip(1,1),ip(1,2));
B11 = B*Bfun1;
% substitute the integration point location if needed
% B11 = double(subs(B,[r s],ip(1,:)));

% save this result as page 1 in the 3D Bmat
Bmat(:,:,1) = B11;
J11 = (1/3)*0.5*J; % J was computed above on line 30
k11 = t * J11 * (B11'*Emat*B11);
%
% Jacobian calculations
[J,Jinv] = jac_fast(coords,ip(2,:));
% B matrix for element 1
B = [1 0 0 0;0 0 0 1;0 1 1 0]*[Jinv zeros(2,2);zeros(2,2) Jinv];
% B12 = B*Bfun(ip(2,1),ip(2,2));
B12 = B*Bfun2;
% substitute the integration point location if needed
% B12 = double(subs(B,[r s],ip(2,:)));
% save this result as page 2 in the 3D Bmat
Bmat(:,:,2) = B12;
J12 = (1/3)*0.5*J; % J was computed above on line 43
k12 = t * J12 * (B12'*Emat*B12);
%
% Jacobian calculations
[J,Jinv] = jac_fast(coords,ip(3,:));
% B matrix for element 1
B = [1 0 0 0;0 0 0 1;0 1 1 0]*[Jinv zeros(2,2);zeros(2,2) Jinv];
% B13 = B*Bfun(ip(3,1),ip(3,2));
B13 = B*Bfun3;
% substitute the integration point location if needed
% B13 = double(subs(B,[r s],ip(3,:)));
% save this result as page 3 in the 3D Bmat
Bmat(:,:,3) = B13;
J13 = (1/3)*0.5*J; % J was computed above on line 56
k13 = t * J13 * (B13'*Emat*B13);
%
% compute k with Gauss quadrature... three points for LST
k = k11 + k12 + k13;

