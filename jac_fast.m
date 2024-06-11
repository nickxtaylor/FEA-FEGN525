function [J,Jinv] = jac_fast(n,ip)
% This function computes the jacobian for a CST and a LST finite element
% Inputs:
% n = array of node coordinates, each row is a node, columns are x, y
%     for CST n must be 3x2, for LST it must be 6x2
% ip = vector of itegration point coordinates, it's 1x2, columns are r, s
% Outputs:
% J = determinant of the Jacobian matrix
% Jinv = inverse of the Jacobian matrix

len = length(n(:,1));
x = n(:,1); y = n(:,2);

switch len
    case 3
        error('myToolbox:myFunction:wrongNumberNodes','Code for 3-node TRI has not been written yet.');
    case 6
        % shape functions for a 6-node TRI
        J11 = @(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,r,s) x2*(4*conj(r) - 1) + x1*(4*conj(r) + 4*conj(s) - 3) - x4*(8*conj(r) + 4*conj(s) - 4) + 4*x5*conj(s) - 4*x6*conj(s);
        J12 = @(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,r,s) y2*(4*conj(r) - 1) + y1*(4*conj(r) + 4*conj(s) - 3) - y4*(8*conj(r) + 4*conj(s) - 4) + 4*y5*conj(s) - 4*y6*conj(s);
        J21 = @(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,r,s) x3*(4*conj(s) - 1) + x1*(4*conj(r) + 4*conj(s) - 3) - x6*(4*conj(r) + 8*conj(s) - 4) - 4*x4*conj(r) + 4*x5*conj(r);
        J22 = @(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,r,s) y3*(4*conj(s) - 1) + y1*(4*conj(r) + 4*conj(s) - 3) - y6*(4*conj(r) + 8*conj(s) - 4) - 4*y4*conj(r) + 4*y5*conj(r);
        j11 = J11(n(1,1),n(2,1),n(3,1),n(4,1),n(5,1),n(6,1),n(1,2),n(2,2),n(3,2),n(4,2),n(5,2),n(6,2),ip(1),ip(2));
        j21 = J21(n(1,1),n(2,1),n(3,1),n(4,1),n(5,1),n(6,1),n(1,2),n(2,2),n(3,2),n(4,2),n(5,2),n(6,2),ip(1),ip(2));
        j12 = J12(n(1,1),n(2,1),n(3,1),n(4,1),n(5,1),n(6,1),n(1,2),n(2,2),n(3,2),n(4,2),n(5,2),n(6,2),ip(1),ip(2));
        j22 = J22(n(1,1),n(2,1),n(3,1),n(4,1),n(5,1),n(6,1),n(1,2),n(2,2),n(3,2),n(4,2),n(5,2),n(6,2),ip(1),ip(2));
        % the Jacobian matrix for a 6-node TRI
        Jmat = [j11 j12;j21 j22];
        J = det(Jmat);
        % inverse of Jmat, used for computing B
        Jinv = inv(Jmat);
    otherwise
        error('myToolbox:myFunction:wrongNumberNodes','Wrong number of nodes... must be 3 or 6.');
end
