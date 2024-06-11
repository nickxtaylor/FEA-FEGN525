function [J,Jinv] = jac(n,ip)
% This function computes the jacobian for a CST and a LST finite element
% Inputs:
% n = array of node coordinates, each row is a node, columns are x, y
%     for CST n must be 3x2, for LST it must be 6x2
% ip = vector of itegration point coordinates, it's 1x2, columns are r, s
% Outputs:
% J = determinant of the Jacobian matrix
% Jinv = inverse of the Jacobian matrix

len = length(n(:,1));
syms r s;
x = sym('x',[len 1]); y = sym('y',[len 1]);
cx = n(:,1); cy = n(:,2);

switch len
    case 3
        % shape functions for a 3-node TRI
        N = [1-r-s r s];
        % the Jacobian matrix for a 3-node TRI, see equation (6.12)
        Jmat = jacobian(N,[r s])'*[x y];
        % substitute the node coordinates and integration point location
        Jmat = double(subs(Jmat,[r s x(1) x(2) x(3) y(1) y(2) y(3)],[ip(1) ip(2) cx(1) cx(2) cx(3) cy(1) cy(2) cy(3)]));
        J = det(Jmat);
        % inverse of Jmat, used for computing B
        Jinv = inv(Jmat);
    case 6
        % shape functions for a 6-node TRI
        N = [(1-r-s)*(1-2*r-2*s), r*(2*r-1), s*(2*s-1), 4*r*(1-r-s), 4*r*s, 4*s*(1-r-s)];
        % the Jacobian matrix for a 3-node TRI, see equation (6.12)
        Jmat = jacobian(N,[r s])'*[x y];
        % substitute the node coordinates and integration point locations
        Jmat = subs(Jmat,[x(1) x(2) x(3) x(4) x(5) x(6)],[cx(1) cx(2) cx(3) cx(4) cx(5) cx(6)]);
        Jmat = subs(Jmat,[y(1) y(2) y(3) y(4) y(5) y(6)],[cy(1) cy(2) cy(3) cy(4) cy(5) cy(6)]);
        Jmat = double(subs(Jmat,[r s],[ip(1) ip(2)]));
        J = det(Jmat);
        % inverse of Jmat, used for computing B
        Jinv = inv(Jmat);
    otherwise
        error('myToolbox:myFunction:wrongNumberNodes','Wrong number of nodes... must be 3 or 6.');
end
