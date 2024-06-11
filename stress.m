function [ip_stress,node_stress] = stress(Emat,e)
% This function calculates stresses at integration points and nodes for a
% 6-node isoparametric triangle
% Inputs:
% Emat = the 3x3 elasticity matrix
% e = a 3x3 array of strain values, columns are integration points 1, 2, 3
%     rows are strain components x-normal, y-normal, xy-shear in that order
% Outputs:
% ip_stress = 3x3 array of stress values, columns are integration points 1, 2, 3
%             rows are stress components x-normal, y-normal, xy-shear in that order
% node_stress = 3x6 array of stress values, columns are nodes in the
%               element order and rows are stress components x-normal, y-normal,
%               xy-shear in that order

ip_stress = zeros(3,3);
% for each integration point
for i=1:3
    ip_stress(:,i) = Emat*e(:,i);
end

node_stress = zeros(3,6);
% for each stress component
for j=1:3
    % fit a plane to the integration point stress values
    v1 = [1/6 1/6 ip_stress(j,1)];
    v2 = [2/3 1/6 ip_stress(j,2)];
    v3 = [1/6 2/3 ip_stress(j,3)];
    vc = cross((v2-v1),(v3-v1)); vcn = vc/norm(vc);
    d = mean([dot(v1,vcn) dot(v2,vcn) dot(v3,vcn)]);
    node_stress(j,1) = (d - vcn(1)*(0) - vcn(2)*(0))/vcn(3);
    node_stress(j,2) = (d - vcn(1)*(1) - vcn(2)*(0))/vcn(3);
    node_stress(j,3) = (d - vcn(1)*(0) - vcn(2)*(1))/vcn(3);
    node_stress(j,4) = (d - vcn(1)*(0.5) - vcn(2)*(0))/vcn(3);
    node_stress(j,5) = (d - vcn(1)*(0.5) - vcn(2)*(0.5))/vcn(3);
    node_stress(j,6) = (d - vcn(1)*(0) - vcn(2)*(0.5))/vcn(3);
end
