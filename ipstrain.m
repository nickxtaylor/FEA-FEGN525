function strain = ipstrain(B,D,gDOF)
% This function calculates strain at the integration points for a
% 6-node isoparametric triangle
% Inputs:
% B = a 3x12 x 3 3D array of B matrices evaluated at the integration
%     points, the third index is the page... i.e., page 1 is a 3x12 B array
%     evaluated at integration point 1, pages 2 and 3 are similarly
%     defined
% D = the solution for global DOF for the model, size is 2*n x 1 where n is
%     the total number of nodes in the model, entries are x1, y1, x2, y2,
%     and so on through xn, yn
% gDOF = a 1x12 vector of nodal DOF associated with this element, order is
%        the same as the order of node ID's used to form the element
% Outputs:
% strain = a 3x3 array of integration point strain values, each element has three integration
%          points which are the columns of the 3x3 strain array, rows are
%          x-normal, y-normal, xy-shear strains

for i=1:3
    strain(:,1) = B(:,:,1)*D(gDOF);
    strain(:,2) = B(:,:,2)*D(gDOF);
    strain(:,3) = B(:,:,3)*D(gDOF);
end