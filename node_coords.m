function [gDOF,coords] = node_coords(node,ele,n)
% This function finds the global DOF and x,y coordinates for each node used
% to form a specific element designated by element ID "n"
% Inputs:
% node = an px3 array of node coords where the first column is the node ID
%        and colums 2-3 are x, y coordinates (mm), p is the total number of
%        nodes in the model
% ele = an mx7 array of element definitions where the first column is the
%       element ID and columns 2-7 are the node ID's of the six nodes used to
%       form the element, m is the total number of ele in the model
% n = the ID of the element of interest
% Outputs:
% gDOF = a 1x12 vector of the global DOF associated with all element DOF in the same order
%        as the node ID's used to form the element
% coords = a 6x2 array of x,y coordinates for all of the nodes used to form
%          this element, rows are in the same order as the node ID's, units
%          are (mm)

[r,c] = size(ele);
gDOF = zeros(1,2*(c-1));
for i=2:c
    j = ele(n,i);
    k = find(node(:,1)==j);
    coords(i-1,:) = node(k,2:end);
    gDOF(2*(i-1)-1:2*(i-1)) = [2*ele(n,i)-1 2*ele(n,i)];
end