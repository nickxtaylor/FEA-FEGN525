function s_ave = averaging(node,ele,snode)
% this function calculates average nodal stresses for all nodes in a mesh when more than one
% element is connected at any given node
% Inputs:
% node = an nx3 array of node coords where the first column is the node ID
%        and colums 2-3 are x, y coordinates, n is the total number of nodes
%        in the mesh
% ele = an mx7 array of element definitions where the first column is the
%       element ID and columns 2-7 are the node ID's of the six nodes used to
%       form the element, m is the total number of ele in the mesh
% snode = 3x6 array of stress values, columns are nodes in the
%         element order and rows are stress components x-normal, y-normal,
%         xy-shear in that order
% Outputs:
% s_ave = px1 x n 3D array of average stress values where...
%         p = the number of stress components calculated at each node, this
%         number is 3 by default, but the user can put addional values in
%         additional rows and the code below will account for this
%         automatically
%         n = the number of "pages" in this 3D array, this is also the number
%         of nodes in the mesh so there is a page for each node, each page
%         is a px1 vector of average stress values at that node

% for each node in the mesh
for i=1:length(node(:,1))
    % find which ele this node is part of
    [r,c] = find(ele(:,2:7)==node(i,1));
    % initialize the average stress based on the number of stress
    % components in snode
    s_ave(:,1,i) = zeros(size(snode(:,1,1)));
    
    % for each ele identified in r (this assumes ID's are sequential)
    for j=1:length(r)
        % add all the stress values for that node
        s_ave(:,1,i) = s_ave(:,1,i) + snode(:,c(j),r(j));
    end
    % average based on the number of ele connected at that node
    s_ave(:,1,i) = s_ave(:,1,i)/length(r);
end

