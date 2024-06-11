function [node,ele] = read_mesh_TRI6(filename)
% This function reads node and element definitions from any Abaqus input
% file, it is assumed (a) that all nodes are defined under a single *Node
% keyword and all ele are defined under a single *Element keyword, and (b)
% that the elements of interest are 6-node triangles
% Inputs:
% filename = a character string defining the name of the file to open, the
%            three-letter extension should be included
% Outputs:
% node = an nx3 array of node coords where the first column is the node ID
%        and colums 2-3 are x, y coordinates, n is the total number of nodes read
%        from the file
% ele = an mx7 array of element definitions where the first column is the
%       element ID and columns 2-7 are the node ID's of the six nodes used to
%       form the element, m is the total number of ele read from the file

% open the file containing the triangle mesh
fid_mesh = fopen(filename);
% arrays for nodes and elements
node = []; ele = [];

% parse node and ele definitions for triangle mesh
ln = fgetl(fid_mesh);
while(~feof(fid_mesh))
   if (length(ln)>=5 & strcmp(ln(1:5),'*Node')),
       % read nodes
       ln=fgetl(fid_mesh);
       while(~strcmp(ln(1),'*')),
           node = [node;cell2mat(textscan(ln,'%f64,%f64,%f64'))];
           ln = fgetl(fid_mesh);
       end
   elseif (length(ln)>=8 & strcmp(ln(1:8),'*Element')),
       % read ele
       ln=fgetl(fid_mesh);
       while(~strcmp(ln(1),'*')),
           ele = [ele;cell2mat(textscan(ln,'%d32,%d32,%d32,%d32,%d32,%d32,%d32'))];
           ln = fgetl(fid_mesh);
       end
   else
       ln=fgetl(fid_mesh);
   end 
end

% close the file containing the mesh
fclose(fid_mesh);