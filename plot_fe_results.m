function [] = plot_fe_results(node,ele,D,snode_avg,scale)
% make basic contour plots of the displacement and von Mises stress results
% for easy comparison to Abaqus
% Inputs:
% node = an px3 array of node coords where the first column is the node ID
%        and colums 2-3 are x, y coordinates (mm), p is the total number of
%        nodes in the model
% ele = an mx7 array of element definitions where the first column is the
%       element ID and columns 2-7 are the node ID's of the six nodes used to
%       form the element, m is the total number of ele in the model
% D = the solution for global DOF for the model, size is 2*n x 1 where n is
%     the total number of nodes in the model, entries are x1, y1, x2, y2,
%     and so on through xn, yn
% snode_avg = px1 x n 3D array of average stress values where...
%             p = the number of stress components calculated at each node, this
%             number is 3 by default, but the user can put addional values in
%             additional rows, n = the number of "pages" in this 3D array, this is also the number
%             of nodes in the mesh so there is a page for each node, each page
%             is a px1 vector of average stress values at that node
% scale = a simple scale factor to magnify the appearanc of deformation in
%         the plot so it's visible... see line 27 below

nnod = length(node(:,1));
nele = length(ele(:,1));

% make an (x,y) array of deformed nodal coordinates
Dnode = reshape(D,2,nnod)'; Dcoord = node(:,2:3) + scale*Dnode;
% compute displacement magnitude at each node
for i=1:length(Dnode),Dmag(i,1) = norm(Dnode(i,:));end
% make a normalized vector of displacement magnitudes for the contour plot
% and find the node location where displacement magnitude is max
Dnorm = normc(Dmag); [~,dindex] = max(Dmag);

% create a vector of just the averaged nodal vM stresss
vMstress = reshape(snode_avg(5,:,:),nnod,1);
% find the max von Mises stress and its node location
[sMax,vmindex] = max(vMstress);
% make a normalized vector of stress magnitudes for the contour plot
vMnorm = normc(vMstress);

% plot nominal and deformed shape with contours of displacement magnitude
subplot(1,2,1); triplot(double(ele(:,2:4)),node(:,2),node(:,3),'k--');
hold on; axis equal; axis off; %axis([-120 90 -220 240]);
title(['Displacement Magnitude, max = ',num2str(max(Dmag)),' (mm) at node ',num2str(dindex)]);
p.Vertices = Dcoord; p.Faces = ele(:,2:4);
p.FaceColor = 'interp'; p.FaceVertexCData = Dnorm;
patch(p);
% mark the max location
plot(Dcoord(dindex,1),Dcoord(dindex,2),'g.','MarkerSize',25);

% plot nominal and deformed shape with contours of vM stress
subplot(1,2,2); triplot(double(ele(:,2:4)),node(:,2),node(:,3),'k--');
hold on; axis equal; axis off; %axis([-120 90 -220 240]);
title(['von Mises Stress, max = ',num2str(sMax),' (MPa) at node ',num2str(vmindex)]);
p.Vertices = Dcoord; p.Faces = ele(:,2:4);
p.FaceColor = 'interp'; p.FaceVertexCData = vMnorm;
patch(p);
% mark the max location
plot(Dcoord(vmindex,1),Dcoord(vmindex,2),'g.','MarkerSize',25);
