% solve an FE model comprising TRI6 elements

% clear all variables from memory, close all figures, clear command win
clear all;close all;clc;

% constants defining model parameters
t = 1; % thickness in (mm)
nu = 0.34; % Poisson's ratio
E = 115e3; % modulus of titanium in (MPa)
c = E/((1+nu)*(1-2*nu)); % material constants
Emat = [(1-nu)*c nu*c 0;nu*c (1-nu)*c 0;0 0 (1-2*nu)*c/2]; % E matrix, plane-strain

% read the mesh from an external file
[node,ele] = read_mesh_TRI6('LST_beam_mesh0.inp');
nnod = length(node(:,1));
nele = length(ele(:,1));
% initialize K and Fc based on 2 DOF per node
K = zeros(2*nnod,2*nnod);
Fc = zeros(2*nnod,1);

% make a cell vector to hold all the N vectors, this is just a 1D cell vector of size equal
% to the number of ele in the mesh, each item in the cell vector is a
% symbolic 1x6 array of shape functions where N = N(r,s)
%N = cell(nele);
% make a cell array to hold all the B arrays, this is just a 1D cell vector of size equal
% to the number of ele in the mesh, each item in the cell vector is a 3D B array with size
% pxq x r, where pxq is 3x12 and is a B array for a single integration point
% the third index r refers to pages in the 3D B array corresponding to the
% integration points... r = 1, 2, 3 to designate the integration point
% B = cell(nele);

% compute [k] for each element and assemble [K]
for i=1:nele
    % define node locations for element
    [gDOF,coords] = node_coords(node,ele,i);
    % get the stiffness and other data for the element
    %[k,N{i},B{i}] = kmat_TRI6(coords,t,Emat);
    [k,B{i}] = kmat_TRI6_fast(coords,t,Emat);
    % assemble this element into global stiffness
    K = assem(K,k,gDOF);
    % print a message to the screen
    clc;disp(['Assembled element ',num2str(i),' of ',num2str(nele),'...']);
end

% use the efixed() function to assign fixed constraints to local DOF, local
% DOF are numbered as u = 1, v = 2, w = 3 at any given node
Kc = efixed(K,[1 12],[1 2]); % beam mesh0
% apply loads to global DOF
Fc(22) = -1; Fc(44) = -1; % beam mesh0

% solve for nodal DOF with Kc and solve for reactions with K
D = inv(Kc)*Fc;
F = K*D;

% compute strain & stress for each element
for i=1:nele
    % compute strain at integration points, stress at itegration points and nodes
    % eip, sip, and snode are all 3D arrays... the third index is the
    % "page" and it refers to an element ID, each page contains ip or node
    % values in columns (3 col for ip, 6 col for nodes), rows are x-normal, y-normal, xy-shear
    eip(:,:,i) = ipstrain(B{i},D,node_coords(node,ele,i));
    [sip(1:3,:,i),snode(1:3,:,i)] = stress(Emat,eip(:,:,i));
    
    % put ip z-stress in row 4 of sip
    for j=1:3, sip(4,j,i) = nu*sip(1,j,i) + nu*sip(2,j,i);end
    % put ip von Mises stress in row 5 of sip
    for j=1:3, sip(5,j,i) = sqrt(sip(1,j,i)^2 + sip(2,j,i)^2 + sip(4,j,i)^2 - sip(1,j,i)*sip(2,j,i) - sip(2,j,i)*sip(4,j,i) - sip(4,j,i)*sip(1,j,i) + 3*sip(3,j,i)^2);end
    
    % put nodal z-stress in row 4 of snode
    for j=1:6, snode(4,j,i) = nu*snode(1,j,i) + nu*snode(2,j,i);end
    % put nodal von Mises stress in row 5 of snode
    for j=1:6, snode(5,j,i) = sqrt(snode(1,j,i)^2 + snode(2,j,i)^2 + snode(4,j,i)^2 - snode(1,j,i)*snode(2,j,i) - snode(2,j,i)*snode(4,j,i) - snode(4,j,i)*snode(1,j,i) + 3*snode(3,j,i)^2);end
    % print a message to the screen
    clc;disp(['Computing strain/stress for element ',num2str(i),' of ',num2str(nele),'...']);
end

% compute average nodal stresses
snode_avg = averaging(node,ele,snode);
% plot deformed results for displacements and von Mises stress
plot_fe_results(node,ele,D,snode_avg,20)

% print some results for comparison to Abaqus and theory
%
% peak values from FE model
vTip_FE = min(D)
sMax_FE = max(max(snode(1,:,:)))
% beam equations
vTip_exact = -2*12*(500)^3/3/E/20^3
sMax_exact = (2)*(500)*(10)*12/20^3



