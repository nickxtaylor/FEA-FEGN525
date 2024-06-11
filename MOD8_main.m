%% FEGN 525 MOD 8 SUMMARY REPORT Main Script

% Nick Taylor
% 10920730

% hosuekeeping
clear all; clc; close all;

%% define constants

% constants defining model parameters
t = 1; % thickness in (mm)
nu = 0.35; % Poisson's ratio
E = 20e3; % modulus of bone (MPa)
c = E/((1+nu)*(1-2*nu)); % material constants
Emat = [(1-nu)*c nu*c 0;nu*c (1-nu)*c 0;0 0 (1-2*nu)*c/2]; % E matrix, plane-strain


%% build geometry and initiate hooke's law variables accordingly

% read the mesh from an external file
[node,ele] = read_mesh_TRI6('femur_2D_TRI6.inp');
nnod = length(node(:,1));
nele = length(ele(:,1));

% initialize K and Fc based on 2 DOF per node
K = zeros(2*nnod,2*nnod);
Fc = zeros(2*nnod,1);

%% build global stiffness matrix

% initiate B cell
B = cell(nele);

% compute [k] for each element and assemble [K]
for i=1:nele

    % define node locations for element
    [gDOF,coords] = node_coords(node,ele,i);

    % get the stiffness and other data for the element
    %[k,N{i},B{i}] = kmat_TRI6(coords,t,Emat);
    [k,B{i}] = kmat_TRI6_fast(coords,t,Emat);

    % assemble this element into global stiffness
    K = assem(K,k,gDOF);

end

%% adjust stiffness according to constraints and compute hookes law

% use the efixed() function to assign fixed constraints to local DOF, local
% DOF are numbered as u = 1, v = 2, w = 3 at any given node
Kc = efixed(K,[68:112],[1 2]); % beam mesh0

% apply loads to global DOF
Fc(184*2-1) = -20; Fc(184*2) = -50; % beam mesh0

% solve for nodal DOF with Kc and solve for reactions with K
D = inv(Kc)*Fc;
F = K*D;

%% compute strain and stress

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
    
end

%% plot

% compute average nodal stresses
snode_avg = averaging(node,ele,snode);

% plot deformed results for displacements and von Mises stress
plot_fe_results(node,ele,D,snode_avg,20)

%%

%%

%%

%%

%%

%%

%%