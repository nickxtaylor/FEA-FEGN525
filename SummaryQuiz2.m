%% solve an FE model comprising arbitrarily oriented truss elements

% clear all variables from memory, close all figures, clear command win
clear all; close all; clc;

%%

r = 2; % radius of round truss member in (mm)
A = pi * r^2; % cross-sectional area of bar in (mm2)
Le = 1000; % truss element length length in (mm)
E = 115e3; % modulus of titanium in (MPa)
P = 1000; % force on end of cantilever beam in (N)
kc = 1e12; % define kc for use in enforcing constrained DOF

%%

% define the 2D transformation matrices for the truss elements
beta1 = 45; c = cosd(beta1); s = sind(beta1);
T1 = [c s 0 0; 0 0 c s];
beta2 = -45; c = cosd(beta2); s = sind(beta2);
T2 = [c s 0 0; 0 0 c s];

%%

% find truss element stiffness value k, same for all three
k = (E * A) / Le;
% form the local stiffness matrix for the each truss element, same for all
kk = [k -k;-k k];
% transform the truss element stiffness matrices from components along local
% element DOF (along element axis in this case) to components aligned with global
% DOF (along x and y axes), the element stiffness matries will be 4x4 after
% this transformation
k1 = T1' * kk * T1;
k2 = T2' * kk * T2;
% global stiffness matrix is 8x8 since we have three elements, four nodes,
% and a total of 8 global DOF (u and v translations at each node)
K = zeros(6, 6);
% for computational simplicity, we will not expand the k's to global size, we
% will simply place each k into K in the locations that are relevant for each
% element; for the first element, k1 occupies the upper left 4x4 region of
% K, operating on global DOF 1,2,3,4
K(1:4,1:4) = K(1:4,1:4) + k1;
% for the second element, k2 occupies a central 4x4 region of K, operating on
% global DOF 3,4,5,6
K(3:6,3:6) = K(3:6,3:6) + k2;
% initialize Kc and then incorporate kc for fixed DOF: u1 = v1 = u3 = v3 = u4 = v4 = 0
Kc = K;
Kc(1, 1) = Kc(1, 1) + kc;
Kc(2, 2) = Kc(2, 2) + kc;
Kc(5, 5) = Kc(5, 5) + kc;
Kc(6, 6) = Kc(6, 6) + kc;

%%

% define column vector Fc of nodal loads, include P in the u direction at node 2
Fc = zeros(6, 1);
Fc(3) = P*cosd(160);
Fc(4) = P*sind(160);

% solve for unknown DOF using Kc and solve for nodal forces using K
D = inv(Kc) * Fc;
F = K*D;

%%

% % compute the stress in each element using (1.20), which is based on the
% % axial deformation in each element; axial deformation can be found by expressing
% % the displacement of node 2 as a vector and taking the dot product with a
% % unit vector along each element, note that for linear analysis we compute
% % stress based on the nominal configuration of the undeformed element
% d2 = [D(3) D(4)]; % displacment vector of node 2
% n1 = T1(1,1:2); % unit vector along element 1 from node 1 to node 2
% s1 = E * dot(d2,n1) / Le;
% n2 = -T2(1,1:2); % unit vector along element 2 from node 3 to node 2
% s2 = E * dot(d2,n2) / Le;
% n3 = -T3(1,1:2); % unit vector along element 3 from node 4 to node 2
% s3 = E * dot(d2,n3) / Le;
% stress = [s1 s2 s3]
% % we can also simply use stress = F/A, note we use the global numbering to
% % pull out components of reaction forces at each ground node, remember
% % these are simple two-force members so the reaction force is automatically
% % aligned with the axis of the element, norm() finds the magnitude of a
% % vector and always returns a positive number, so you have to reason the
% % sign manually
% s1 = norm([F(1) F(2)])/A;
% s2 = -norm([F(5) F(6)])/A;
% s3 = norm([F(7) F(8)])/A;
% stress = [s1 s2 s3]





