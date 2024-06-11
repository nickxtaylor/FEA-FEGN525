function K = efixed(K,nodeID,eDOF)
% this function applies kc = 1e12 to any local DOF that are to be fixed to
% zero displacement, K is the global stiffness matrix, nodeID is a vector of node ID's
% that will be affected, and eDOF is a vector of the local DOF to be fixed,
% the updated K is returned with the kc values implemented
% 
% For example, if one wants to fix local DOF 2 (y-translation) at nodes 7
% and 9, then the function call would look something like this...
% Kc = efixed(K,[7 9],[2])

for i=1:length(nodeID)
    for j=1:length(eDOF)
        n = nodeID(i)*2-2 + eDOF(j);
        K(n,n) = K(n,n) + 1e12;
    end
end