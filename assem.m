function K = assem(K,k,gDOF)
% this function places entries from element stiffness k into the global
% stiffness K based on a vector of global DOF numbers (gDOF) corresponding
% to the DOF for the particular element, the updated K is returned
for i=1:length(gDOF)
    for j=1:length(gDOF)
        K(gDOF(i),gDOF(j)) = K(gDOF(i),gDOF(j)) + k(i,j);
    end
end