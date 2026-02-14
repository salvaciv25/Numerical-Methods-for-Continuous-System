function [SYSMAT,RHS]=lifting(SYSMAT,RHS,DirNod,DirVal)
% SYSMAT: matrix of the linear system
% RHS: right hand side: vector of zeros
% dirichletNodes: indices of the Dirichlet nodes
% uD: values of the Dirichlet boundary condition
    
for i = 1:length(DirNod)

    dnode = DirNod(i);
    SYSMAT(dnode,:) = 0;
    SYSMAT(dnode,dnode) = 1;

    RHS(dnode) = DirVal(i);
    

end
end