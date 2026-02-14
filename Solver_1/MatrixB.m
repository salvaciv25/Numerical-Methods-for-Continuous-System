function [stiff_mat, transpmat, stab_mat]=MatrixB(N_elem,Nodes,triang,b,c,Area,beta,Diff,tau)

h = 0;
stiff_mat = sparse(441,441);
transpmat = sparse(Nodes,Nodes);
stab_mat = sparse(Nodes,Nodes);

%norm of K
K = sum(Diff)/2;

%norm of Beta
B = sqrt((sum(beta.^2)));

if B == 0
        B = 1;
end

%% Compute general matrices H, B and S
for iel = 1:N_elem
    
    %length of the triangle
    for s = 1:3
            hnew = sqrt((b(iel,s)*Area(iel)*2)^2 + (c(iel,s)*Area(iel)*2)^2);
            if hnew > h
               h = hnew;
            end
    end
    
    % Building local matrices and adding them to the general ones
    for iloc = 1:3
        iglob = triang(iel,iloc);

        for jloc = 1:3
            jglob = triang(iel,jloc);
            stiff_mat(iglob,jglob) = stiff_mat(iglob,jglob) + (b(iel,iloc)*b(iel,jloc)*Diff(1) + c(iel,iloc)*c(iel,jloc)*Diff(2))*(Area(iel));
            transpmat(iglob,jglob) = transpmat(iglob,jglob) + (b(iel,jloc)*beta(1) + c(iel,jloc)*beta(2))*Area(iel)*(1/3);
            stab_mat(iglob,jglob) = stab_mat(iglob,jglob) + ((tau*h*Area(iel))/(4*B*K))*(b(iel,iloc)*beta(1) + c(iel,iloc)*beta(2))*(b(iel,jloc)*beta(1) + c(iel,jloc)*beta(2));
        end
    end
end


end