function [u1, u2, p, res, x, y, u_ex] = matrix_build(u1_b, u2_b, p_b, f1_b, f2_b, xy, mesh, dirtot, mu, c, tau, bu)
%%% mes = N_elem
%%% dirnod = bou
%%% mesb = N_elem_b
%%% i = iel
%%% coo = Nodes
%%% l = mesh(i, 1) = triang(iel,1); m = mesh(i, 2)= triang(iel,2); n = mesh(i, 3)= triang(iel,3);
% tau parameter for the GLS stabilized
% c constant parameter that multiplies the stiffness matrix
% mu viscocity coefficient
% bu indicator for the bubble method
%% continuare da linea 106
    if bu ~= 0
        bu = 1;
    end

%mesh's data
    N_elem = size(mesh, 1);
    N_dir= size(dirtot, 1); 

   

    N_elem_b= bu * N_elem;      %mesb

    Nodes = size(xy, 1) + N_elem_b; %coo
   
    bnod = 3 + bu; 

    x = xy(:, 1);
    y = xy(:, 2);

    % baricenters nodes
    if bu == 1
        for iel = 1:N_elem

            mesh(iel,4)= Nodes - N_elem + iel;

            x(Nodes - N_elem + iel) = (x(mesh(iel,1)) + x(mesh(iel,2)) + x(mesh(iel,3))) / 3;
            y(Nodes - N_elem + iel) = (y(mesh(iel,1)) + y(mesh(iel,2)) + y(mesh(iel,3))) / 3;
        end
    end

    %B.C., rhs and exact solution
    [u1, u2, p, f1, f2] = deal(zeros(Nodes, 1));
    
    for j = 1:Nodes
        u1(j) = u1_b(x(j),y(j));
        u2(j) = u2_b(x(j),y(j));
        p(j) = p_b(x(j),y(j));
        f1(j) = f1_b(x(j),y(j));
        f2(j) = f2_b(x(j),y(j));
    end
    u_ex = [u1(1:Nodes - N_elem_b); u2(1:Nodes - N_elem_b); p(1:Nodes - N_elem_b)];
    
    [u1_rhs, u2_rhs, p_rhs] = deal(zeros(Nodes, 1));
   
    
    % initialize matrices 
    [K, M, B1, B2, C] = deal(sparse(Nodes, Nodes));

    A = sparse(3*Nodes,3*Nodes);
    b_loc=zeros(N_elem,3);
    c_loc=zeros(N_elem,3);
    
    [u1_rhs, u2_rhs, p_rhs]= deal(zeros(Nodes,1));
    bu * N_elem

    %% Compute the general matrices
    for iel = 1:N_elem

        % useful data
        p1=xy(mesh(iel,1),:);
        p2=xy(mesh(iel,2),:);
        p3=xy(mesh(iel,3),:);
        
        Area = 0.5 * abs(det([1 p1(1) p1(2); 1 p2(1) p2(2); 1 p3(1) p3(2)]));

        nodes=mesh(iel,1:3);
     
       for inod=1:3

           n1 = mod_n(inod+1,3);
           n2 = mod_n(inod+2,3);

           %a_loc(iel,inod)= x(nodes(n1))*y(nodes(n2))...    
           % - x(nodes(n2))*xy(nodes(n1),1); da controllare
           b_loc(iel,inod)=(y(nodes(n1)) - y(nodes(n2))) / (2 * Area);            
           c_loc(iel,inod)=(x(nodes(n2)) - x(nodes(n1))) / (2 * Area);
       end
       
       %important coefficient: devi calcolarli a mano
        I00 = Area;
        I10 = 2 * Area / 6;
        I20 = 2 * Area / 24;
        I30 = 2 * Area / 120;
        I11 = 2 * Area / 60;
        I21 = 2 * Area / 360;
        I02 = 2 * Area / 180;
        I12 = 2 * Area / 1260;
        I03 = 2 * Area / 5040;
        
        %local matrices and rhs
        f = [f1(mesh(iel,1)) + f1(mesh(iel,2)) + f1(mesh(iel,3)); f2(mesh(iel,1)) + f2(mesh(iel,2)) + f2(mesh(iel,3))] / 3;
        
        [ K_l, C_l, B1_l, B2_l, M_l] = deal(zeros(bnod, bnod));
        
        M_l(1:3, 1:3) = Area * (eye(3 , 3) + ones(3 , 3)) / 12;
        
        [u1_rhs_loc, u2_rhs_loc, p_rhs_loc] = deal(zeros(bnod, 1));


        
        %% P1/P1 GLS
        for iloc = 1:3
            for jloc = 1:3
                K_l(iloc,jloc) = (b_loc(iel,iloc)*b_loc(iel,jloc) + c_loc(iel,iloc)*c_loc(iel,jloc))*I00;
                C_l(iloc,jloc) = - tau * (b_loc(iel,iloc)*b_loc(iel,jloc) + c_loc(iel,iloc)*c_loc(iel,jloc))*I00;
            end
            B1_l(1:3,iloc) = b_loc(iel,:)' * I10;
            B2_l(1:3,iloc) = c_loc(iel,:)' * I10;

            u1_rhs_loc(iloc) = f(1) * I10;
            u2_rhs_loc(iloc) = f(2) * I10;
            p_rhs_loc(iloc) = (- tau) * (f(1) * b_loc(iel,iloc) + f(2) * c_loc(iel,iloc)) * I00;

        end

        %% Bubble
        if bu == 1
            for iloc = 1:3
                for jloc = 1:3
                    K_l(iloc,4) = K_l(iloc,4) + 27 * (b_loc(iel,iloc)*b_loc(iel,jloc) + c_loc(iel,iloc)*c_loc(iel,jloc)) * I20;

                    if jloc ==  iloc

                        B1_l(iloc,4)= B1_l(iloc,4) + 27 * b_loc(iel,jloc) *I30;
                        B2_l(iloc,4)= B2_l(iloc,4) + 27 * c_loc(iel,jloc) *I30;
                    else
                        B1_l(iloc,4)= B1_l(iloc,4) + 27 * b_loc(iel,jloc) *I11;
                        B2_l(iloc,4)= B2_l(iloc,4) + 27 * c_loc(iel,jloc) *I11;
                    end
                end

                    K_l (4,iloc) = K_l(iloc,4);

                    M_l(iloc,4) = 27 * I21;
                    M_l(4,iloc) = M_l(iloc,4);
             end

                u1_rhs_loc(4) = 27 * f(1) * I30;
                u2_rhs_loc(4) = 27 * f(2) * I30;

                B1_l(4,:) = 27 * [b_loc(iel,1), b_loc(iel,2), b_loc(iel,3), 0] * I30;
                B2_l(4,:) = 27 * [c_loc(iel,1), c_loc(iel,2), c_loc(iel,3), 0] * I30;

                for iloc = 1:3
                    for jloc = iloc:3
                        if jloc == iloc
                            K_l(4,4) = K_l(4,4) + (b_loc(iel,iloc)*b_loc(iel,jloc) + c_loc(iel,iloc)*c_loc(iel,jloc)) * (27^2) * I02;
                        else
                            K_l(4,4) = K_l(4,4) + (b_loc(iel,iloc)*b_loc(iel,jloc) + c_loc(iel,iloc)*c_loc(iel,jloc)) * (27^2) * I21;
                        end
                    end

                    B1_l(4,4) = B1_l(4,4) +  b_loc(iel,iloc) * (27^2) * I12;
                    B2_l(4,4) = B2_l(4,4) +  c_loc(iel,iloc) * (27^2) * I12;

                end
                M_l(4,4) = (27^2)*I03;
         end

                   
 %%% computing global matrices and rhs
        for iloc = 1:bnod
              iglob = mesh(iel,iloc);
           
              for jloc = 1:bnod
              jglob = mesh(iel,jloc);
            
              K(iglob,jglob) = K(iglob, jglob) + K_l(iloc,jloc);
              C(iglob,jglob) = C(iglob,jglob) + C_l(iloc,jloc);
              B1(iglob,jglob) = B1(iglob,jglob) + B1_l(iloc,jloc);
              B2(iglob,jglob) = B2(iglob,jglob) +B2_l(iloc,jloc);
              M(iglob,jglob) = M(iglob, jglob) + M_l(iloc,jloc);
              end

           
            u1_rhs(iglob) = u1_rhs(iglob) + u1_rhs_loc(iloc);
            u2_rhs(iglob) = u2_rhs(iglob) + u2_rhs_loc(iloc);
            p_rhs(iglob) = p_rhs(iglob) + p_rhs_loc(iloc);        
        end
    end

    A1 = c * M + mu * K;
    A2 = c * M + mu * K;

    T1 = B1';
    T2 = B2';

    %boundary conditions
    for i = 1:N_dir
        j = dirtot(i);

        A1(j, :) = 0;
        A2(j, :) = 0;
        T1(j, :) = 0;
        T2(j, :) = 0;
        B1(j, :) = 0;
        B2(j, :) = 0;
        C(j, :) = 0;

        A1(j, j) = 1;
        A2(j, j) = 1;
        C(j, j) = 1;

        u1_rhs(j) = u1(j);
        u2_rhs(j) = u2(j);
        p_rhs(j) = p(j);
    end
        
    %general matrix
    A = [A1, zeros(Nodes, Nodes), T1; zeros(Nodes, Nodes), A2, T2; B1, B2, C];    
    A = A(1:3 * Nodes - N_elem_b, 1:3 * Nodes - N_elem_b); %eliminate p-baricenters

    rhs = [u1_rhs; u2_rhs; p_rhs(1:Nodes - N_elem_b)]; %rhs without p-baricenters

    %% solution
    u = A\rhs;

    % take out the baricenters from A,rhs and the solution
    Nodes = Nodes - N_elem_b;

    A = [A1(1:Nodes, 1:Nodes), zeros(Nodes, Nodes), T1(1:Nodes, 1:Nodes); zeros(Nodes, Nodes), A2(1:Nodes, 1:Nodes), T2(1:Nodes, 1:Nodes); B1(1:Nodes, 1:Nodes), B2(1:Nodes, 1:Nodes), C(1:Nodes, 1:Nodes)];

    rhs = [u1_rhs(1:Nodes); u2_rhs(1:Nodes); p_rhs(1:Nodes)];

    u1 = u(1 : Nodes);
    u2 = u(Nodes + N_elem_b + 1 : 2 * (Nodes + N_elem_b) - N_elem_b);
    p = u(2 * (Nodes + N_elem_b) + 1: 3 * (Nodes + N_elem_b) - N_elem_b);    

    x = x(1:Nodes);
    y = y(1:Nodes);

    %residual error
    res = norm(A * u_ex - rhs);
end