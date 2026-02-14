clear

% initial data
c= 1;
mu = 1;
tau = 0.1;

u1 = @(x,y) -cos(2*pi*x)*sin(2*pi*y)+sin(2*pi*y);
u2 = @(x,y) sin(2*pi*x)*cos(2*pi*y)-sin(2*pi*x);
p  = @(x,y) 2*pi*(cos(2*pi*y)-cos(2*pi*x));
f1 = @(x,y) c*u1(x,y)-4*mu*pi^2*sin(2*pi*y)*(2*cos(2*pi*x)-1)+4*pi^2*sin(2*pi*x);
f2 = @(x,y) c*u2(x,y)+4*mu*pi^2*sin(2*pi*x)*(2*cos(2*pi*y)-1)-4*pi^2*sin(2*pi*y);

%Residual's Matrix
res= zeros(3,5);
res_u1 = zeros(3,5);
 res_u2 = zeros(3,5);
 res_p = zeros(3,5);

for i = 1:5

    %input variables
    j = i-1;
    folder1 = 'mesh';
    folder2 = sprintf('mesh%i', j);
    dirtot = load(fullfile(folder1, folder2, 'dirtot.dat'));
    mesh = load(fullfile(folder1, folder2, 'mesh.dat'));
    xy = load(fullfile(folder1, folder2, 'xy.dat'));

    %% P1/P1 
    tau = 0;
    bu = 0;
    [u1_p1, u2_p1, p_p1, res(1,i),res_u1(1,i), res_u2(1,i), res_p(1,i), x, y, u_ex] = Stokes_solver(u1, u2, p, f1, f2, xy, mesh, dirtot, mu, c, tau, bu);

    %% P1/P1 stabilized
    tau = 0.1;
    bu = 0;
    [u1_st, u2_st, p_st, res(2,i), res_u1(2,i), res_u2(2,i), res_p(2,i)] = Stokes_solver(u1, u2, p, f1, f2, xy, mesh, dirtot, mu, c, tau, bu);

    %% P1-Bubble/P1
    tau = 0;
    bu = 1;
    [u1_bu, u2_bu, p_bu, res(3,i), res_u1(3,i), res_u2(3,i), res_p(3,i)] = Stokes_solver(u1, u2, p, f1, f2, xy, mesh, dirtot, mu, c, tau, bu);

    u1_ex = u_ex(1:size(xy, 1));
    u2_ex = u_ex(size(xy, 1) + 1:2 * size(xy, 1));
    p_ex = u_ex(2 * size(xy, 1) + 1: 3 * size(xy, 1));

    figure(3 * i + 1)
    subplot(2,2,1)
    trisurf(delaunay(x,y), x, y, u1_p1, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u1 of P1/P1');

    subplot(2,2,2)
    trisurf(delaunay(x,y), x, y, u1_st, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u1 of GLS')

    subplot(2,2,3)
    trisurf(delaunay(x,y), x, y, u1_bu, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u1 of P1-bubble/P1');


    subplot(2,2,4)
    trisurf(delaunay(x,y), x, y, u1_ex, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u1 exact');

    figure(3*i + 2)
    subplot(2,2,1)
    trisurf(delaunay(x,y), x, y, u2_p1, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u2 of P1/P1');

    subplot(2,2,2)
    trisurf(delaunay(x,y), x, y, u2_st, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u2 of GLS')

    subplot(2,2,3)
    trisurf(delaunay(x,y), x, y, u2_bu, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u2 of P1-bubble/P1');


    subplot(2,2,4)
    trisurf(delaunay(x,y), x, y, u2_ex, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('u2 exact');

    figure(3*i + 3)
    subplot(2,2,1)
    trisurf(delaunay(x,y), x, y, p_p1, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('p of P1/P1');

    subplot(2,2,2)
    trisurf(delaunay(x,y), x, y, p_st, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('p of GLS')

    subplot(2,2,3)
    trisurf(delaunay(x,y), x, y, p_bu, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('p of P1-bubble/P1');


    subplot(2,2,4)
    trisurf(delaunay(x,y), x, y, p_ex, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('p exact');
end

T = array2table(res, 'VariableNames', {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4'});
T.Properties.RowNames = {'P1/P1', 'P1/P1 GLS', 'P1-bubble/P1'};

disp(T);

T1 = array2table(res_u1, 'VariableNames', {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4'});
T1.Properties.RowNames = {'P1/P1', 'P1/P1 GLS', 'P1-bubble/P1'};

disp(T1);

T2 = array2table(res_u2, 'VariableNames', {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4'});
T2.Properties.RowNames = {'P1/P1', 'P1/P1 GLS', 'P1-bubble/P1'};

disp(T2);

T3 = array2table(res_p, 'VariableNames', {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4'});
T3.Properties.RowNames = {'P1/P1', 'P1/P1 GLS', 'P1-bubble/P1'};

disp(T3);