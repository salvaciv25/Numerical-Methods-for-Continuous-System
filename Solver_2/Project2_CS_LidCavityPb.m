clear;

%in this case c=0
bu = 0;
c = 0;
mu = 1;
tau = 1e-4;

%boundary functions
u1_b = @(x,y) H(x,y);
u2_b = @(x,y) 0;
p_b  = @(x,y) 0; 

f1_b = @(x,y) 0; 
f2_b = @(x,y) 0; 

for i = 1:5
    
    j = i - 1;
    folder1 = 'mesh';
    folder2 = sprintf('mesh%i', j);
    mesh = load(fullfile(folder1, folder2, 'mesh.dat'));
    xy = load(fullfile(folder1, folder2, 'xy.dat'));

    x = xy(:,1);
    y = xy(:,2);

    %compute boundary

    dirtot = 0;
    l = 1;
    for k = 1:size(xy,1)
        if x(k) == 0 || x(k) == 1 || y(k) == 0 || y(k) == 1
            dirtot(l) = k;
            l = l + 1;
        end
    end

    dirtot = dirtot';
    
    [u1, u2, p, res, x, y] = Stokes_solver(u1_b, u2_b, p_b, f1_b, f2_b, mesh, xy, dirtot, mu, c, tau, bu);

    figure(2 * i)

    subplot(2,2,1)
    trisurf(delaunay(x,y), x, y, u1, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('u1 GLS mesh%i', j));

    subplot(2,2,2)
    trisurf(delaunay(x,y), x, y, u2, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('u2 GLS mesh%i', j));

    subplot(2,2,3)
    trisurf(delaunay(x,y), x, y, p, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('p GLS mesh%i', j));   

    [u1, u2, p, res, x, y] = Stokes_solver(u1_b, u2_b, p_b, f1_b, f2_b, xy, mesh, dirtot, mu, c, 0, 1);

    figure(2 * i + 1)

    subplot(2,2,1)
    trisurf(delaunay(x,y), x, y, u1, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('u1 bubble mesh%i', j));

    subplot(2,2,2)
    trisurf(delaunay(x,y), x, y, u2, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('u2 bubble mesh%i', j));

    subplot(2,2,3)
    trisurf(delaunay(x,y), x, y, p, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('p bubble mesh%i', j));   
end