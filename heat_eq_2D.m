clear
clc
Code="Y";
while Code=="Y" || Code=="y"
    
    Nx = 50; % Number of grids in x, y-direction
    Ny = 50;
    Lx = 10; % Sheet size in x, y-direction
    Ly = 10;
    
    x = linspace(0,Lx,Nx); %Compose grid in x, y-direction
    y = linspace(0,Ly,Ny);
    x = (x(1:end-1)+x(2:end))/2;
    y = (y(1:end-1)+y(2:end))/2;
    [X,Y] = meshgrid(x,y);
    
    dx = x(2)-x(1); % Since the sheet is uniform, dx = x_n - x_n-1
    dy = y(2)-y(1); % Since the sheet is uniform, dy = y_n - y_n-1
    
    u_0(:,:) = getINT(X,Y);
    
    alpha = menu('Choose the type of material used/wanted:','Aluminum',...
        'Copper','Gold','Iron', 'Custom'); %alpha for specific materials, animation is accelerated by 10E7 for convenience
    if alpha == 5
        alpha = input("Enter thermal conductivity value: [range: 1-5] --> ");
    else
        choice = load('material_type.mat','-ascii');
        alpha  = choice(alpha);
    end
    
    dt = (dx^2+dy^2)/alpha;
    
    graph_surf = surf(X,Y,u_0); % 3-D Plot
    graph_axes = gca;
    graph_axes.ZLim = [-7,7];
    graph_axes.CLim = [-7,7];
    colorbar;
    title('2-D Heat equation');
    
    tspan = linspace(0,1,1/dt+70); % Modeling
    [t,u] = ode15s(@(t,x)getRHS(x,alpha,dx,Nx,Ny),tspan,u_0(:));
    Tn = length(t);
    u = reshape(u,Tn,Nx-1,Ny-1);
    
    filename = 'heat_eq.gif';
    for i=1:Tn
        Z = u(i,:,:);
        Z = squeeze(Z);
        graph_surf.ZData = Z;
        drawnow;
        frame = getframe(gcf);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
        if i==1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.7);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.7);
        end
    end
    Code=input('Do you want to run again? (Y/N): ', 's');
end

function u_0 = getINT(X,Y)
% Used to formulate the Initial Coditions
%
%   Input: X, Y
%   Output: u_0
%
%See also input.

Type = input('Please choose: 1) Gaussian Distributions 2)Sine and Cosine waves --> ');
while Type~=1 && Type~=2
    Type = input('Error!! Please choose a valid type number: 1 or 2: ');
end
if Type==1
    u_0(:,:) = peaks; % Uses the Gaussian distributions due to its 3D useful demonstration
else
    sin_A = input('Please enter sine wave amplitude: ');
    cos_A = input('Please enter cosine wave amplitude: ');
    u_0(:,:) = cos_A*cos(X)+sin_A*sin(Y);
end
end

function du = getRHS(u,alpha,dx,Nx,Ny)
% Copyright 2015-2016 The MathWorks, Inc.

% Reshape the date in 2D
u = reshape(u,Nx-1,Ny-1);

% Copy it to new array with boundaries.
u_bounded = zeros(Nx+1,Ny+1);
u_bounded(2:end-1,2:end-1) = u;

% Neumann Boundary condition
% set the zero-gradient = no heat flux
u_bounded(1,:) = u_bounded(2,:);
u_bounded(end,:) = u_bounded(end-1,:);
u_bounded(:,1) = u_bounded(:,2);
u_bounded(:,end) = u_bounded(:,end-1);

% Get the second derivatives
u = u_bounded;
du = alpha/dx^2*(u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1)...
    + u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end));

du = du(:);
end
