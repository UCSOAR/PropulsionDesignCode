% MOC to compute a minimum length bell nozzle contour
% Author: Mia Mrljic
% Date: 2024-11-02
% Outputs: Contour Plot, Characteristic Plot, 3D Plot, Nozzle STL

%% INPUTS
n = 20; %  number of characteristics (refinement parameter)
Ve = 800; % exit velocity, from CEA
Rt = 1; % throat radius for scaling
Ri = 2; % inlet radius
gamma = 1.4; % Cp/Cv, from CEA
Tamb = 200; % ambient temperature
MM = 28.96/1000; % molar mass, from CEA
scaling = 'ND'; % ND for non-dimensional, D for dimensional
ThetaInlet = 45; % inlet cone angle


%% DO NOT TOUCH BELOW THIS LINE

%% TEST
node=0.5*n*(4+n-1);

[V,Theta,mach,mu,ThetaMax,Me,gamma] = MOCminimum(n,Ve,gamma,Tamb,MM);
bellWall = plotCharacteristics(Theta,mu,ThetaMax,n,scaling,Rt,Me,gamma); 
[convergentWall,conicalWall] = convergentContour(Rt,Ri,ThetaInlet,scaling);

figure(2)
hold on
plot([conicalWall(:,1),conicalWall(:,1)], [conicalWall(:,2),-conicalWall(:,2)], color='b')
plot([convergentWall(:,1),convergentWall(:,1)], [convergentWall(:,2),-convergentWall(:,2)], color='m')
plot([bellWall(:,1),bellWall(:,1)], [bellWall(:,2),-bellWall(:,2)], color='r')
title(sprintf('%.0f Degree Inlet Bell Nozzle Contour for Mach=%.2f and Cp/Cv=%.2f',ThetaInlet,Me,gamma))
axis equal
hold off

z = [conicalWall(:,1); convergentWall(:,1); bellWall(:,1)];
r = [conicalWall(:,2); convergentWall(:,2); bellWall(:,2)];
a = linspace(0, 2*pi, 60);

[R,A] = ndgrid(r,a);
Z = repmat(z, 1, length(a));
[X,Y,Z] = pol2cart(A,R,Z);
figure(3)
mesh(X, Y, Z, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'k');  % Gray face, no edges
lighting gouraud
camlight headlight
axis equal
grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('%.0f Degree Inlet Bell Nozzle Revolved for Mach=%.2f and Cp/Cv=%.2f',ThetaInlet,Me,gamma));

% STL
% Convert the mesh data to a triangulation object for STL export
[F, V] = surf2patch(X, Y, Z, 'triangles');
TR = triangulation(F, V);

% Export to STL
stlwrite(TR, 'nozzle.stl');
disp('STL file "nozzle.stl" created successfully.');


%% FUNCTIONS
% convergent contour Rao nozzle
function [convergentWall,conicalWall] = convergentContour(Rt,Ri,ThetaInlet,scaling)
    if scaling == 'ND'
        Rt = 1;
        Ri = Ri/Rt;
    end

    % throat section
    Phi = (-90-ThetaInlet):1:-90;
    x = 1.5*Rt*cosd(Phi);
    y = 1.5*Rt*sind(Phi) + 2.5*Rt;

    % conical section
    Xt = (Ri-y(1,1))/tand(180-ThetaInlet)+x(1,1);

    % wall contour
    convergentWall = [x', y'];
    conicalWall = [Xt, Ri; x(1,1), y(1,1)];
end


% Prandtl Meyer function for ideal gas
function v = PrandtlMeyer(M,gamma)
    v=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atan(sqrt(M^2-1));
end

% inverse Prandtl Meyer function for ideal gas
function mach = InvPrandtlMeyer(V)
    V = V.*pi/180;
    Vo = pi*(sqrt(6)-1)/2;
    y = (V./Vo).^(2/3);
    one = ones(1,size(y,2));
    mach = (one+1.3604.*y+0.0962.*(y.^2)-0.5127.*(y.^3))./(one-0.6722.*y-0.3278.*(y.^2));
end

% minimum length bell nozzle
function [V,Theta,mach,mu,ThetaMax,Me,gamma] = MOCminimum(n,Ve,gamma,Tamb,MM)
    % total nodes
    nodes=0.5*n*(4+n-1);
    
    % speed of sound
    Vs = sqrt(gamma*8.314*Tamb/MM);
    
    % exit mach number
    Me = Ve/Vs;
    
    % exit Prandtl-Meyer
    Vme = PrandtlMeyer(Me,gamma);
    
    % max contour angle
    ThetaMax = Vme/2*180/pi;
    
    % subdivide domain
    dTheta = ThetaMax/n;
    
    % set characteristics
    Cplus = zeros(1,nodes);
    Cminus = zeros(1,nodes);
    Theta = zeros(1,nodes);
    V = zeros(1,nodes);
    
    % set nodes
    ifirst = 1:n+1;
    
    icenter = [n+2];
    value = n;
    while icenter(end) < nodes-1
        icenter = [icenter,icenter(end)+value];
        value = value - 1;
    end
    
    iwall = [n+1+n];
    value = n-1;
    while iwall < nodes
        iwall = [iwall,iwall(end)+value];
        value = value-1;
    end
    
    iint = [];
    for i=1:nodes
        if (ismember(i,ifirst) == 0) && (ismember(i,icenter) == 0) && (ismember(i,iwall) == 0)
            iint = [iint,i];
        end
    end
    
    charCounter = 0;
    for i = 1:nodes
        % first right running characteristics
        if ismember(i,ifirst) == 1
            Theta(i) = i.*dTheta;
            V(i) = Theta(i);
            Cminus(i) = V(i) + Theta(i);
            Cplus(i) = Theta(i) - V(i);
        end
        Theta(n+1) = Theta(n);
        V(n+1) = Theta(n+1);
        Cminus(n+1) = V(n+1) + Theta(n+1);
        Cplus(n+1) = Theta(n+1) - V(n+1);
        
        % centerline
        if ismember(i,icenter) == 1
            Theta(i) = 0;
            Cminus(i) = Cminus(i-(n-charCounter));
            V(i) = Cminus(i) - Theta(i);
            Cplus(i) = Theta(i) - V(i);
        end
        
        % intermediate
        if ismember(i,iint) == 1
            Cplus(i) = Cplus(i-1);
            Cminus(i) = Cminus(i-(n-charCounter));
            Theta(i) = (1/2)*(Cplus(i) + Cminus(i));
            V(i) = (1/2)*(Cminus(i) - Cplus(i));
        end
        
        % wall
        if ismember(i,iwall) == 1
            Theta(i) = Theta(i-1);
            V(i) = V(i-1);
            Cplus(i) = Cplus(i-1);
            Cminus(i) = Cminus(i-1);
            charCounter = charCounter + 1;
        end
    end
    
    % mach angle
    mach = InvPrandtlMeyer(V);
    mu = asind(1./mach);
end


function bellWall = plotCharacteristics(Theta,mu,ThetaMax,n,scaling,Rt,Me,gamma)    
    % calculate nodes
    nodes=0.5*n*(4+n-1);
    
    % contour coordinates non-dimensional
    % initialize scaling
    if scaling == 'ND'
        R = 1;
    else
        R = Rt;
    end
    
    % initialize coorindates
    x = zeros(1,nodes);
    y = zeros(1,nodes);
    figure(1)
    bellWall = [0,R];
    
    % first characteristic
    i = 1;
    while (i<=n+1)
        % throat node
        if i==1
            x(i)=-R/(tand(Theta(i)-mu(i)));
            y(i)=0;
            plot([0 x(i)],[R 0]);
            hold on
        else 
            % wall node
            if i==n+1
                x(i)=(y(i-1)-R-x(i-1)*tand((Theta(i-1)+Theta(i)+mu(i-1)+mu(i))*0.5))/(tand(0.5*(ThetaMax+Theta(i)))-tand((Theta(i-1)+Theta(i)+mu(i-1)+mu(i))*0.5));
                y(i)=R+x(i)*tand(0.5*(ThetaMax+Theta(i)));
                bellWall = [bellWall;x(i),y(i)];
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
                plot([0 x(i)],[R y(i)]);
                hold on
            % intermediate nodes
            else
                x(i)=(R-y(i-1)+x(i-1)*tand(0.5*(mu(i-1)+Theta(i-1)+mu(i)+Theta(i))))/(tand(0.5*(mu(i-1)+Theta(i-1)+mu(i)+Theta(i)))-tand(Theta(i)-mu(i)));
                y(i)=tand(Theta(i)-mu(i))*x(i)+R;
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
                plot([0 x(i)],[R y(i)]);
                hold  on
            end
        end
        i=i+1;
        hold on
    end
    % remaining characteristics
    h=i;
    k=0;
    i=h;
    for j=1:n-1
        while (i<=h+n-k-1)
            % centerline node
            if (i==h)
                x(i)=x(i-n+k)-y(i-n+k)/(tand(0.5*(Theta(i-n+k)+Theta(i)-mu(i-n+k)-mu(i))));
                y(i)=0;
                plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
                hold on
            else 
                % wall node
                if (i==h+n-k-1)
                    x(i)=(x(i-n+k)*tand(0.5*(Theta(i-n+k)+Theta(i)))-y(i-n+k)+y(i-1)-x(i-1)*tand((Theta(i-1)+Theta(i)+mu(i-1)+mu(i))*0.5))/(tand(0.5*(Theta(i-n+k)+Theta(i)))-tand((Theta(i-1)+Theta(i)+mu(i-1)+mu(i))*0.5));
                    y(i)=y(i-n+k)+(x(i)-x(i-n+k))*tand(0.5*(Theta(i-n+k)+Theta(i)));
                    bellWall = [bellWall;x(i),y(i)];
                    plot([x(i-1) x(i)],[y(i-1) y(i)]);
                    hold on
                    plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
                    hold on
                % intermediate nodes
                else
                    s1= tand(0.5*(Theta(i)+Theta(i-1)+mu(i)+mu(i-1)));
                    s2= tand(0.5*(Theta(i)+Theta(i-n+k)-mu(i)-mu(i-n+k)));
                    x(i)=(y(i-n+k)-y(i-1)+s1*x(i-1)-s2*x(i-n+k))/(s1-s2);
                    y(i)=y(i-1)+(x(i)-x(i-1))*s1;
                    plot([x(i-1) x(i)],[y(i-1) y(i)]);
                    hold on
                    plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
                    hold on
                end
            end
            i=i+1;
        end
        k=k+1;
        h=i;
        i=h;
        hold on
    end
    
    % display
    title(sprintf('Characteristic lines for Mach=%.2f and Cp/Cv=%.2f',Me,gamma))
    if scaling == 'ND'
        xlabel('x/x0');
        ylabel('y/y0');
    else
        xlabel('x');
        ylabel('y');
    end
    axis equal
    xlim([0 x(nodes)+0.5])
    ylim([0 y(nodes)+0.5])
end