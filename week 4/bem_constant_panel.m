function [G,Gn]=bem_constant_panel(panels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Comments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function gives the integration of the Green's function G and its 
%normal derivative. The only parameter is the matrix 'panels' which gives 
%the coordinates of each extremity of each panel all along the boundary
%meshing. We will give a simple ewemple of what this matrix could be, but
%firstly we will present the organisation of 'panels'. This matrix must
%have 4 columns. The 2 first columns define the coordinates (x,z) of the
%first extremity of the current panel with respect of the route direction.
%Of course the 2 second columns define the coordinates of the other
%extremity of the panel. Hence each line of the matrix represents one
%panel. Then this matrix is organized with respect of the route direction
%of the boundary line treatment, namely for all 1<i<=N (N number of lines),
%the first extremity coordinates of line i (panels(i,1:2)) is equal to the 
%end coordinates of line i-1 (panels(i-1,3:4)), and of course 
%panels(1,1:2)=panels(N,3:4) because the boundary is closed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Exemple of matrix 'panels'%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We consider here a side length 2 square centered in the origin as decribed
%in the following figure.
%                           .
%                          /|\ normal
%     first panel           |
%     x=-1 .------------------------------------.x=1
%     z=1  |                                    |z=1
%          |                                    |
%          |                                    |
%          |                                   /|\ route direction
%          |                 z                  |
%          |                 |_ x               |
%         \|/               0                   |
%          |                                    |
%          |                                    |
%          |                                    |
%          |                                    |
%          |                                    |
%     x=-1 .------------------------------------.x=1
%     z=-1                                       z=-1
%
%What is really important to define is the direction we run all the panels
%so that the normals of all panels are oriented towards the same way.
%Moreover, we begin the treatment with the panel at the corner top on the
%left. We discretize here the domain with a constant panel length (not an
%obligation) and we generate n panels over each side of the square.
%
%We define 2 vectors x and z for the discretization of the domain
%z = 1:-2/n:-1;
%x=-z;
%m = length(x) - 1;

%The definition of the meshing is the part where we need to be very
%careful because of the orientation of the normals. 
%panels = [x(1)*ones(n,1),z(1:n).',x(1)*ones(n,1),z(2:n+1).';...
%        x(1:m).',z(n+1)*ones(m,1),x(2:m+1).',z(n+1)*ones(m,1);...
%         x(m+1)*ones(n,1),-z(1:n).',x(m+1)*ones(n,1),-z(2:n+1).';...
%        -x(1:m).',z(1)*ones(m,1),-x(2:m+1).',z(1)*ones(m,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Implementation of matrices G and Gn%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We can now calculate the matrices G and Gn. To do that, we call the 
%functions SlantedG and SlantedGnormal which respectively integrates The
%Green's function and its normal derivative.
q = size(panels,1);
G = zeros(q,q);
Gn = zeros(q,q);

for j=1:q
    for k=1:q
        G(j,k) = SlantedG(panels(k,:),panels(j,:));
        Gn(j,k) = SlantedGnormal(panels(k,:),panels(j,:));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Calculation of Phi_n knowing Phi (optional)%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do not forget to put Phi_n as a result of the function
%Here we have Phi(x,z)=x

%Phi=zeros(q,1);
%for i=1:q,
%    x=(panels(i,1)+panels(i,3))/2;
%    z=(panels(i,2)+panels(i,4))/2;
%    Phi(i)=x;
%end

%Phi_n = inv(G)*(Gn*Phi-(1/2).*Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function SlantedG%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = SlantedG(xpanel,xsource)
% f = SlantedG(xsource,xpanel,panel_length)
% xsource is the end points of the source panel and xpanel is the end
% points of the panel over which we are integrating the Green's function. 
% these are with respect to the cartestian co-ordinate system. 

% we have two co-ordinate systems, the original (which is the one
% which xsource and xpanel and in) and the new co-ordinate system
% in which the integration is performed. 

%let's calculate the length of the panel
panel_length = sqrt((xpanel(1) - xpanel(3))^2 + (xpanel(2) - xpanel(4))^2);

%int the new co-ordinates the end points of the panel are at plus or minus
%the length of the panel divided by 2.
xi1 = -panel_length /2;
xi2 = panel_length /2;

%we now convert the source midpoint into the new co-ordinate system.
%first of all we find the angle.
deltax = xpanel(3)-xpanel(1);
deltaz = xpanel(4)-xpanel(2);

theta = atan2(deltaz,deltax);
 
%now we find the distance between the midpoints of the panels
mid_panel = [(xpanel(1) + xpanel(3))/2,(xpanel(2) + xpanel(4))/2];
mid_source = [(xsource(1) + xsource(3))/2,(xsource(2) + xsource(4))/2];

%now we calculate the source points with respect to the new co-ordinate 
%system.
xnew = (mid_source(1) - mid_panel(1))*cos(theta) + (mid_source(2) - ...
    mid_panel(2))*sin(theta);
znew = (mid_source(1) - mid_panel(1))*sin(theta) - (mid_source(2) - ...
    mid_panel(2))*cos(theta);

out = Greens_xi(xi1,xi2,xnew,znew);

function out = Greens_xi(xi1,xi2,x,z)
% xi1 and xi2 the end point of the panel in the new co-ordinates
% and x and z and the co-ordinates of the mid point of the source
% in the new co-ordinates. 

out = Integrand(xi2,x,z) - Integrand(xi1,x,z);

function f = Integrand(xi,x,z)
%this is the value of the integral of 1/(4*pi)*log((x1 - x2)^2+(z1 - z2)^2) 
if x == xi
      f = -1/(2*pi)*xi;
else
    if z == 0
        f = 1/(4*pi)*log((x-xi).^2+(z)^2)*(xi-x) - 1/(2*pi)*xi;
    else
        f = 1/(4*pi)*log((x-xi).^2+(z)^2)*(xi-x) - 1/(2*pi)*xi - ...
            1/(2*pi)*atan((x-xi)/(z))*z;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%function SlantedGnormal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = SlantedGnormal(xpanel,xsource)
% f = SlantedG(xsource,xpanel)
% xsource is the end points of the source panel and xpanel is the end
% points of the panel over which we are integrating the normal derivative 
%of Green's function. These are with respect to the cartestian co-ordinate
%system. The normal direction is taken to point outward on the the right 
%hand side of the panel which stretches from (x1,z1) to (x2,z2)

% we have two co-ordinate systems, the original (which is the one
% which xsource and xpanel and in) and the new co-ordinate system
% in which the integration is performed. 

%lets calculate the length of the panel
panel_length = sqrt((xpanel(1) - xpanel(3))^2 + (xpanel(2) - xpanel(4))^2);

%int the new co-ordinates the end points of the panel are at plus or minus
%the length of the panel divided by 2.
xi1 = -panel_length /2;
xi2 = panel_length /2;

%we now convert the source midpoint into the new co-ordinate system.
%first of all we find the angle.
deltax = xpanel(3)-xpanel(1);
deltaz = xpanel(4)-xpanel(2);

theta = atan2(deltaz,deltax);

%now we find the distance between the midpoints of the panels
mid_panel = [(xpanel(1) + xpanel(3))/2,(xpanel(2) + xpanel(4))/2];
mid_source = [(xsource(1) + xsource(3))/2,(xsource(2) + xsource(4))/2];

%now we calculate the source points with respect to the new co-ordinate 
%system.
xnew = (mid_source(1) - mid_panel(1))*cos(theta) + (mid_source(2) - ...
    mid_panel(2))*sin(theta);
znew = (mid_source(1) - mid_panel(1))*sin(theta) - (mid_source(2) - ...
    mid_panel(2))*cos(theta);

% the normal direction is taken to point outward on the the right hand 
%side of the panel which stretches from (x1,z1) to (x2,z2).
% This is the value of the integral of -1/(2*pi) * (atan(x2 - x1)/z)
%znew = - znew;
if abs(znew) <= 1e-12, % due to rounding error, we assume 1e-12 is close 
                       % enough to zero
   out = 0;
else
   out = 1/(2*pi)*(atan((xnew-xi2)/znew)-atan((xnew-xi1)/znew));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


