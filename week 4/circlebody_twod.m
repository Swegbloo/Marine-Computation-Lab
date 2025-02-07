function [panelsL,panelsR,panels0,panelsB,panelsF,vol,W] = ...
    circlebody_twod(D,a,l,r,h,Delta,scale)
% INPUT
% calculates the nodes for a single plate
% D, the extent of the plate
% a - the offsets relative to the origin of the centre of the plate
% l and r the left and right ends of the domain
% h is the water depth
% Delta is the (approximate) panel length
% Scale is an optional command which scales the circle (so that scale = 0
% corresponds to a dock)

%OUTPUT
% panelsL -- panels for left control surface
% panelsR -- panels for right control surface
% panels0 -- panels for the body
% panelsB -- body panels
% panelsF -- free-surface panels
% V is the volume and W is the wetted surface area
if nargin == 8
    scale = 1;
end


% first we calculate the number of nodes
Nb = max(round(pi/2*D/Delta),1); % nodes on plate
Nfl = max(round((-l + a - D/2)/Delta),0); % nodes on left free surface (these can be zero)
Nfr = max(round((r - a - D/2)/Delta),0); % nodes on right free surface (these can be zero)
Nc  = max(round(h/Delta),1); % nodes on vertical surfaces
Nbed = max(round((r+-l)/Delta),1); % nodes on vertical surfaces

[bodypanels_plate,xl,xr,vol,W]=body_shape(Nb,D,a,scale);



%%%%%%% Creat intermediate free surface panels %%%%%%%%%%%%%%%%%%%%%%%%%%%


dxr=(xr-r)/(Nfr);
dxl=(l-xl)/(Nfl);
xr=r:dxr:xr;  % Nfr+1 nodes, Nfr panels
xl=xl:dxl:l; % Nfl+1 nodes, Nfl panels
zf=0;
z1 = 0:-h/Nc:-h;
xbed=l:(r+-l)/Nbed:r;
z2=-h:h/Nc:0;

% Panels L
panelsL=[xbed(1)*ones(Nc,1),z1(1:Nc).',xbed(1)*ones(Nc,1),z1(2:Nc+1).'];
% Panels R
panelsR=[xbed(Nbed+1)*ones(Nc,1),z2(1:Nc).',xbed(Nbed+1)*ones(Nc,1),z2(2:Nc+1).'];
% Panels Rigid Surface
panels0=[xbed(1:Nbed).',z1(Nc+1)*ones(Nbed,1),xbed(2:Nbed+1).',z1(Nc+1)*ones(Nbed,1)];
panelsB=[bodypanels_plate];
% Panels F
panelsF=[ xr(1:Nfr).',zf*ones(Nfr,1),xr(2:Nfr+1).',zf*ones(Nfr,1);...
    xl(1:Nfl).',zf*ones(Nfl,1),xl(2:Nfl+1).',zf*ones(Nfl,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [bodypanels,xl,xr,vol,W]=body_shape(Nb,D,a,scale)

        z_body=zeros(1,Nb);    % initialisation

        
        %your body shape
        for i=1:Nb+1,
            z_body(i) = a+(D/2)*(cos((i-1)*pi/Nb)-j*sin((i-1)*pi/Nb));
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bodypanels=[real(z_body(1:Nb))',scale*imag(z_body(1:Nb))',...
            real(z_body(2:Nb+1))',scale*imag(z_body(2:Nb+1))'];
        xl=bodypanels(Nb,3);
        xr=bodypanels(1,1);
        
        %Mass and waterplane area of the body 
        if(scale==1) 
          % semicircle
          vol = (pi*(D/2)*(D/2));
          W = D;
        elseif(scale==0)
          % flat plate
          vol=D;
          W = D;
        end

    end

    function [N]=bodypanels(eqbodynode)

        %load('nodes.in', Nb, Nb1, Nb2)
        if eqbodynode==1,
            Nb=input('Enter the number of panels on each body: ')
        else
            Nbr=input('Enter Nbr, number of panels on the right body: ')
            Nbl=input('Enter Nbl, number of panels on the left body: ')
        end

        % free-surface nodes
        [Nf]=input('Enter Nfl, Nfm, Nfr free-surface panel numbers: ')

        Nfl=Nf(1);
        Nfm=Nf(2);
        Nfr=Nf(3);

        Nc=input('Enter Nc, the number of panels on the control surfaces: ')
        Nbed=input('Enter the number of panels on the sea bed: ')
        N=[Nc,Nbed,Nc,Nfl,Nb,Nfm,Nb,Nfr]
    end


end