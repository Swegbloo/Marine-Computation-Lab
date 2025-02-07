function [Phi,Phi_n,refl,trans,MI,f] = fixed_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsF)
% [Phi,Phi_n,refl,trans] = fixed_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsF)
%
% or
%
% [Phi,Phi_n,refl,trans,MI,f] = fixed_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsF)
% 
% This program returns the potential Phi and its normal derivative 
% The domain consists of a left and right hand vertical part (which must
% have the same length) and a surface, seafloor and body. The boundary
% condition on the body and seafloor is that the normal derivative of
% potential vanishes, the conition on the surface is Phi_n = alpha*Phi.
%
%
% Input
% alpha - the square frequency,
% number_modes - the number of modes taken into account for the evanescents modes,
% N - the array of panel numbers over all segements,
% The panels are
% panelsL is the left hand vertial (1)
% panelsR is the right hand panels (2)
% panels0 is the surface on which the normal vanishes
% panelsF is the free surface
%
% Output
%
% Phi is the potential 
% Phi_n is normal derivative of potential
% refl is the coefficients of the expansion of the potential to the left
% trans is the coefficients of the expansion of the potential to the right
% MI is the matrix lhs 
% f is the rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%z=0     Phi_n=alpha*Phi                                Phi_n=alpha*Phi
%-----.-------------------|                           |----------.---------
%     |       (4) a       /                        __/   (4)     |
%     |                  /          ____ (3)      /              |
%     |x=-1/2            \         /    \        /               |x=1/2
%     |                   \____   /     |  _____/                |
%     |                        \_/       \/                      |
%  (1)|                         Phi_n=0                          |
%     |                                                          |
%     |                                                          |Phi_n=Q2*Phi
%     |                                                          |
%     |Phi_n=Q1*Phi                                              |
%     |                                                          |
%     |                                                          |
%     |                                                          |
%     |                                                          |(2)
%     |                                                          |
%     |                                                          |
%     |                                                          |
%     |                                                          |
%     |                                                          |
%-----.----------------------------------------------------------.---------
% z=-h                  (3)             Phi_n=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The equations we are going to solve are | 1/2*Phi = H*Phi - G*Phi_n
%                                        | Phi_n = Q*Phi + F
%
%Hence we can find the potential over the boundaries, solving the simple
%equation (1/2I-H+GQ)Phi=-GF

% Extract information from panel arrays
% Output panels grouped in blocks according to the boundary condition
panels=[panelsL;
    panelsR;
    panels0;
    panelsF];
NL=size(panelsL,1);
NR=size(panelsR,1);
N0=size(panels0,1);
NF=size(panelsF,1);

%Firstly, we can generate the matrices G and H
[G,H]=bem_constant_panel(panels);

%implementation of the matrix Q (boundary conditions matrix) and the vector
%F (incident wave vector)
[Q_left,F,S,R,K]=QF_matrices(panelsL,alpha,number_modes); %method 1

Q_right=flipud(fliplr(Q_left));

NTOT=NL+NR+N0+NF;

%Implementation of the full matrix Q
Q_big=[Q_left,zeros(NL,NTOT-NL);
    zeros(NL,NL),Q_right,zeros(NL,NTOT-2*NL);
    zeros(N0,NTOT);
    zeros(NF,NTOT-NF),alpha*eye(NF)];

F_big=[F;zeros(NTOT-NL,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solving
MI = (1/2*eye(NL+NR+N0+NF)-H+G*Q_big);
f=G*F_big;
Phi=MI\f;
Phi_n=Q_big*Phi-F_big;


% calculating r (noting we have to remove the incident wave) and t
refl = R*(Phi(1:NL) - 1/(-2*K(1))*F);
trans = flipud(flipud(fliplr(R))*Phi(NL+1:NL+NR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% function QF_matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Q,F,S,R,K]=QF_matrices(panelsL,alpha,number_modes)
        %
        %This program returns the matrix Q and the vector F, which represent
        %respectively the boundary conditions over all the segments, and the
        %incident wave.
        % Inputs
        % panels - the matrix defining the meshing of the boudary of the domain,
        % alpha  - the square frequency,
        % number_modes  - the number of real roots of dispersion equation
        % NL,NR,N0,NF - number of panels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NL=size(panelsL,1);
        l = -panelsL(1,1); % position of the left hand edge is -l
        h = - panelsL(end,4);

        %We define the calculation points in the middle of each panel
        mid_points_z = (panels(1:NL,2) + panels(1:NL,4))/2;

        %Obtaining of the roots of the dispersion equation
        K=dispersion_free_surface(alpha,number_modes,h);
        M=length(K);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculation of the matrix Q :

        % M is the number of evanescent modes
        A=zeros(1,M);
        for j=1:M,
            A(j)=0.5*(cos(K(j)*h)*sin(K(j)*h)+K(j)*h)/(K(j)*(cos(K(j)*h)^2));
        end

        %We begin with the boundary condition on the left hand of the domain
        S = zeros(NL,M);
        R = zeros(M,NL);

        for k = 1 : NL,
            for j = 1 : M,
                S(k,j)=-K(j)*cos(K(j)*(mid_points_z(k)+h))/(cos(K(j)*h));
            end
        end

        for j = 1 : M,
            for k = 1 : NL
                R(j,k) = 1/(K(j)*A(j))*(sin(K(j)*(panels(k,2)+h))-sin(K(j)*(panels(k,4)+h)))...
                    /cos(K(j)*h);
            end
        end
        Q = S*R;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Calculation of the vector F :
        F=zeros(NL,1);

        for k=1:NL,
            F(k) = -2*exp(l*K(1))*K(1)*cos(K(1)*(mid_points_z(k)+h))/cos(K(1)*h);
        end



    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
