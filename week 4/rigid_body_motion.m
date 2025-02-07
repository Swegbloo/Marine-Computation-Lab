function [xi,R,T,poten,fpoten]=...
    rigid_body_motion(alpha,h,l,r,panelsL,panelsR,panels0,panelsB,panelsF,...
    m,Im,W,Iw,number_modes,mode)
% program to calculate the motion of a floating body
%--------------------------------------------------------------------------
% Variables:
% alpha  = frequency squared
% h      =  water depth
% l      =  left end of bem domain
% r      =  right end of bem domain
% panelsL to panelsF  = contain specifications of the panel data on all boundaries
%
% m      = mass of the structure
% Im     = moment of inertia of the structure
% W      = hydrostatic restoring coefficent
% number_modes = number of modes used outside
% mode  = 'heave' or 'pitch' for separate or 'coupled' if both are needed.
%
% xi = displacement of body first component is heave, second is pitch.
% R = reflection
% T = transmission
% poten = potential on the body panels
% fpoten = complex potential on free surface



NL=size(panelsL,1);
NR=size(panelsR,1);
N0=size(panels0,1);
NB=size(panelsB,1);
NF=size(panelsF,1);

% Calculate the geometrical parameters for computing forces
n3_dS = zeros(1,size(panelsB,1));
n5_dS = zeros(1,size(panelsB,1));
for inc=1:size(panelsB,1)
    delta_x = panelsB(inc,1)-panelsB(inc,3);
    delta_z = panelsB(inc,2)-panelsB(inc,4);
    x_mean = (panelsB(inc,1)+panelsB(inc,3))/2;
    z_mean = (panelsB(inc,2)+panelsB(inc,4))/2;
    n3_dS(inc) = delta_x;
    n5_dS(inc) = -z_mean*delta_z-x_mean*delta_x;
end

roots = dispersion_free_surface(alpha,number_modes,h);

% calculate the motion of the structure and fluid
if strcmp(mode ,'heave')
    % solve for the diffraction potential
    %panels0=[panels0;panelsB];
    [PhiS,PhiS_n,refl,trans] =...
        fixed_body_twod(alpha,number_modes,panelsL,panelsR,[panels0;panelsB],panelsF);

    [Phi3R,Phi3R_n,rad3_left,rad3_right] = ...
        forced_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsB,panelsF,mode);
    % solve for the diffraction potential
    Fad3=n3_dS*Phi3R(NL+NR+N0+1:NL+NR+N0+NB); % matrix multiplication
    X3=-sqrt(-1)*sqrt(alpha)*(n3_dS*PhiS(NL+NR+N0+1:NL+NR+N0+NB));
    xi3 = X3/(W-alpha*(m+Fad3));
    poten = sqrt(-1)*sqrt(alpha)*xi3*Phi3R(NL+NR+N0+1:NL+NR+N0+NB)...
        + PhiS(NL+NR+N0+1:NL+NR+N0+NB);
    fpoten = sqrt(-1)*sqrt(alpha)*xi3*Phi3R(NL+NR+N0+NB+1:NL+NR+N0+NB+NF)...
        + PhiS(NL+NR+N0+NB+1:NL+NR+N0+NB+NF);
    R = (refl+sqrt(-1)*sqrt(alpha)*xi3*rad3_left);
    T = (trans+sqrt(-1)*sqrt(alpha)*xi3*rad3_right);
    xi = [xi3;0];
    
elseif strcmp(mode ,'pitch')
    % solve for the diffraction potential
    %panels0=[panels0;panelsB];
    [PhiS,PhiS_n,refl,trans] =...
        fixed_body_twod(alpha,number_modes,panelsL,panelsR,[panels0;panelsB],panelsF);
    [Phi5R,Phi5R_n,rad5_left,rad5_right] = ...
        forced_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsB,panelsF,mode);
    % solve for the diffraction potential
    Fad5=n5_dS*Phi5R(NL+NR+N0+1:NL+NR+N0+NB); % matrix multiplication
    X5=-sqrt(-1)*sqrt(alpha)*(n5_dS*PhiS(NL+NR+N0+1:NL+NR+N0+NB));
    xi5 = X5/(Iw-alpha*(Im+Fad5));
    poten = sqrt(-1)*sqrt(alpha)*xi5*Phi5R(NL+NR+N0+1:NL+NR+N0+NB)...
        + PhiS(NL+NR+N0+1:NL+NR+N0+NB);
    fpoten = sqrt(-1)*sqrt(alpha)*xi5*Phi5R(NL+NR+N0+NB+1:NL+NR+N0+NB+NF)...
        + PhiS(NL+NR+N0+NB+1:NL+NR+N0+NB+NF);
    R = (refl+sqrt(-1)*sqrt(alpha)*xi5*rad5_left);
    T = (trans+sqrt(-1)*sqrt(alpha)*xi5*rad5_right);
    xi = [0,xi5];


elseif strcmp(mode ,'coupled')
    % Solve for the radiation and diffraction potentials
    % heave radiation potential
    [Phi3R,Phi3R_n,rad3_left,rad3_right] = ...
        forced_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsB,panelsF,'heave');
    % pitch radiation potential
    [Phi5R,Phi5R_n,rad5_left,rad5_right] = ...
        forced_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsB,panelsF,'pitch');
    % solve for the diffraction potential
    panels0=[panels0;panelsB];
    [PhiS,PhiS_n,refl,trans] =...
        fixed_body_twod(alpha,number_modes,panelsL,panelsR,panels0,panelsF);


    % Calculate the motions of the structure
    % Heave
    Fad3 = n3_dS*Phi3R(NL+NR+N0+1:NL+NR+N0+NB); % matrix multiplication
    X3 = -sqrt(-1)*sqrt(alpha)*(n3_dS*PhiS(NL+NR+N0+1:NL+NR+N0+NB));
    xi3 = X3/(W-alpha*(m+Fad3));
    % Pitch
    Fad5 = n5_dS*Phi5R(NL+NR+N0+1:NL+NR+N0+NB); % matrix multiplication
    X5 = -sqrt(-1)*sqrt(alpha)*(n5_dS*PhiS(NL+NR+N0+1:NL+NR+N0+NB));
    xi5 = X5/(Iw-alpha*(Im+Fad5));

    % Coupled system (the subscript 3 is used to store the total values)

    poten = sqrt(-1)*sqrt(alpha)*(xi3*Phi3R(NL+NR+N0+1:NL+NR+N0+NB)...
        + xi5*Phi5R(NL+NR+N0+1:NL+NR+N0+NB)) + PhiS(NL+NR+N0+1:NL+NR+N0+NB);
    fpoten = sqrt(-1)*sqrt(alpha)*(xi3*Phi3R(NL+NR+N0+NB+1:NL+NR+N0+NB+NF)...
        + xi5*Phi5R(NL+NR+N0+NB+1:NL+NR+N0+NB+NF)) + PhiS(NL+NR+N0+NB+1:NL+NR+N0+NB+NF);
    % we calculate the reflection and transmission coefficients.
    R = (refl+sqrt(-1)*sqrt(alpha)*(xi3*rad3_left+xi5*rad5_left));
    T = (trans+sqrt(-1)*sqrt(alpha)*(xi3*rad3_right+xi5*rad5_right));
    xi = [xi3;xi5];
end
