% clc
close all
clear all


alpha = 0:0.01:3;
number_modes = 5;
Rc = zeros(size(alpha));Rp=Rc;Rh=Rc;
Tc = ones(size(alpha));Tp=Tc;Th=Tc;
% we know the solution for alpha = 0 is R=0, T=1;

L=1;%radius of the circle or half the width
left_length = -5*L;right_length = 5*L;
h=1;Delta = 0.04; 
k_number = 100;k_max = 10*pi;x_number =400;%k_number = 1000;
number_modes=10;dof='one';mode='pitch';
density=0.5; % density of the body relative to water
scale=0;   % scale = 0 - rigid plate, scale = 1 - circular cylinder
fsindx=x_number/2;% fsindx is the index for the node at which we measure the free-surface elevation
 xx = zeros(length(alpha),2);
%Geometric configuration of domain
if(scale==0) % thin rigid plate of width 2L
    l=-L;
    r=L;
    x_left = left_length:L/x_number:l;
    x_right = r:L/x_number:right_length;
    [panelsL,panelsR,panels0,panelsB,panelsF] = ...
      circlebody_twod(2*L,0,l,r,h,Delta,scale);
elseif(scale==1) % semi-submerged cylinder of radius L
    DeltaCyl=10*Delta; % space around cylinder to reduce singular effects
    l=-L-DeltaCyl;
    r=L+DeltaCyl;
    x_left = left_length:L/x_number:l;
    x_right = r:L/x_number:right_length;
    [panelsL,panelsR,panels0,panelsB,panelsF,Nfl,Nfr] = ...
      circlebody_twod(2*L,0,l,r,h,Delta,scale);
    for inc=1:Nfr
       x_rin(inc) = (panelsF(inc,1)+panelsF(inc,3))/2; 
    end
    for inc=1:Nfl
       x_lin(inc) = (panelsF(Nfr+inc,1)+panelsF(Nfr+inc,3))/2; 
    end
end

if(scale==0) % thin rigid plate geoemtrical parameters
  vol=(pi*L^2); % plate has no volume so we set the same as the semicircle
  m = density*vol;
  Im = density*2*(L^3)/3;
  Ia = 2*(L^3)/3;
  I3v = 0;
  W = 2*L;
  Iw = Ia +I3v;
elseif(scale==1) % semi-submerged cylinder geometrical parameters
  vol = (pi*L^2);
  m = density*vol;
  Im = density*pi*(L^4)/2;
  Ia = 2*(L^3)/3;
  I3v = -2*(L^3)/3;
  W = 2*L;
  Iw = Ia +I3v;
end


for j = 2:length(alpha)
   j
[xi,R,T,poten,fpoten]=...
    rigid_body_motion(alpha(j),h,l,r,panelsL,panelsR,panels0,panelsB,panelsF,...
    m,Im,W,Iw,number_modes,'coupled');
Rc(j) = R(1);Tc(j) = T(1);
% 
% [xi,R,T,poten,fpoten]=...
%     rigid_body_motion(alpha(j),h,l,r,panelsL,panelsR,panels0,panelsB,panelsF,...
%     m,Im,W,Iw,number_modes,'heave');
% Rh(j) = R(1);Th(j) = T(1);
% 
% [xi,R,T,poten,fpoten]=...
%     rigid_body_motion(alpha(j),h,l,r,panelsL,panelsR,panels0,panelsB,panelsF,...
%     m,Im,W,Iw,number_modes,'pitch');
%  Rp(j) = R(1);Tp(j) = T(1);   

 xx(j,:) = xi;

end

set(gca,'FontSize',18,'LineWidth',2)
plot(alpha,abs(Rc),'r-',alpha,abs(Tc),'r--', ...
alpha,abs(Rh),'b-',alpha,abs(Th),'b--', ...
alpha,abs(Rp),'k-',alpha,abs(Tp),'k--','LineWidth',2)
xlabel('\alpha','FontSize',18)
ylabel('|R|, |T|','FontSize',18)
axis([0,alpha(end),0,1])

print -djpeg wave_bem_example_floating_RT