% clear all
% close all
format long
font = 'Arial';
fontsize = 30;
linewidth = 3;

theta = 90;
theta = pi/180*theta;
Y20 = 0.25*sqrt(15/2/pi)*sin(theta)^2;
pct = 0.99; %pct of peak for peak finding
plots = 1;  %if 1, show amplitudes and frequencies, else show fractional errors
tcut = 200;
tmin = -20;
tmax = 50;
N = 1000;

af = 0.69;
Mf = 0.95;

t = (linspace(tmin,tmax,N))';

%the following is just a fit I came up with, the better solution is in isco.m
Omisco = (-0.09193*af + 0.09756)./(af.^2 - 2.423*af + 1.437);
omqnm = 1-0.63*(1-af).^0.3;  %Echeverria
Q = 2*(1-af).^(-0.45);
tau = Q./omqnm;
omqnm = omqnm./Mf;
tau = tau.*Mf;

Tau = 2*tau;
Omqnm = omqnm/2;

%Amplitude
% Ap = 0.9073*(1-Mf).^0.7989;
% Ap = 0.91*(1-Mf).^0.8;
Ap = 1.068*(1-Mf).^0.8918;
% Ap = 1.07*(1-Mf).^0.892;

%Frequency
Omref = Omisco.*(1-af);
%Omref = Omisco;
%Omref = 0;
Omp = ((Omqnm.^4 + Omref.^4)/2).^(1/4);
Omm = ((Omqnm.^4 - Omref.^4)/2).^(1/4);

kappap = (Omp.^4-Omm.^4).^(1/4);
kappam = (Omp.^4+Omm.^4).^(1/4);
    
%BOB memory
Om = (Omp^4 + Omm^4*tanh(t/Tau)).^(1/4);
%NB:if Omref = 0, then Om = Omqnm*((1 + tanh(t/Tau))/2)^(1/4);
membob = Ap^2*tau*Om.^2/Omm^4;
membob = membob+membob(2)-2*membob(1);
membob = (sin(theta))^2*(17 + (cos(theta))^2)/(192*pi)*membob;
h20bob = membob/Y20;
%h20bob = h20bob + 4/7*sqrt(5*pi/6)*eta/R; %Eq. 5.16 of http://arxiv.org/pdf/0812.0069v2.pdf
membob = Y20*h20bob;

plot(t,membob,'LineWidth',linewidth)
xlabel('$t/M$','FontSize',fontsize,'Interpreter','latex')
ylabel('$rh_{\rm mem}/M$','FontSize',fontsize,'Interpreter','latex')
xlim([tmin tmax])
set(gca,'FontSize',fontsize,'LineWidth',linewidth)
