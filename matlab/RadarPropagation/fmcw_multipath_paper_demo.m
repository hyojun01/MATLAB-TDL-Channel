%% FMCW multipath: power fading and ghost targets
% Liu et al., Eqs. (1)-(15). Amplitudes are normalized (Kr*A0/2 = 1).
clear; close all; clc

%% Parameter
c = 299792458;                              % velocity of light 
f0 = 77e9;                                  % center frequency
B = 500e6;                                  % chirp Bandwidth
T = 40.96e-6;                               % chirp Period
mu = B/T;                                   % chirp Slope
rho = 0.7;                                  % reflectivity
N = 128;                                    % the number of chirps

%% 1) Vertical multipath power fading: Eq. (15)
d = linspace(40,150,3000);
dr = 0.5;
dt = 0.8;
r = sqrt(d.^2 + (dt-dr)^2);
rt = sqrt(d.^2 + (dt+dr)^2);
deltaR1 = (rt-r)/2;
F = abs(1 + rho*exp(1j*4*pi*f0*deltaR1/c)).^2;          % Eq. (15)
powerFadingDb = 20*log10(F);                            % |s/y|=F, so 10*log10(Ps/Py)=20*log10(F)

%% 2) Horizontal multipath beat signal and range-Doppler map: Eq. (11)
M = 256; 
fs = M/T; 
tau = (0:M-1)'/fs;
n = 0:N-1;
d0 = 20;
dr = 1.75;
dt = 8.75;
v = 20;                                                 % along-road values [m], [m/s]
r = sqrt(d0^2 + (dt-dr)^2); 
rt = sqrt(d0^2 + (dt+dr)^2);
range = [r, (r+rt)/2, rt];                              % Eq. (3)
velocity = [v*d0/r, v*d0*(1/r+1/rt)/2, v*d0/rt];        % Eq. (4)
deltaR = range(2:3)-range(1);
deltaV = velocity(2:3)-velocity(1);

fb = 2*mu*range(1)/c - 2*velocity(1)*f0/c;              % Eqs. (9),(10)
fv = 2*velocity(1)*f0/c;
y = exp(1j*2*pi*fb*tau).*exp(-1j*2*pi*fv*n*T) ...
    *exp(1j*4*pi*f0*range(1)/c);
df = 2*mu*deltaR/c; 
dfv = 2*deltaV*f0/c;                % Eqs. (12),(13)
s = y.*(1 ...
    + 2*rho*exp(1j*2*pi*df(1)*tau).*exp(-1j*2*pi*dfv(1)*n*T) ...
        *exp(1j*4*pi*f0*deltaR(1)/c) ...
    + rho^2*exp(1j*2*pi*df(2)*tau).*exp(-1j*2*pi*dfv(2)*n*T) ...
        *exp(1j*4*pi*f0*deltaR(2)/c));

% Windowed 2D FFT (toolbox-free Hann windows).
wf = 0.5-0.5*cos(2*pi*(0:M-1)'/(M-1));
ws = 0.5-0.5*cos(2*pi*(0:N-1)/(N-1));
RD = fftshift(fft(fft(s.*(wf*ws),M,1),N,2),2);
RD = fliplr(RD(1:M/2,:));
RDdb = 20*log10(abs(RD)/max(abs(RD),[],'all') + eps);
rangeAxis = (0:M/2-1)*fs/M*c/(2*mu);
fdAxis = (-N/2:N/2-1)/(N*T);
velocityAxis = fliplr(-fdAxis*c/(2*f0));

%% Results
figure(Color='w'); tiledlayout(1,2,TileSpacing='compact')
nexttile; plot(d,powerFadingDb,LineWidth=1.2); grid on; ylim([-50 10])
xlabel('Along-road distance d (m)'); ylabel('Received power ratio (dB)')
title('Power fading, \rho = 0.7')

nexttile; imagesc(velocityAxis,rangeAxis,RDdb); axis xy
xlim([min(velocity)-2 max(velocity)+2]); ylim([range(1)-2 range(3)+2]); clim([-45 0])
xlabel('Radial velocity (m/s)'); ylabel('Range (m)'); colorbar
title('Real target and two ghosts'); hold on
plot(velocity(1),range(1),'ko',MarkerSize=8,LineWidth=1.5)
plot(velocity(2:3),range(2:3),'kx',MarkerSize=9,LineWidth=1.5); hold off
legend('Real','Ghost 2/3 and 4',Location='southwest')

disp(table(["Real";"Ghost 2/3";"Ghost 4"],range.',velocity.', ...
    VariableNames={'Target','Range_m','Velocity_mps'}))
