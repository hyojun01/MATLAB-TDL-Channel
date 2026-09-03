clear; close all; clc;

%% Parameter
c = 299792458;                              % velocity of light
fc = 24e9;                                  % center frequency
B = 100e6;                                  % chirp bandwidth
T = 40.96e-6;                               % chirp period
S = B/T;                                    % chirp slope
rho = 0.3;                                  % reflectivity
N = 128;                                    % the number of chirps

%% Modeling Beat signal with multipath effect
M = 512;                                    % the number of samples in chirp
Fs = M/T;                                   % sample rate
tau = (0:M-1)'/Fs;                          
n = 0:N-1;
d0 = 25;
dr = 10;
dt = 8.4087335455;
v = 0;
r = sqrt(d0^2 + (dt-dr)^2);
rp = sqrt(d0^2 + (dt+dr)^2);
range = [r, (r+rp)/2, rp];
velocity = [v*d0/r, v*d0*(1/r+1/rp)/2, v*d0/rp];
deltaR = range(:)-range(1);
deltaV = velocity(:)-velocity(1);

fb = 2*S*range(1)/c - 2*velocity(1)*fc/c;
fv = 2*velocity(1)*fc/c;
Kr_dB = fspl(r,fc/c);
Kr = 10^(Kr_dB/20);
y = Kr*exp(1j*2*pi*fb*tau).*exp(-1j*2*pi*fv*n*T) ...
    *exp(1j*4*pi*fc*range(1)/c);
dfb = 2*S*deltaR/c;
dfv = 2*deltaV*fc/c;
s = y.*(1 ...
        + 2*rho*exp(1j*2*pi*dfb(2)*tau).*exp(-1j*2*pi*dfv(2)*n*T) ...
            *exp(1j*4*pi*fc*deltaR(2)/c) ...
        + rho^2*exp(1j*2*pi*dfb(3)*tau).*exp(-1j*2*pi*dfv(3)*n*T) ...
            *exp(1j*4*pi*fc*deltaR(3)/c));

% Windowed 2D FFT (toolbox-free Hann windows).
wf = 0.5-0.5*cos(2*pi*(0:M-1)'/(M-1));
ws = 0.5-0.5*cos(2*pi*(0:N-1)/(N-1));
RD = fftshift(fft(fft(s.*(wf*ws),M,1),N,2),2);
RD = fliplr(RD(1:M/2,:));
RDdb = 20*log10(abs(RD)/max(abs(RD),[],'all') + eps);
rangeAxis = (0:M/2-1)*Fs/M*c/(2*S);
fdAxis = (-N/2:N/2-1)/(N*T);
velocityAxis = fliplr(-fdAxis*c/(2*fc));

%% Results
% All plots use the current d0, dr and dt parameters above.
Y = fft(y(:,1).*wf,M);
Smp = fft(s(:,1).*wf,M);
frequencyAxis = (0:M/2-1)*Fs/M/1e6;
reference = max(abs(Y(1:M/2)));
Ydb = 20*log10(abs(Y(1:M/2))/reference + eps);
Sdb = 20*log10(abs(Smp(1:M/2))/reference + eps);

tauShot = 128;
nShot = 64;
d0Axis = linspace(max(1,d0/4),max(150,2*d0),10000);
rAxis = sqrt(d0Axis.^2 + (dt-dr)^2);
rpAxis = sqrt(d0Axis.^2 + (dt+dr)^2);
range2Axis = [rAxis; (rAxis+rpAxis)/2; rpAxis];
velocity2Axis = [v*d0Axis./rAxis; ...
                 v*d0Axis.*(1./rAxis+1./rpAxis)/2; ...
                 v*d0Axis./rpAxis];
deltaRAxis = [range2Axis(1,:)-range2Axis(1,:); ...
              range2Axis(2,:)-range2Axis(1,:); ...
              range2Axis(3,:)-range2Axis(1,:)];
deltaVAxis = [velocity2Axis(1,:)-velocity2Axis(1,:); ...
              velocity2Axis(2,:)-velocity2Axis(1,:); ...
              velocity2Axis(3,:)-velocity2Axis(1,:)];
dfbAxis = 2*S*deltaRAxis/c;
dfbAxis = zeros(3,10000);
dfvAxis = 2*deltaVAxis*fc/c;
dfvAxis = zeros(3,10000);
amplitudeFactor = abs(1 + 2*rho*exp(1j*2*pi*dfbAxis(2,:)*tauShot) ...
                            .*exp(-1j*2*pi*dfvAxis(2,:)*nShot*T) ...
                            .*exp(1j*4*pi*fc*deltaRAxis(2,:)/c) + ...
                      + rho^2*exp(1j*2*pi*dfbAxis(3,:)*tauShot) ...
                            .*exp(-1j*2*pi*dfvAxis(3,:)*nShot*T) ...
                            .*exp(1j*4*pi*fc*deltaRAxis(3,:)/c));
amplitudeDb = 20*log10(amplitudeFactor + eps);

fig = figure(Color='w',Name='FMCW Multipath Effects');
tl = tiledlayout(fig,2,2,TileSpacing='compact',Padding='compact');

ax1 = nexttile(tl);
timeReference = mean(abs(y(:,1)));
plot(ax1,tau*1e6,abs(y(:,1))/timeReference,'--',LineWidth=1.2)
hold(ax1,'on')
plot(ax1,tau*1e6,abs(s(:,1))/timeReference,LineWidth=1.2)
hold(ax1,'off')
grid(ax1,'on')
xlabel(ax1,'Fast time (\mus)')
ylabel(ax1,'Amplitude (LOS reference)')
title(ax1,sprintf('Time-domain fading, d_0 = %.1f m',d0))
legend(ax1,'LOS','Multipath',Location='best')

ax2 = nexttile(tl);
plot(ax2,frequencyAxis,Ydb,'--',LineWidth=1.2)
hold(ax2,'on')
plot(ax2,frequencyAxis,Sdb,LineWidth=1.2)
hold(ax2,'off')
grid(ax2,'on')
xlabel(ax2,'Beat frequency (MHz)')
ylabel(ax2,'Magnitude (dB, LOS reference)')
title(ax2,'Frequency-domain fading')
legend(ax2,'LOS','Multipath',Location='best')
ylim(ax2,[-50 10])

ax3 = nexttile(tl);
plot(ax3,d0Axis,amplitudeDb,LineWidth=1.2)
grid(ax3,'on')
xline(ax3,d0,'--',sprintf('d_0 = %.1f m',d0), ...
    LabelOrientation='horizontal')
xlabel(ax3,'d_0 (m)')
ylabel(ax3,'Amplitude ratio (dB)')
title(ax3,'Amplitude variation versus d_0')
ylim(ax3,[-25 12])

ax4 = nexttile(tl);
imagesc(ax4,velocityAxis,rangeAxis,RDdb)
axis(ax4,'xy')
clim(ax4,[-45 0])
xlim(ax4,[min(velocity)-2 max(velocity)+2])
ylim(ax4,[range(1)-2 range(3)+2])
xlabel(ax4,'Radial velocity (m/s)')
ylabel(ax4,'Range (m)')
title(ax4,sprintf('Range-Doppler map, d_0 = %.1f m',d0))
colorbar(ax4)
hold(ax4,'on')
plot(ax4,velocity(1),range(1),'ko',MarkerSize=8,LineWidth=1.5)
plot(ax4,velocity(2:3),range(2:3),'kx',MarkerSize=9,LineWidth=1.5)
hold(ax4,'off')
legend(ax4,'Real','Ghost',Location='southwest')
