%% Multipath fading in an FMCW received signal
% The beat-domain multiplier is not reused here. Each propagation path is
% applied to the received complex-baseband FMCW signal and summed coherently.
clear; close all; clc

%% Parameters
c = physconst('LightSpeed');
fc = 24e9;                                  % RF center frequency
B = 50e6;                                  % sweep bandwidth
T = 40.96e-6;                               % sweep time
Fs = 8*B;                                   % oversampled complex-baseband rate
rho = 0.3;                                  % reflection coefficient

d0 = 25;                                  
dr = 10;                                   
dt = 8.4087335455;                         
v = 0;                                     

%% FMCW waveform: two sweeps are generated; the second is observed
waveform = phased.FMCWWaveform( ...
    SampleRate=Fs, ...
    SweepTime=T, ...
    SweepBandwidth=B, ...
    SweepDirection='Up', ...
    SweepInterval='Positive', ...
    NumSweeps=2);
xt = waveform();
M = round(Fs*T);
t = (0:M-1)'/Fs;

%% Three equivalent paths: LOS, two single reflections, double reflection
r = sqrt(d0^2 + (dt-dr)^2);
rp = sqrt(d0^2 + (dt+dr)^2);
targetRange = [r,(r+rp)/2,rp];
targetVelocity = [v*d0/r,v*d0*(1/r+1/rp)/2,v*d0/rp];
pathWeight = [1,2*rho,rho^2];

radarPosition = [0;0;0];
radarVelocity = [0;0;0];
rxPath = complex(zeros(size(xt,1),3));

for k = 1:3
    channel = phased.FreeSpace( ...
        SampleRate=Fs, ...
        OperatingFrequency=fc, ...
        TwoWayPropagation=true);
    targetPosition = [targetRange(k);0;0];
    targetMotion = [-targetVelocity(k);0;0];
    rxPath(:,k) = channel(xt,radarPosition,targetPosition, ...
        radarVelocity,targetMotion);
end

rxLosAll = rxPath(:,1);
rxMultipathAll = rxPath*pathWeight.';
observe = M+(1:M);                          % remove the initial delay transient
xrLos = rxLosAll(observe);
xrMultipath = rxMultipathAll(observe);

%% Frequency-domain representation of the received chirp
% window = 0.5-0.5*cos(2*pi*(0:M-1)'/(M-1));
Nfft = 2^nextpow2(2*M);
Xlos = fftshift(fft(xrLos,Nfft));
Xmp = fftshift(fft(xrMultipath,Nfft));
fBaseband = (-Nfft/2:Nfft/2-1)'*Fs/Nfft;
spectrumReference = max(abs(Xlos));
XlosDb = 20*log10(abs(Xlos)/spectrumReference + eps);
XmpDb = 20*log10(abs(Xmp)/spectrumReference + eps);

timeReference = mean(abs(xrLos));
fadingDb = 20*log10( ...
    sqrt(mean(abs(xrMultipath).^2))/sqrt(mean(abs(xrLos).^2)));

%% Results
fig = figure(Color='w',Name='FMCW Received-Signal Multipath Fading');
tl = tiledlayout(fig,1,2,TileSpacing='compact',Padding='compact');

ax1 = nexttile(tl);
plot(ax1,t*1e6,abs(xrLos)/timeReference,'--',LineWidth=1.2)
hold(ax1,'on')
plot(ax1,t*1e6,abs(xrMultipath)/timeReference,LineWidth=1.2)
hold(ax1,'off')
grid(ax1,'on')
xlabel(ax1,'Time within one chirp (\mus)')
ylabel(ax1,'Amplitude (LOS reference)')
title(ax1,sprintf('Time-domain fading: %.2f dB',fadingDb))
legend(ax1,'LOS','Multipath',Location='best')

ax2 = nexttile(tl);
plot(ax2,fBaseband/1e6,XlosDb,'--',LineWidth=1.2)
hold(ax2,'on')
plot(ax2,fBaseband/1e6,XmpDb,LineWidth=1.2)
hold(ax2,'off')
grid(ax2,'on')
xlabel(ax2,'Baseband frequency (MHz)')
ylabel(ax2,'Magnitude (dB, LOS reference)')
title(ax2,'Frequency-domain fading')
legend(ax2,'LOS','Multipath',Location='best')
ylim(ax2,[-50 5])

fprintf('Range resolution: %.3f m\n',c/(2*B))
fprintf('Excess equivalent ranges: %.6f m, %.6f m\n', ...
    targetRange(2)-targetRange(1),targetRange(3)-targetRange(1))
fprintf('Received-signal RMS fading: %.2f dB\n',fadingDb)
