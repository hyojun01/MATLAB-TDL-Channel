c = 299792458;
fc = 3.5e9;
lambda = c/fc;
BW = 200e6;

tiledlayout(2,1)

% frequency response
delayDifference = [50e-9, 100e-9, 150e-9];
amplDifference = [0.8, 0.3, 0.1];
f = linspace(-BW/2,BW/2,8001);
H = 1 + amplDifference(1)*exp(-1j*2*pi*f*delayDifference(1)) + ...
        amplDifference(2)*exp(-1j*2*pi*f*delayDifference(2)) + ...
        amplDifference(3)*exp(-1j*2*pi*f*delayDifference(3));

nexttile
plot(f/1e6,20*log10(abs(H)+eps),'LineWidth',1.3)
grid on
xlabel("Baseband frequency offset (MHz)")
ylabel("|H(f)| (dB)")
title("Frequency fading from two delayed paths")
ylim([-60 10])

% Spatial/phase fading
pathDifference = linspace(0,3*lambda,3000);
h = 1 + 0.9*exp(-1j*2*pi*pathDifference/lambda);

nexttile
plot(pathDifference/lambda,20*log10(abs(h)+eps),'LineWidth',1.3)
grid on
xlabel("Relative path-length change / wavelength")
ylabel("|h| (dB)")
title("Small-scale fading from changing relative phase")
ylim([-30 10])