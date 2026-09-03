%% Static Multipath Channel Demo
% This script models a deterministic, time-invariant multipath channel.
% Change the path-delay and attenuation vectors in the configuration section.
% The number of vector elements determines the number of paths.

clearvars;
close all;
clc;

%% User configuration
sampleRateHz = 200e6;
carrierFrequencyHz = 3e9;

% Path 1 is the LOS reference path. Other entries are excess delays.
% The number of entries determines the number of paths.
pathDelaysSec = [0, 20e-9, 40e-9];

% Positive values denote attenuation relative to the LOS path.
pathAttenuationsDb = [0, 4.4369749923, 20.9151498112];

% Optional phase added by reflection or scattering. Set all entries to zero
% if only the propagation phase caused by each path delay is required.
pathReflectionPhasesDeg = [0, 0, 0];

% Complex-baseband test pulse parameters.
signalDurationSec = 1.5e-6;
pulseCenterSec = 0.40e-6;
pulseStdSec = 40e-9;
basebandToneHz = 10e6;

%% Validate and organize the path configuration
validateattributes(sampleRateHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, "sampleRateHz");
validateattributes(carrierFrequencyHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, "carrierFrequencyHz");
validateattributes(pathDelaysSec, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, "pathDelaysSec");
validateattributes(pathAttenuationsDb, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, "pathAttenuationsDb");
validateattributes(pathReflectionPhasesDeg, {'numeric'}, ...
    {'vector', 'real', 'finite'}, mfilename, "pathReflectionPhasesDeg");

pathDelaysSec = pathDelaysSec(:).';
pathAttenuationsDb = pathAttenuationsDb(:).';
pathReflectionPhasesDeg = pathReflectionPhasesDeg(:).';
pathCount = numel(pathDelaysSec);

if numel(pathAttenuationsDb) ~= pathCount || ...
        numel(pathReflectionPhasesDeg) ~= pathCount
    error("ChannelConfig:PathVectorSizeMismatch", ...
        "pathDelaysSec, pathAttenuationsDb, and " + ...
        "pathReflectionPhasesDeg must have the same number of elements.");
end

if abs(basebandToneHz) >= sampleRateHz/2
    error("ChannelConfig:ToneOutsideNyquist", ...
        "The baseband tone must lie strictly inside the Nyquist interval.");
end

% Keep all path-related vectors in increasing delay order.
[pathDelaysSec, sortOrder] = sort(pathDelaysSec);
pathAttenuationsDb = pathAttenuationsDb(sortOrder);
pathReflectionPhasesDeg = pathReflectionPhasesDeg(sortOrder);

%% Construct fixed complex path gains
% A positive attenuation in dB becomes a voltage gain below one.
pathMagnitudes = 10.^(-pathAttenuationsDb/20);

% The explicit path delay shifts the baseband envelope. The same physical
% delay also creates a carrier-phase shift, which is stored in the complex
% path gain. Reflection/scattering can add another fixed phase.
propagationPhasesRad = -2*pi*carrierFrequencyHz*pathDelaysSec;
reflectionPhasesRad = deg2rad(pathReflectionPhasesDeg);
pathPhasesRad = propagationPhasesRad + reflectionPhasesRad;
complexPathGains = pathMagnitudes.*exp(1j*pathPhasesRad);

pathPhaseDeg = mod(rad2deg(angle(complexPathGains)) + 180, 360) - 180;
channelConfiguration = table((1:pathCount).', pathDelaysSec.'*1e9, ...
    pathAttenuationsDb.', pathMagnitudes.', pathPhaseDeg.', ...
    VariableNames=["Path", "Delay_ns", "Attenuation_dB", ...
    "Magnitude", "FixedPhase_deg"]);

disp("Static multipath configuration:");
disp(channelConfiguration);

%% Generate a complex-baseband test pulse
sampleCount = floor(signalDurationSec*sampleRateHz) + 1;
timeSec = (0:sampleCount-1).'/sampleRateHz;

gaussianEnvelope = exp(-0.5*((timeSec - pulseCenterSec)/pulseStdSec).^2);
transmittedSignal = gaussianEnvelope.* ...
    exp(1j*2*pi*basebandToneHz*(timeSec - pulseCenterSec));

%% Apply the static multipath channel
% Each row contains the path gains at one time sample. Repeating the same
% row for the entire signal makes the channel strictly time invariant.
staticPathGains = repmat(complexPathGains, sampleCount, 1);

channel = comm.ChannelFilter( ...
    SampleRate=sampleRateHz, ...
    PathDelays=pathDelaysSec);

receivedSignal = channel(transmittedSignal, staticPathGains);

pathGainVariation = max(abs(diff(staticPathGains, 1, 1)), [], "all");
fprintf("Maximum path-gain change over time: %.3g\n", pathGainVariation);

%% Calculate the theoretical static frequency response
frequencyPointCount = 8192;
frequencyHz = linspace(-sampleRateHz/2, sampleRateHz/2, ...
    frequencyPointCount).';

channelFrequencyResponse = ...
    exp(-1j*2*pi*frequencyHz*pathDelaysSec)*complexPathGains.';
magnitudeFloor = 10^(-80/20);
frequencyResponseDb = 20*log10(max(abs(channelFrequencyResponse), ...
    magnitudeFloor));

%% Obtain the implemented impulse response for visualization
impulseLeadSamples = 32;
impulseTailSamples = 64;
maximumDelaySamples = ceil(max(pathDelaysSec)*sampleRateHz);
impulseSampleCount = impulseLeadSamples + maximumDelaySamples + ...
    impulseTailSamples;

impulseInput = zeros(impulseSampleCount, 1);
impulseInput(impulseLeadSamples + 1) = 1;
impulsePathGains = repmat(complexPathGains, impulseSampleCount, 1);

impulseChannel = comm.ChannelFilter( ...
    SampleRate=sampleRateHz, ...
    PathDelays=pathDelaysSec);
implementedImpulseResponse = impulseChannel(impulseInput, impulsePathGains);
impulseTimeSec = ((0:impulseSampleCount-1).' - impulseLeadSamples)/sampleRateHz;

%% Visualize the static channel
figure(Name="Static Multipath Channel", Color="white");
channelLayout = tiledlayout(2, 2, TileSpacing="compact", Padding="compact");
title(channelLayout, "Static Multipath Channel: Fixed Delays and Path Gains");

nexttile;
stem(pathDelaysSec*1e9, -pathAttenuationsDb, "filled", LineWidth=1.2);
grid on;
xlabel("Excess delay (ns)");
ylabel("Relative path gain (dB)");
title("Power-delay profile settings");

nexttile;
plot(frequencyHz/1e6, frequencyResponseDb, LineWidth=1.3);
grid on;
xlabel("Baseband frequency (MHz)");
ylabel("|H(f)| (dB)");
title("Time-invariant frequency response");

nexttile;
plot(impulseTimeSec*1e9, abs(implementedImpulseResponse), LineWidth=1.2);
hold on;
stem(pathDelaysSec*1e9, pathMagnitudes, "filled", LineWidth=1.0);
hold off;
grid on;
xlabel("Delay (ns)");
ylabel("Magnitude");
title("Implemented channel impulse response");
legend("Fractional-delay filter response", "Specified paths", ...
    Location="best");

nexttile;
plot(timeSec*1e6, abs(transmittedSignal), LineWidth=1.2);
hold on;
plot(timeSec*1e6, abs(receivedSignal), LineWidth=1.2);
hold off;
grid on;
xlabel("Time (\mus)");
ylabel("Complex-baseband envelope");
title("Input and channel output");
legend("Transmitted", "Received", Location="best");
