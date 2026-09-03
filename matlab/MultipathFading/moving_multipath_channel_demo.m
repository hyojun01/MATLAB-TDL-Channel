%% Static Multipath Channel with One Moving Path
% This script extends a static tapped-delay-line channel with one reflected
% path whose complex gain rotates at a specified Doppler frequency.
% Transmitter, receiver, and all other reflectors remain stationary.

clearvars;
close all;
clc;

%% User configuration
sampleRateHz = 5e6;
carrierFrequencyHz = 3.5e9;
simulationDurationSec = 25e-3;

% Static paths. The first entry is the LOS reference path.
staticPathDelaysSec = [0, 70e-9, 180e-9];
staticPathAttenuationsDb = [0, 3, 8];
staticReflectionPhasesDeg = [0, 0, 0];

% One path reflected by a moving object.
movingPathDelaySec = 280e-9;
movingPathAttenuationDb = 3;
movingReflectionPhaseDeg = 0;
movingPathDopplerHz = 200;

% A constant-envelope complex-baseband tone probes temporal fading.
% Set this value to zero to probe the RF carrier frequency itself.
basebandToneHz = 0;

%% Validate the configuration
validateattributes(sampleRateHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, "sampleRateHz");
validateattributes(carrierFrequencyHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, ...
    "carrierFrequencyHz");
validateattributes(simulationDurationSec, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, ...
    "simulationDurationSec");
validateattributes(staticPathDelaysSec, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, ...
    "staticPathDelaysSec");
validateattributes(staticPathAttenuationsDb, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, ...
    "staticPathAttenuationsDb");
validateattributes(staticReflectionPhasesDeg, {'numeric'}, ...
    {'vector', 'real', 'finite'}, mfilename, ...
    "staticReflectionPhasesDeg");
validateattributes(movingPathDelaySec, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, ...
    "movingPathDelaySec");
validateattributes(movingPathAttenuationDb, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, ...
    "movingPathAttenuationDb");
validateattributes(movingReflectionPhaseDeg, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, ...
    "movingReflectionPhaseDeg");
validateattributes(movingPathDopplerHz, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "movingPathDopplerHz");
validateattributes(basebandToneHz, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "basebandToneHz");

staticPathDelaysSec = staticPathDelaysSec(:).';
staticPathAttenuationsDb = staticPathAttenuationsDb(:).';
staticReflectionPhasesDeg = staticReflectionPhasesDeg(:).';
staticPathCount = numel(staticPathDelaysSec);

if numel(staticPathAttenuationsDb) ~= staticPathCount || ...
        numel(staticReflectionPhasesDeg) ~= staticPathCount
    error("ChannelConfig:StaticPathVectorSizeMismatch", ...
        "All static-path vectors must have the same number of elements.");
end

if abs(basebandToneHz) >= sampleRateHz/2
    error("ChannelConfig:ToneOutsideNyquist", ...
        "The baseband tone must lie strictly inside the Nyquist interval.");
end

%% Combine and sort the static and moving paths
pathDelaysSec = [staticPathDelaysSec, movingPathDelaySec];
pathAttenuationsDb = [staticPathAttenuationsDb, ...
    movingPathAttenuationDb];
pathReflectionPhasesDeg = [staticReflectionPhasesDeg, ...
    movingReflectionPhaseDeg];
pathDopplerHz = [zeros(1, staticPathCount), movingPathDopplerHz];
pathTypes = [repmat("Static", 1, staticPathCount), "Moving"];

[pathDelaysSec, sortOrder] = sort(pathDelaysSec);
pathAttenuationsDb = pathAttenuationsDb(sortOrder);
pathReflectionPhasesDeg = pathReflectionPhasesDeg(sortOrder);
pathDopplerHz = pathDopplerHz(sortOrder);
pathTypes = pathTypes(sortOrder);

movingPathIndex = find(pathTypes == "Moving", 1);
staticPathIndices = find(pathTypes == "Static");
pathCount = numel(pathDelaysSec);

%% Construct the time-varying path-gain matrix
sampleCount = floor(simulationDurationSec*sampleRateHz) + 1;
timeSec = (0:sampleCount-1).'/sampleRateHz;

pathMagnitudes = 10.^(-pathAttenuationsDb/20);
initialPropagationPhasesRad = -2*pi*carrierFrequencyHz*pathDelaysSec;
reflectionPhasesRad = deg2rad(pathReflectionPhasesDeg);
initialPathPhasesRad = initialPropagationPhasesRad + ...
    reflectionPhasesRad;
initialPathGains = pathMagnitudes.*exp(1j*initialPathPhasesRad);

% Static paths repeat the same complex coefficient at every time sample.
pathGains = repmat(initialPathGains, sampleCount, 1);

% Positive Doppler corresponds to counterclockwise complex-phase rotation.
movingPathGain = pathMagnitudes(movingPathIndex).*exp(1j*( ...
    initialPathPhasesRad(movingPathIndex) + ...
    2*pi*movingPathDopplerHz*timeSec));
pathGains(:, movingPathIndex) = movingPathGain;

channelConfiguration = table((1:pathCount).', pathTypes.', ...
    pathDelaysSec.'*1e9, pathAttenuationsDb.', pathDopplerHz.', ...
    VariableNames=["Path", "Type", "Delay_ns", "Attenuation_dB", ...
    "Doppler_Hz"]);

disp("Multipath configuration:");
disp(channelConfiguration);

%% Filter a constant-envelope baseband probe signal
transmittedSignal = exp(1j*2*pi*basebandToneHz*timeSec);

channel = comm.ChannelFilter( ...
    SampleRate=sampleRateHz, ...
    PathDelays=pathDelaysSec);
receivedSignal = channel(transmittedSignal, pathGains);

%% Calculate the theoretical narrowband channel coefficient
% The explicit channel delay contributes an additional phase at the
% baseband tone frequency. The carrier-frequency phase is already included
% in each complex path gain.
basebandDelayPhases = exp(-1j*2*pi*basebandToneHz*pathDelaysSec);
effectiveChannel = sum(pathGains.*basebandDelayPhases, 2);
predictedReceivedSignal = effectiveChannel.*transmittedSignal;

channelInformation = info(channel);
filterCoefficientCount = size( ...
    channelInformation.ChannelFilterCoefficients, 2);
settlingSampleCount = filterCoefficientCount + ...
    ceil(max(pathDelaysSec)*sampleRateHz) + 2;
analysisIndices = (settlingSampleCount:sampleCount).';

idealOutputError = max(abs(receivedSignal(analysisIndices) - ...
    predictedReceivedSignal(analysisIndices)));

% Account for the finite fractional-delay FIR implementation. Its DC gain
% can differ slightly from the ideal value of one for a fractional delay.
filterCoefficientSums = sum( ...
    channelInformation.ChannelFilterCoefficients, 2).';
implementedEffectiveChannel = sum(pathGains.*( ...
    basebandDelayPhases.*filterCoefficientSums), 2);
implementedPrediction = implementedEffectiveChannel.*transmittedSignal;
implementedOutputError = max(abs(receivedSignal(analysisIndices) - ...
    implementedPrediction(analysisIndices)));

staticGainVariation = max(abs(diff( ...
    pathGains(:, staticPathIndices), 1, 1)), [], "all");
movingMagnitudeVariation = max(abs(diff(abs(movingPathGain))));

%% Report the physical and numerical interpretation
speedOfLightMps = 299792458;
wavelengthM = speedOfLightMps/carrierFrequencyHz;
totalPathLengthRateMps = -wavelengthM*movingPathDopplerHz;

staticEffectiveChannel = sum(initialPathGains(staticPathIndices).* ...
    basebandDelayPhases(staticPathIndices));
movingEffectiveMagnitude = pathMagnitudes(movingPathIndex);
expectedMinimumMagnitude = abs(abs(staticEffectiveChannel) - ...
    movingEffectiveMagnitude);
expectedMaximumMagnitude = abs(staticEffectiveChannel) + ...
    movingEffectiveMagnitude;

fprintf("Static-path gain variation: %.3g\n", staticGainVariation);
fprintf("Moving-path magnitude variation: %.3g\n", ...
    movingMagnitudeVariation);
fprintf("Channel-filter/ideal-path maximum difference: %.3g\n", ...
    idealOutputError);
fprintf("Channel-filter/implemented-model maximum error: %.3g\n", ...
    implementedOutputError);
fprintf("Equivalent total reflected-path length rate: %.3f m/s\n", ...
    totalPathLengthRateMps);
fprintf("Expected channel-magnitude range: %.4f to %.4f\n", ...
    expectedMinimumMagnitude, expectedMaximumMagnitude);

if abs(movingPathDopplerHz) > eps(max(1, abs(movingPathDopplerHz)))
    fprintf("Doppler rotation period: %.3f ms\n", ...
        1e3/abs(movingPathDopplerHz));
else
    fprintf("Doppler rotation period: Inf (the moving path is static)\n");
end

%% Visualize the time-selective channel
magnitudeFloor = 10^(-80/20);
effectiveChannelDb = 20*log10(max(abs(effectiveChannel), ...
    magnitudeFloor));
receivedEnvelopeDb = 20*log10(max(abs(receivedSignal), ...
    magnitudeFloor));

%% Calculate the time-frequency channel response
frequencyPointCount = 801;
frequencyHz = linspace(-sampleRateHz/2, sampleRateHz/2, ...
    frequencyPointCount);

% Downsample only the display grid. The channel simulation keeps its full
% sample rate, while the heatmap uses a compact set of time snapshots.
timeFrequencyPointCount = min(501, numel(analysisIndices));
timeFrequencyIndices = unique(round(linspace(analysisIndices(1), ...
    analysisIndices(end), timeFrequencyPointCount))).';
delayFrequencyPhases = exp(-1j*2*pi*pathDelaysSec.'*frequencyHz);
timeFrequencyResponse = pathGains(timeFrequencyIndices, :)* ...
    delayFrequencyPhases;
timeFrequencyResponseDb = 20*log10(max(abs(timeFrequencyResponse), ...
    magnitudeFloor));

% Select four snapshots across one Doppler period when the simulation is
% long enough. Otherwise, spread the snapshots over the available time.
dopplerTolerance = eps(max(1, abs(movingPathDopplerHz)));
analysisStartTimeSec = timeSec(analysisIndices(1));
availableTimeSec = timeSec(end) - analysisStartTimeSec;
if abs(movingPathDopplerHz) > dopplerTolerance
    dopplerPeriodSec = 1/abs(movingPathDopplerHz);
    if availableTimeSec >= 0.75*dopplerPeriodSec
        snapshotTimesSec = analysisStartTimeSec + ...
            [0, 0.25, 0.50, 0.75]*dopplerPeriodSec;
    else
        snapshotTimesSec = linspace(analysisStartTimeSec, ...
            timeSec(end), 4);
    end
else
    snapshotTimesSec = linspace(analysisStartTimeSec, ...
        timeSec(end), 4);
end

snapshotIndices = round(snapshotTimesSec*sampleRateHz) + 1;
snapshotIndices = min(max(snapshotIndices, analysisIndices(1)), ...
    sampleCount);
snapshotFrequencyResponse = pathGains(snapshotIndices, :)* ...
    delayFrequencyPhases;
snapshotFrequencyResponseDb = 20*log10(max( ...
    abs(snapshotFrequencyResponse), magnitudeFloor));
snapshotLabels = compose("t = %.2f ms", ...
    timeSec(snapshotIndices)*1e3);

% Limit plotted point counts without changing the simulated data.
displayPointLimit = 6000;
displayStride = max(1, ceil(sampleCount/displayPointLimit));
timeDisplayIndices = (1:displayStride:sampleCount).';
analysisDisplayIndices = analysisIndices(1:displayStride:end);

%% Visualize the time- and frequency-selective channel
channelFigure = figure(Name="One Moving Multipath Component", ...
    Color="white");
channelLayout = tiledlayout(channelFigure, 3, 2, ...
    TileSpacing="compact", Padding="compact");
title(channelLayout, ...
    "Fixed Tx/Rx: Static Paths Plus One Doppler-Shifted Path");

pathAxes = nexttile(channelLayout);
stem(pathAxes, pathDelaysSec(staticPathIndices)*1e9, ...
    -pathAttenuationsDb(staticPathIndices), "filled", ...
    LineWidth=1.2, DisplayName="Static paths");
hold(pathAxes, "on");
stem(pathAxes, pathDelaysSec(movingPathIndex)*1e9, ...
    -pathAttenuationsDb(movingPathIndex), "filled", ...
    LineWidth=1.2, DisplayName="Moving path");
hold(pathAxes, "off");
grid(pathAxes, "on");
xlabel(pathAxes, "Excess delay (ns)");
ylabel(pathAxes, "Relative path gain (dB)");
title(pathAxes, "Path-delay configuration");
legend(pathAxes, Location="best");

movingGainAxes = nexttile(channelLayout);
plot(movingGainAxes, timeSec(timeDisplayIndices)*1e3, ...
    real(movingPathGain(timeDisplayIndices)), LineWidth=1.1, ...
    DisplayName="In-phase");
hold(movingGainAxes, "on");
plot(movingGainAxes, timeSec(timeDisplayIndices)*1e3, ...
    imag(movingPathGain(timeDisplayIndices)), LineWidth=1.1, ...
    DisplayName="Quadrature");
hold(movingGainAxes, "off");
grid(movingGainAxes, "on");
xlabel(movingGainAxes, "Time (ms)");
ylabel(movingGainAxes, "Moving-path gain");
title(movingGainAxes, sprintf("Moving tap: Doppler = %.1f Hz", ...
    movingPathDopplerHz));
legend(movingGainAxes, Location="best");

fadingAxes = nexttile(channelLayout);
plot(fadingAxes, timeSec(analysisDisplayIndices)*1e3, ...
    effectiveChannelDb(analysisDisplayIndices), LineWidth=1.4, ...
    DisplayName="Theoretical phasor sum");
hold(fadingAxes, "on");
plot(fadingAxes, timeSec(analysisDisplayIndices)*1e3, ...
    receivedEnvelopeDb(analysisDisplayIndices), "--", ...
    LineWidth=1.0, DisplayName="Channel-filter output");
hold(fadingAxes, "off");
grid(fadingAxes, "on");
xlabel(fadingAxes, "Time (ms)");
ylabel(fadingAxes, "Magnitude (dB)");
title(fadingAxes, "Temporal fading at the probe frequency");
legend(fadingAxes, Location="best");

phasorAxes = nexttile(channelLayout);
plot(phasorAxes, real(effectiveChannel(analysisDisplayIndices)), ...
    imag(effectiveChannel(analysisDisplayIndices)), LineWidth=1.2, ...
    DisplayName="Total channel");
hold(phasorAxes, "on");
plot(phasorAxes, real(staticEffectiveChannel), ...
    imag(staticEffectiveChannel), "kx", MarkerSize=10, ...
    LineWidth=2, DisplayName="Static-path sum");
plot(phasorAxes, real(effectiveChannel(analysisIndices(1))), ...
    imag(effectiveChannel(analysisIndices(1))), "o", ...
    MarkerSize=6, LineWidth=1.2, DisplayName="Initial point");
hold(phasorAxes, "off");
axis(phasorAxes, "equal");
grid(phasorAxes, "on");
xlabel(phasorAxes, "In-phase");
ylabel(phasorAxes, "Quadrature");
title(phasorAxes, "Rotating phasor around the static-path sum");
legend(phasorAxes, Location="best");

snapshotAxes = nexttile(channelLayout);
plot(snapshotAxes, frequencyHz/1e6, ...
    snapshotFrequencyResponseDb.', LineWidth=1.1);
grid(snapshotAxes, "on");
xlabel(snapshotAxes, "Baseband frequency (MHz)");
ylabel(snapshotAxes, "|H(t,f)| (dB)");
title(snapshotAxes, "Frequency-response snapshots");
legend(snapshotAxes, snapshotLabels, Location="best");

timeFrequencyAxes = nexttile(channelLayout);
imagesc(timeFrequencyAxes, frequencyHz/1e6, ...
    timeSec(timeFrequencyIndices)*1e3, timeFrequencyResponseDb);
axis(timeFrequencyAxes, "xy");
xlabel(timeFrequencyAxes, "Baseband frequency (MHz)");
ylabel(timeFrequencyAxes, "Time (ms)");
title(timeFrequencyAxes, "Time-frequency channel fading");
colormap(timeFrequencyAxes, parula(256));
channelColorbar = colorbar(timeFrequencyAxes);
channelColorbar.Label.String = "Magnitude (dB)";
colorLowerDb = max(-80, min(timeFrequencyResponseDb, [], "all"));
colorUpperDb = max(timeFrequencyResponseDb, [], "all");
if (colorUpperDb - colorLowerDb) < 1
    colorLowerDb = colorUpperDb - 1;
end
clim(timeFrequencyAxes, [colorLowerDb, colorUpperDb]);
