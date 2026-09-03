%% FPGA AWGN simulation: 61.44 MS/s complex I/Q, 50 MHz measurement bandwidth
% Run this script from the matlab folder. See README_AWGN.md for Korean notes.
% Required: Signal Processing Toolbox and Communications Toolbox for the
% demo waveform. fpga_awgn_channel itself requires only base MATLAB.
% No clear/close all/rng calls: unrelated workspace, figures and RNG survive.

%% 1. Host settings -- edit here
cfg = struct;
cfg.sampleRateHz = 61.44e6;
cfg.bandwidthHz = 50e6;              % -25 MHz ... +25 MHz, NOT +/-50 MHz
cfg.noisePowerDbm = -64;             % TOTAL noise inside the 50 MHz band
cfg.powerAtUnitRmsDbm = 0;           % EXAMPLE calibration; replace for hardware
cfg.enable = true;
cfg.fixedPoint = true;
cfg.iqWordLength = 16;               % signed, 1 sign bit + 15 fractional bits
cfg.iqFractionLength = 15;
cfg.gainWordLength = 24;             % unsigned, 24 fractional bits
cfg.gainFractionLength = 24;
cfg.peakSigma = 6;

numSamples = 2^18;                  % 4.267 ms
frameLength = 4096;                 % processing chunk, NOT FPGA clock rate
demoSignalPowerDbm = -34;           % about 30 dB SNR in 50 MHz
seed = 20260827;
makePlots = true;
saveResults = true;                 % writes only matlab/results/awgn_*

% To use your own stream, replace [] with a normalized complex vector.
% Example for signed 16-bit integer I/Q codes:
% externalIQ = complex(double(iCodes), double(qCodes))/2^15;
% Example for a MAT file: loadedIQ = load('capture.mat','iq');
% externalIQ = loadedIQ.iq;
% External data is NOT rescaled or resampled; it must already be 61.44 MS/s.
externalIQ = [];

%% 2. Input stream -- pulse-shaped 16-QAM is just a test stimulus
if isempty(externalIQ)
    signalStream = RandStream('mt19937ar', 'Seed', seed);
    samplesPerSymbol = 2;
    symbolRateHz = cfg.sampleRateHz/samplesPerSymbol;  % 30.72 Msymbol/s
    rolloff = cfg.bandwidthHz/symbolRateHz-1;         % 0.627604...
    assert(rolloff >= 0 && rolloff <= 1, ...
        'Demo RRC rolloff must be in [0,1]; use externalIQ for other bands.');
    filterSpan = 12;
    pulse = rcosdesign(rolloff, filterSpan, samplesPerSymbol, 'sqrt');
    symbolCount = ceil(numSamples/samplesPerSymbol)+2*filterSpan;
    symbols = qammod(randi(signalStream, [0 15], symbolCount, 1), ...
        16, 'UnitAveragePower', true);
    shaped = upfirdn(symbols, pulse, samplesPerSymbol, 1);
    % Remove startup/tail effects using guard symbols on both ends.
    first = filterSpan*samplesPerSymbol + (numel(pulse)-1)/2 + 1;
    xInput = shaped(first:first+numSamples-1);
    desiredPower = 10^((demoSignalPowerDbm-cfg.powerAtUnitRmsDbm)/10);
    xInput = xInput*sqrt(desiredPower/mean(abs(xInput).^2));
else
    validateattributes(externalIQ, {'double', 'single'}, ...
        {'vector', 'nonempty', 'finite'});
    xInput = double(externalIQ(:));
    numSamples = numel(xInput);
end
assert(numSamples >= 8192, 'Use at least 8192 samples for PSD measurement.');

%% 3. FPGA end-of-path injection, one frame at a time
% Keep this object alive through every frame, including host level changes.
noiseStream = RandStream('mt19937ar', 'Seed', seed+1, ...
    'NormalTransform', 'Inversion');
yFpga = complex(zeros(numSamples, 1));
xPort = complex(zeros(numSamples, 1));
noiseIdeal = complex(zeros(numSamples, 1));
numFrames = ceil(numSamples/frameLength);
frameInfo = cell(numFrames, 1);
for frame = 1:numFrames
    idx = (frame-1)*frameLength+1:min(frame*frameLength, numSamples);
    [yFpga(idx), frameInfo{frame}, noiseIdeal(idx), xPort(idx)] = ...
        fpga_awgn_channel(xInput(idx), cfg, noiseStream);
end
frameInfo = vertcat(frameInfo{:});
host = frameInfo(1);
yIdeal = xPort+noiseIdeal;
noiseEffective = yFpga-xPort;        % Includes noise rounding and clipping
implementationError = yFpga-yIdeal;

%% 4. Independent MATLAB awgn reference, explicit signal-power assumption
% SNR here is a direct POWER RATIO over Fs, not Eb/No or Es/No.
% Never use awgn(...,'measured'): hardware noise must not track input power.
signalPower = mean(abs(xPort).^2);
awgnReference = struct('available', false);
if cfg.enable && signalPower > 0 && exist('awgn', 'file') == 2
    sampleSnrDb = 10*log10(signalPower/host.requestedVariance);
    referenceStream = RandStream('mt19937ar', 'Seed', seed+2);
    [referenceOutput, referenceVariance] = awgn(xPort, sampleSnrDb, ...
        10*log10(signalPower), referenceStream);
    awgnReference.available = true;
    awgnReference.variance = referenceVariance;
    awgnReference.measuredVariance = mean(abs(referenceOutput-xPort).^2);
    assert(abs(referenceVariance/host.requestedVariance-1) < 1e-10, ...
        'The explicit awgn reference variance disagrees with the host value.');
end

%% 5. Two-sided PSD and power measurements
fftLength = 8192;
window = hann(fftLength, 'periodic');
[psd, frequencyHz] = pwelch([xPort, yFpga, noiseIdeal, noiseEffective], ...
    window, fftLength/2, fftLength, cfg.sampleRateHz, 'centered', 'psd');
binHz = cfg.sampleRateHz/fftLength;
% Partial-bin integration; wrap the Nyquist bin when bandwidth equals Fs.
bandBinHz = zeros(size(frequencyHz));
for shift = [-cfg.sampleRateHz, 0, cfg.sampleRateHz]
    centers = frequencyHz+shift;
    bandBinHz = bandBinHz + max(0, ...
        min(centers+binHz/2, cfg.bandwidthHz/2) ...
        - max(centers-binHz/2, -cfg.bandwidthHz/2));
end
bandPower = sum(psd.*bandBinHz, 1);
fullNoisePower = [mean(abs(noiseIdeal).^2), mean(abs(noiseEffective).^2)];
toDbm = @(power) cfg.powerAtUnitRmsDbm+10*log10(power);
requestedBandDbm = cfg.noisePowerDbm;
requestedFullDbm = host.totalNoisePowerDbm;
requestedDensity = host.noiseDensityDbmHz;
if ~cfg.enable
    requestedBandDbm = -Inf;
    requestedFullDbm = -Inf;
    requestedDensity = -Inf;
end
noiseSummary = table( ...
    [requestedBandDbm; toDbm(bandPower(3)); toDbm(bandPower(4))], ...
    [requestedFullDbm; toDbm(fullNoisePower(1)); toDbm(fullNoisePower(2))], ...
    [requestedDensity; toDbm(bandPower(3)/cfg.bandwidthHz); ...
        toDbm(bandPower(4)/cfg.bandwidthHz)], ...
    'VariableNames', {'InBand_dBm', 'FullSampleBand_dBm', 'Density_dBm_Hz'}, ...
    'RowNames', {'Requested', 'IdealGaussian', 'FPGA_effective'});
snrInBandDb = 10*log10(bandPower(1)/bandPower(4));
snrSampleDb = 10*log10(signalPower/fullNoisePower(2));
iqVariance = [var(real(noiseIdeal), 1), var(imag(noiseIdeal), 1)];
iqCorrelation = NaN;
noiseAcf = zeros(33, 1);
if cfg.enable
    correlationMatrix = corrcoef(real(noiseIdeal), imag(noiseIdeal));
    iqCorrelation = correlationMatrix(1,2);
    noiseAcf = xcorr(noiseIdeal, 32, 'biased');
    noiseAcf = noiseAcf(33:end)/noiseAcf(33);
end

fprintf('\n--- FPGA AWGN / complex I/Q ---\n');
fprintf('Fs = %.2f MS/s, measurement band = %.2f MHz (+/- %.2f MHz)\n', ...
    cfg.sampleRateHz/1e6, cfg.bandwidthHz/1e6, cfg.bandwidthHz/2e6);
fprintf('Calibration: mean(|x|^2)=1 -> %.2f dBm (replace for your hardware)\n', ...
    cfg.powerAtUnitRmsDbm);
fprintf('Full-sample-band correction: %.6f dB\n', ...
    10*log10(cfg.sampleRateHz/cfg.bandwidthHz));
fprintf('Requested complex variance = %.9g; sigma per I or Q = %.9g\n', ...
    host.requestedVariance, host.sigmaIQ);
fprintf('Sigma register = %u; applied sigma = %.9g\n', ...
    host.gainCode, host.appliedSigmaIQ);
fprintf('Status: %s\n', strjoin(unique(string({frameInfo.status})), ', '));
fprintf('Noise quality/headroom guideline for this input: %.2f to %.2f dBm\n', ...
    host.minRecommendedNoiseDbm, min([frameInfo.maxRecommendedNoiseDbm]));
fprintf('Clipped samples: input=%d, noise=%d, output=%d\n', ...
    sum([frameInfo.inputClippedSamples]), sum([frameInfo.noiseClippedSamples]), ...
    sum([frameInfo.outputClippedSamples]));
disp(noiseSummary);
fprintf('Measured SNR: %.3f dB in %.2f MHz; %.3f dB over Fs\n', ...
    snrInBandDb, cfg.bandwidthHz/1e6, snrSampleDb);
fprintf('Ideal I/Q variances: %.9g / %.9g; I-Q correlation: %.5f\n', ...
    iqVariance, iqCorrelation);
fprintf('RMS implementation error = %.9g (normalized I/Q)\n', ...
    sqrt(mean(abs(implementationError).^2)));
if awgnReference.available
    fprintf('awgn reference: variance %.9g, measured %.9g\n', ...
        awgnReference.variance, awgnReference.measuredVariance);
end

%% 6. Diagnostic figures (raw sample statistics, not a demodulator/BER test)
if makePlots
    awgnFigure = figure('Name', 'FPGA AWGN verification', 'Color', 'w', ...
        'Position', [80 80 1280 820]);
    % Four results in a 2-by-2 subplot figure.
    sgtitle(awgnFigure, sprintf('Complex AWGN | %.2f MS/s | %.0f MHz | %.0f dBm in band', ...
        cfg.sampleRateHz/1e6, cfg.bandwidthHz/1e6, requestedBandDbm));

    ax = subplot(2, 2, 1, 'Parent', awgnFigure);
    if cfg.enable
        plot(ax, frequencyHz/1e6, toDbm(max(psd, realmin)), 'LineWidth', 1);
    else
        plot(ax, frequencyHz/1e6, toDbm(max(psd(:,1:2), realmin)), 'LineWidth', 1);
    end
    hold(ax, 'on');
    xline(ax, -cfg.bandwidthHz/2e6, ':', 'HandleVisibility', 'off');
    xline(ax, cfg.bandwidthHz/2e6, ':', 'HandleVisibility', 'off');
    if cfg.enable
        yline(ax, requestedDensity, '--k', 'Target density', ...
            'HandleVisibility', 'off');
    end
    xlabel(ax, 'Baseband frequency (MHz)'); ylabel(ax, 'PSD (dBm/Hz)');
    title(ax, 'Two-sided power spectral density'); grid(ax, 'on');
    if cfg.enable
        legend(ax, 'Input at FPGA port', 'FPGA output', 'Ideal noise', ...
            'Effective injected noise', 'Location', 'southwest');
    else
        legend(ax, 'Input at FPGA port', 'FPGA output', 'Location', 'southwest');
    end
    xlim(ax, [-cfg.sampleRateHz cfg.sampleRateHz]/2e6);
    if cfg.enable
        ylim(ax, [requestedDensity-15, max(toDbm(max(psd(:), realmin)))+5]);
    end

    ax = subplot(2, 2, 2, 'Parent', awgnFigure);
    preview = 1:min(256, numSamples);
    plot(ax, (preview-1)/cfg.sampleRateHz*1e6, ...
        [real(noiseIdeal(preview)), real(noiseEffective(preview))]);
    xlabel(ax, 'Time (microseconds)'); ylabel(ax, 'I amplitude (normalized)');
    title(ax, 'Injected I noise: floating point vs fixed point');
    legend(ax, 'Ideal', 'FPGA effective'); grid(ax, 'on');

    ax = subplot(2, 2, 3, 'Parent', awgnFigure);
    if cfg.enable
        edges = linspace(-5, 5, 81);
        if cfg.fixedPoint
            % Center bins on code groups; arbitrary edges can create false
            % histogram spikes when quantized samples fall on boundaries.
            normalizedLsb = 2^(-cfg.iqFractionLength)/host.sigmaIQ;
            codesPerBin = max(1, 2*floor(0.125/normalizedLsb/2)+1);
            binWidth = codesPerBin*normalizedLsb;
            halfBinCount = max(1, ceil(5/binWidth));
            edges = (-halfBinCount-0.5:halfBinCount+0.5)*binWidth;
        end
        histogram(ax, real(noiseEffective)/host.sigmaIQ, edges, ...
            'Normalization', 'pdf', 'FaceAlpha', 0.45);
        hold(ax, 'on');
        histogram(ax, imag(noiseEffective)/host.sigmaIQ, edges, ...
            'Normalization', 'pdf', 'FaceAlpha', 0.45);
        gridPoints = linspace(-5, 5, 400);
        plot(ax, gridPoints, exp(-gridPoints.^2/2)/sqrt(2*pi), ...
            'k', 'LineWidth', 1.5);
        legend(ax, 'Effective I', 'Effective Q', 'N(0,1)');
    else
        text(ax, 0.15, 0.5, 'AWGN OFF: no injected noise');
    end
    xlabel(ax, 'Noise / requested sigma'); ylabel(ax, 'Probability density');
    title(ax, 'Gaussian distribution and quantization'); grid(ax, 'on');

    ax = subplot(2, 2, 4, 'Parent', awgnFigure);
    stem(ax, 1:32, real(noiseAcf(2:end)), 'filled');
    hold(ax, 'on');
    yline(ax, 5/sqrt(numSamples), ':', '5/sqrt(N)', 'HandleVisibility', 'off');
    yline(ax, -5/sqrt(numSamples), ':', 'HandleVisibility', 'off');
    xlabel(ax, 'Lag (complex samples)'); ylabel(ax, 'Re\{R_n[k]/R_n[0]\}');
    title(ax, 'Ideal noise whiteness (zero lag excluded)'); grid(ax, 'on');
    diagnosticAxes = findall(awgnFigure, 'Type', 'axes');
    for axesIndex = 1:numel(diagnosticAxes)
        diagnosticAxes(axesIndex).Toolbar.Visible = 'on';
    end
end

%% 7. Results for subsequent comparison with FPGA captures
awgnResult = struct('cfg', cfg, 'noiseSummary', noiseSummary, ...
    'frameInfo', frameInfo, 'snrInBandDb', snrInBandDb, ...
    'snrSampleDb', snrSampleDb, 'iqVariance', iqVariance, ...
    'iqCorrelation', iqCorrelation, 'awgnReference', awgnReference, ...
    'frequencyHz', frequencyHz, 'psd', psd, ...
    'xInput', xInput, 'xPort', xPort, 'yIdeal', yIdeal, 'yFpga', yFpga, ...
    'noiseIdeal', noiseIdeal, 'noiseEffective', noiseEffective);
if saveResults
    resultFolder = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(resultFolder)
        mkdir(resultFolder);
    end
    save(fullfile(resultFolder, 'awgn_result.mat'), 'awgnResult');
    writetable(noiseSummary, fullfile(resultFolder, 'awgn_power_summary.csv'), ...
        'WriteRowNames', true);
    if makePlots
        exportFigureWithoutToolbars(awgnFigure, ...
            fullfile(resultFolder, 'awgn_verification.png'));
    end
    fprintf('Saved results: %s\n', resultFolder);
end


function exportFigureWithoutToolbars(figureHandle, outputPath)
% Hide controls only in the exported image; restore them even on errors.
exportAxes = findall(figureHandle, 'Type', 'axes');
exportToolbars = [exportAxes.Toolbar];
originalVisibility = get(exportToolbars, {'Visible'});
restoreToolbars = onCleanup(@() set(exportToolbars, ...
    {'Visible'}, originalVisibility));
set(exportToolbars, 'Visible', 'off');
exportgraphics(figureHandle, outputPath, 'Resolution', 150);
end
