function [y, info, noiseIdeal, xPort] = fpga_awgn_channel(x, cfg, noiseStream)
%FPGA_AWGN_CHANNEL Complex AWGN reference with optional fixed-point effects.
%   [Y,INFO,NIDEAL,XPORT] = FPGA_AWGN_CHANNEL(X,CFG,STREAM) adds noise
%   independent of X's power. X is a normalized complex I/Q vector. Reuse
%   the same RandStream object across frames; do not reseed each frame.
%
%   CFG fields (omitted fields use the defaults below):
%     sampleRateHz          61.44e6  Complex samples per second.
%     bandwidthHz           50e6     TOTAL two-sided measurement bandwidth.
%     noisePowerDbm         -64      Noise power inside bandwidthHz.
%     powerAtUnitRmsDbm     0        Calibration: mean(abs(x).^2)=1 -> dBm.
%     enable                true     OFF bypasses the noise path.
%     fixedPoint            true     Quantize input, gain, noise and output.
%     iqWordLength          16       Signed I and Q word length.
%     iqFractionLength      15       I and Q fractional bits.
%     gainWordLength        24       Unsigned sigma register word length.
%     gainFractionLength    24       Sigma register fractional bits.
%     peakSigma             6        Headroom planning margin, not a bound.
%
%   For complex white noise across the full sampled bandwidth Fs:
%     Pn = 10^((noisePowerDbm-powerAtUnitRmsDbm)/10) * Fs/B;
%     sigmaIQ = sqrt(Pn/2); n = sigmaIQ*(randnI + 1j*randnQ).
%   There is NO additional factor of 2 in complex two-sided PSD integration.
%   powerAtUnitRmsDbm absorbs RF gain, impedance and I/Q voltage conventions.
%
%   Fixed-point behavior: round to nearest (ties away from zero), saturate;
%   signed noise samples and input are added with an extended accumulator,
%   then saturated to the I/Q output width. Returned samples are doubles
%   representing fixed-point values, not MATLAB fi objects.
%
%   This is NOT synthesizable RTL or a bit-exact Gaussian PRNG model.
%   randn remains floating point; pipeline latency/valid signals are absent.
%   The 6-sigma limit is a conservative headroom estimate for the given
%   input frame, NOT an absolute limit or a guarantee of zero overflow.
%
%   See also simulate_fpga_awgn, RandStream, awgn.

defaults = struct('sampleRateHz', 61.44e6, 'bandwidthHz', 50e6, ...
    'noisePowerDbm', -64, 'powerAtUnitRmsDbm', 0, ...
    'enable', true, 'fixedPoint', true, ...
    'iqWordLength', 16, 'iqFractionLength', 15, ...
    'gainWordLength', 24, 'gainFractionLength', 24, 'peakSigma', 6);
validateattributes(cfg, {'struct'}, {'scalar'}, mfilename, 'cfg');
names = fieldnames(cfg);
for k = 1:numel(names)
    assert(isfield(defaults, names{k}), 'fpga_awgn:UnknownSetting', ...
        'Unknown configuration field: %s.', names{k});
    defaults.(names{k}) = cfg.(names{k});
end
cfg = defaults;
validateattributes(x, {'double', 'single'}, {'vector', 'nonempty', 'finite'}, ...
    mfilename, 'x');
assert(isa(noiseStream, 'RandStream'), 'fpga_awgn:InvalidStream', ...
    'Pass a RandStream object and retain it across frames.');
validateattributes(cfg.sampleRateHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.bandwidthHz, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
assert(cfg.bandwidthHz <= cfg.sampleRateHz, 'fpga_awgn:Bandwidth', ...
    'Complex two-sided bandwidth must not exceed sampleRateHz.');
validateattributes(cfg.noisePowerDbm, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(cfg.powerAtUnitRmsDbm, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(cfg.enable, {'logical'}, {'scalar'});
validateattributes(cfg.fixedPoint, {'logical'}, {'scalar'});
validateattributes(cfg.iqWordLength, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2, '<=', 32});
validateattributes(cfg.iqFractionLength, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0, '<', cfg.iqWordLength});
validateattributes(cfg.gainWordLength, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, '<=', 32});
validateattributes(cfg.gainFractionLength, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0, '<=', cfg.gainWordLength});
validateattributes(cfg.peakSigma, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

originalSize = size(x);
x = double(x(:));
variance = 10^((cfg.noisePowerDbm-cfg.powerAtUnitRmsDbm)/10) ...
    * cfg.sampleRateHz/cfg.bandwidthHz;
assert(isfinite(variance) && variance > 0, 'fpga_awgn:PowerRange', ...
    'Requested noise power cannot be represented in double precision.');
sigmaIQ = sqrt(variance/2);
iqLsb = 2^(-cfg.iqFractionLength);
gainLsb = 2^(-cfg.gainFractionLength);
gainMaxCode = 2^cfg.gainWordLength-1;
gainCodeUnclipped = round(sigmaIQ/gainLsb);
gainCode = min(gainCodeUnclipped, gainMaxCode);

info = struct;
info.cfg = cfg;
info.noiseDensityDbmHz = cfg.noisePowerDbm-10*log10(cfg.bandwidthHz);
info.totalNoisePowerDbm = cfg.noisePowerDbm ...
    + 10*log10(cfg.sampleRateHz/cfg.bandwidthHz);
info.requestedVariance = variance;      % E[abs(n)^2], I + Q together
info.sigmaIQ = sigmaIQ;
info.gainCode = gainCode;
info.appliedSigmaIQ = sigmaIQ;
info.inputClippedSamples = 0;
info.noiseClippedSamples = 0;
info.outputClippedSamples = 0;
info.gainClipped = false;
info.minRecommendedNoiseDbm = -Inf;
info.maxRecommendedNoiseDbm = Inf;
info.status = 'NOISE LEVEL OK';

if cfg.fixedPoint
    [xPort, info.inputClippedSamples] = quantizeIQ(x, cfg);
    info.appliedSigmaIQ = gainCode*gainLsb;
    info.gainClipped = gainCodeUnclipped > gainMaxCode;
    lowerRail = -2^(cfg.iqWordLength-1)*iqLsb;
    upperRail = (2^(cfg.iqWordLength-1)-1)*iqLsb;
    components = [real(xPort); imag(xPort)];
    headroom = min(upperRail-max(components), min(components)-lowerRail);
    maxSigma = min(max(0, headroom)/cfg.peakSigma, gainMaxCode*gainLsb);
    % Four LSB rms is a chosen noise-quality guideline, not a hard limit.
    minSigma = max(4*iqLsb, gainLsb);
    info.minRecommendedNoiseDbm = cfg.powerAtUnitRmsDbm ...
        + 10*log10(2*minSigma^2*cfg.bandwidthHz/cfg.sampleRateHz);
    info.maxRecommendedNoiseDbm = cfg.powerAtUnitRmsDbm ...
        + 10*log10(2*maxSigma^2*cfg.bandwidthHz/cfg.sampleRateHz);
else
    xPort = x;
end

noiseIdeal = complex(zeros(size(x)));
y = xPort;
if cfg.enable
    % Interleave I,Q draws so frame partitioning does not change the sequence.
    % Do not normalize each frame's measured noise power to the target.
    gaussianIQ = randn(noiseStream, 2, numel(x));
    gaussian = complex(gaussianIQ(1,:).', gaussianIQ(2,:).');
    noiseIdeal = sigmaIQ*gaussian;
    if cfg.fixedPoint
        [noiseFixed, info.noiseClippedSamples] = quantizeIQ( ...
            info.appliedSigmaIQ*gaussian, cfg);
        [y, info.outputClippedSamples] = quantizeIQ(xPort+noiseFixed, cfg);
    else
        y = xPort+noiseIdeal;
    end
end

if ~cfg.enable
    info.status = 'AWGN OFF';
    info.appliedSigmaIQ = 0;
elseif info.gainClipped
    info.status = 'GAIN SATURATED';
elseif info.noiseClippedSamples > 0 || info.outputClippedSamples > 0 ...
        || cfg.noisePowerDbm > info.maxRecommendedNoiseDbm
    info.status = 'HEADROOM WARNING';
elseif cfg.noisePowerDbm < info.minRecommendedNoiseDbm
    info.status = 'BELOW RESOLUTION GUIDELINE';
end
if info.inputClippedSamples > 0
    info.status = [info.status ' / INPUT CLIPPED'];
end
info.noiseLevelOK = strcmp(info.status, 'NOISE LEVEL OK');
y = reshape(y, originalSize);
noiseIdeal = reshape(noiseIdeal, originalSize);
xPort = reshape(xPort, originalSize);
end

function [q, clippedSamples] = quantizeIQ(x, cfg)
scale = 2^cfg.iqFractionLength;
minCode = -2^(cfg.iqWordLength-1);
maxCode = 2^(cfg.iqWordLength-1)-1;
iCode = round(real(x)*scale);
qCode = round(imag(x)*scale);
clippedSamples = nnz(iCode < minCode | iCode > maxCode ...
    | qCode < minCode | qCode > maxCode);
q = complex(min(max(iCode, minCode), maxCode), ...
    min(max(qCode, minCode), maxCode))/scale;
end
