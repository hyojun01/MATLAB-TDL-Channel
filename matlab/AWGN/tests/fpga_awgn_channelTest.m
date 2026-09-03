classdef fpga_awgn_channelTest < matlab.unittest.TestCase
    % Public behavior: scaling, bypass, streaming, calibration and saturation.
    methods (TestClassSetup)
        function addSourceFolder(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function absoluteNoisePowerAndIQStatistics(testCase)
            cfg = struct('fixedPoint', false);
            stream = RandStream('mt19937ar', 'Seed', 101);
            [y, info] = fpga_awgn_channel(complex(zeros(2^18,1)), cfg, stream);
            measuredBandDbm = 10*log10(mean(abs(y).^2)*50e6/61.44e6);
            iqCorrelation = corrcoef(real(y), imag(y));
            lagOne = mean(y(2:end).*conj(y(1:end-1)))/mean(abs(y).^2);
            testCase.verifyEqual(measuredBandDbm, -64, 'AbsTol', 0.06);
            testCase.verifyEqual(info.noiseDensityDbmHz, -140.98970004336, 'AbsTol', 1e-9);
            testCase.verifyEqual(var(real(y),1)/info.requestedVariance, 0.5, 'AbsTol', 0.01);
            testCase.verifyEqual(var(imag(y),1)/info.requestedVariance, 0.5, 'AbsTol', 0.01);
            testCase.verifyLessThan(abs(mean(y))/info.sigmaIQ, 0.015);
            testCase.verifyLessThan(abs(iqCorrelation(1,2)), 0.015);
            testCase.verifyLessThan(abs(lagOne), 0.015);
        end

        function fixedPointDefaultPower(testCase)
            stream = RandStream('mt19937ar', 'Seed', 102);
            [y, info] = fpga_awgn_channel(complex(zeros(2^18,1)), struct, stream);
            measuredBandDbm = 10*log10(mean(abs(y).^2)*50e6/61.44e6);
            testCase.verifyEqual(measuredBandDbm, -64, 'AbsTol', 0.08);
            testCase.verifyTrue(info.noiseLevelOK);
            testCase.verifyEqual(info.outputClippedSamples, 0);
            testCase.verifyEqual(real(y)*2^15, round(real(y)*2^15), 'AbsTol', 0);
        end

        function offIsExactPortBypassAndDoesNotAdvanceRng(testCase)
            stream = RandStream('mt19937ar', 'Seed', 103);
            before = stream.State;
            x = complex([-12000, 0, 23000], [17000, -500, 2])/2^15;
            [y, info, n, xPort] = fpga_awgn_channel(x, struct('enable',false), stream);
            testCase.verifyEqual(y, x, 'AbsTol', 0);
            testCase.verifyEqual(xPort, x, 'AbsTol', 0);
            testCase.verifyEqual(n, complex(zeros(size(x))), 'AbsTol', 0);
            testCase.verifyEqual(stream.State, before);
            testCase.verifyEqual(info.status, 'AWGN OFF');
        end

        function unevenFramesMatchOneContinuousStream(testCase)
            wholeStream = RandStream('mt19937ar', 'Seed', 104, 'NormalTransform','Inversion');
            chunkStream = RandStream('mt19937ar', 'Seed', 104, 'NormalTransform','Inversion');
            x = complex(linspace(-0.1,0.1,2003)', 0.05*ones(2003,1));
            whole = fpga_awgn_channel(x, struct, wholeStream);
            first = fpga_awgn_channel(x(1:513), struct, chunkStream);
            second = fpga_awgn_channel(x(514:end), struct, chunkStream);
            testCase.verifyEqual([first; second], whole, 'AbsTol', 0);
        end

        function noiseDoesNotFollowInputAmplitude(testCase)
            streamA = RandStream('mt19937ar', 'Seed', 105);
            streamB = RandStream('mt19937ar', 'Seed', 105);
            x = complex(0.05*ones(4096,1), 0.02*ones(4096,1));
            [~, ~, nA] = fpga_awgn_channel(x, struct, streamA);
            [~, ~, nB] = fpga_awgn_channel(0.01*x, struct, streamB);
            testCase.verifyEqual(nA, nB, 'AbsTol', 0);
        end

        function calibrationAndBandwidthChangeVariance(testCase)
            cfg = struct('fixedPoint',false, 'bandwidthHz',25e6, 'powerAtUnitRmsDbm',10);
            stream = RandStream('mt19937ar', 'Seed', 106);
            [~, info] = fpga_awgn_channel(complex(zeros(8,1)), cfg, stream);
            expected = 10^(-74/10)*61.44e6/25e6;
            testCase.verifyEqual(info.requestedVariance, expected, 'AbsTol', expected*1e-12);
        end

        function highLevelSaturatesRatherThanWraps(testCase)
            stream = RandStream('mt19937ar', 'Seed', 107);
            x = complex(0.8*ones(2^15,1), 0.7*ones(2^15,1));
            [y, info] = fpga_awgn_channel(x, struct('noisePowerDbm',-3), stream);
            testCase.verifyGreaterThan(info.outputClippedSamples, 0);
            testCase.verifyLessThanOrEqual(max([real(y);imag(y)]), 1-2^-15);
            testCase.verifyGreaterThanOrEqual(min([real(y);imag(y)]), -1);
            testCase.verifyFalse(info.noiseLevelOK);
            testCase.verifyEqual(info.status, 'HEADROOM WARNING');
        end

        function tooSmallGainIsReported(testCase)
            stream = RandStream('mt19937ar', 'Seed', 108);
            [y, info] = fpga_awgn_channel(complex(zeros(1024,1)), ...
                struct('noisePowerDbm',-180), stream);
            testCase.verifyEqual(info.gainCode, 0);
            testCase.verifyEqual(y, complex(zeros(1024,1)), 'AbsTol', 0);
            testCase.verifyEqual(info.status, 'BELOW RESOLUTION GUIDELINE');
        end

        function gainOverflowIsReported(testCase)
            stream = RandStream('mt19937ar', 'Seed', 109);
            [~, info] = fpga_awgn_channel(complex(zeros(1024,1)), ...
                struct('noisePowerDbm',10), stream);
            testCase.verifyTrue(info.gainClipped);
            testCase.verifyEqual(info.gainCode, 2^24-1);
            testCase.verifyEqual(info.status, 'GAIN SATURATED');
        end

        function impossibleBandwidthIsRejected(testCase)
            stream = RandStream('mt19937ar', 'Seed', 110);
            testCase.verifyError(@() fpga_awgn_channel(complex(zeros(8,1)), ...
                struct('bandwidthHz',70e6), stream), 'fpga_awgn:Bandwidth');
        end
    end
end
