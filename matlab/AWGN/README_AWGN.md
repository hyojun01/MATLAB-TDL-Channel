# FPGA AWGN 채널 MATLAB 시뮬레이션

## 실행

MATLAB의 Current Folder를 이 `matlab` 폴더로 바꾸고 실행합니다.

```matlab
simulate_fpga_awgn
testResults = runtests('tests/fpga_awgn_channelTest.m');
assertSuccess(testResults);
```

`simulate_fpga_awgn.m`은 입력 생성 → 프레임 단위 AWGN 추가 → 전력/PSD 측정 → 그래프/결과 저장을 수행합니다. `fpga_awgn_channel.m`은 실제 I/Q 입력에도 재사용할 수 있는 채널 함수입니다. 데모와 PSD 측정에는 Communications Toolbox, Signal Processing Toolbox가 필요하고, 채널 함수 자체는 기본 MATLAB만 사용합니다.

## 전력 기준: 가장 먼저 확인할 부분

- 샘플레이트: **61.44 MS/s의 복소 샘플**. I와 Q가 각각 61.44 MS/s입니다.
- 측정 대역폭: **총 50 MHz**, 복소 기저대역에서 **−25~+25 MHz**입니다.
- `cfg.noisePowerDbm`: 이 50 MHz 안에 포함되는 I+Q 전체 잡음 전력입니다.
- `cfg.powerAtUnitRmsDbm`: `mean(abs(x).^2) = 1`인 디지털 복소 신호에 대응하는 물리 전력입니다. 기본 **0 dBm은 예시이며 장비 사양이 아닙니다.** 실제 DAC/RF 경로 이득과 임피던스 및 I/Q 전압 정의를 반영한 보정값으로 교체해야 합니다.

전력 보정값을 `P_ref`, 50 MHz 안의 목표 잡음 전력을 `P_B`라고 하면:

```text
noiseDensity_dBmHz = P_B - 10*log10(B)
noisePowerOverFs_dBm = P_B + 10*log10(Fs/B)
varianceComplex = 10^((P_B - P_ref)/10) * Fs/B
sigmaI = sigmaQ = sqrt(varianceComplex/2)
n[k] = sigmaI * (gI[k] + j*gQ[k])
y[k] = x[k] + n[k]
```

`gI`, `gQ`는 서로 독립적인 평균 0, 분산 1의 가우시안 난수입니다. 신호가 작아지거나 0이 되어도 잡음 설정은 바뀌지 않습니다. 프레임마다 잡음 전력을 강제로 정규화하지 않습니다.

기본값 −64 dBm / 50 MHz에서 이론값은 다음과 같습니다.

| 항목 | 값 |
|---|---:|
| 50 MHz 안의 잡음 전력 | −64 dBm |
| 잡음 밀도 | −140.9897 dBm/Hz |
| Fs/B 전력 보정 | 약 +0.894812 dB |
| 61.44 MHz 전체 잡음 전력 | 약 −63.105188 dBm |

첨부 그림의 잡음 밀도 숫자를 그대로 사용하지 않고, 사용자 조건인 50 MHz로 다시 계산합니다. 위 근삿값보다 스크립트 출력값을 우선합니다.

이 구현은 독립 잡음 샘플을 생성하므로 **잡음이 ±30.72 MHz 전체에 존재**합니다. 그중 ±25 MHz 안의 전력이 설정값이 되도록 분산을 조정합니다. 잡음 자체를 ±25 MHz로 제한하는 필터는 포함하지 않습니다. 그런 필터를 추가하면 시간 샘플이 상관되며, 실제 필터의 등가 잡음 대역폭(ENBW)으로 다시 계산해야 합니다. 입력 데모는 RRC 펄스 성형 16-QAM이며 이론 대역폭은 50 MHz입니다. 유한 길이 FIR이므로 대역 밖 누설이 완전히 0은 아닙니다.

## 수정할 설정

스크립트 맨 위의 값을 바꿉니다.

| 설정 | 의미 |
|---|---|
| `cfg.enable` | AWGN ON/OFF |
| `cfg.noisePowerDbm` | 50 MHz 안의 총 잡음 전력 |
| `cfg.powerAtUnitRmsDbm` | 정규화 I/Q와 물리 전력의 보정값 |
| `cfg.fixedPoint` | `false`: 이상적인 부동소수점, `true`: 양자화·포화 포함 |
| `cfg.iqWordLength`, `cfg.iqFractionLength` | 기본 signed 16-bit, fractional 15-bit |
| `cfg.gainWordLength`, `cfg.gainFractionLength` | 기본 unsigned 24-bit, fractional 24-bit인 sigma 계수 |
| `numSamples`, `frameLength`, `seed` | 관측 길이, 처리 프레임 길이, 재현 가능한 난수 초기값 |
| `externalIQ` | 자신의 복소 I/Q 벡터. 비우면 데모 자동 생성 |
| `saveResults`, `makePlots` | 결과 저장 및 그래프 표시 |

예를 들어 정수 I/Q 캡처를 사용하려면 `externalIQ = []` 줄을 다음과 같이 교체합니다.

```matlab
externalIQ = complex(double(iCodes), double(qCodes))/2^15;
```

입력은 이미 61.44 MS/s여야 합니다. 자동 리샘플링이나 전력 정규화는 하지 않습니다. PSD 측정을 위해 최소 8192개 샘플을 사용합니다. 채널 함수는 더 짧은 프레임과 행/열 벡터도 지원합니다.

## FPGA 모델에서 포함한 것과 제외한 것

```text
Host: P_B, Fs, B, P_ref -> sigma -> 양자화된 sigma 레지스터
FPGA 개념: Gaussian I/Q -> sigma 곱셈 -> 잡음 양자화 -> I/Q 덧셈 -> 포화 -> ON/OFF
```

입력 I/Q, 잡음 계수, 잡음 I/Q, 출력 I/Q의 양자화를 모사합니다. 반올림은 nearest/ties away from zero이고 overflow는 wrap이 아닌 saturation입니다. 덧셈 중간 결과에는 여유 비트를 가정합니다. AWGN OFF일 때 출력은 **양자화된 입력 포트 `xPort`와 정확히 동일**하며 난수 상태가 진행하지 않습니다. OFF에서도 입력 양자화 자체는 남습니다.

`NOISE LEVEL OK`는 선택한 수치 형식과 **관측한 입력 피크**에 따른 판정입니다. 상한은 기본 6σ headroom 및 계수 레지스터 범위로 계산합니다. 가우시안 분포에는 엄밀한 최대값이 없으므로 6σ는 무포화 보장이 아닙니다. 기본 4 LSB RMS 미만은 가우시안 분포 재현이 거칠어지는 영역으로 경고합니다. 이 하한은 설계 지침이며 물리적인 절대 한계가 아닙니다. 실제 입력 최대치, FPGA PRNG 출력 범위, DAC/RF 경로 한계가 정해지면 별도 상한 검증이 필요합니다.

난수는 MATLAB `randn` 기준 모델입니다. **HDL 합성용 코드나 특정 FPGA Gaussian PRNG의 비트 일치 모델이 아닙니다.** PRNG 주기/꼬리 분포, 내부 LUT·곱셈기 비트폭, 파이프라인 지연, AXI valid/ready, 클록과 병렬 샘플 수, DAC/RF 응답은 포함하지 않습니다. 구현 시 실제 생성기와 고정소수점 데이터 경로로 교체하고 이 결과와 통계적으로 비교해야 합니다. `frameLength`는 MATLAB 처리 단위이며 FPGA 클록을 뜻하지 않습니다.

연속 스트림 처리에서는 `RandStream`을 한 번 생성해 모든 프레임에 재사용합니다. 프레임마다 seed를 재설정하면 잡음 패턴이 반복됩니다.

## 결과 읽기

- `noiseSummary`: 목표 / 이상 가우시안 / FPGA 유효 잡음의 대역 내 전력, 전체 전력, 밀도.
- `noiseEffective = yFpga - xPort`: 양자화와 포화 영향을 포함한 실제 주입 성분. 입력 자체의 양자화 오차와 구분합니다.
- `awgnResult`: 입력, 출력, 잡음, 설정, PSD, 측정 통계. MATLAB workspace 및 `results/awgn_result.mat`에 저장합니다.
- `results/awgn_power_summary.csv`: 측정표.
- `results/awgn_verification.png`: PSD, 시간 파형, I/Q 분포, 자기상관 그래프.

결과 파일은 재실행 시 덮어씁니다. PSD는 복소 양측 스펙트럼을 Hz 단위로 적분하며, 단측 스펙트럼의 2배 보정을 추가하지 않습니다. 측정 오차는 유한 표본과 Welch 추정 때문에 발생합니다. 포화되거나 잡음이 LSB 아래이면 FPGA 유효 잡음은 이상적인 AWGN과 달라집니다. 이 예제는 복조기/BER 시험이 아니라 채널 통계 및 데이터 경로 검증입니다.

## 참고

- [MathWorks awgn: 명시적 신호 전력과 전체 잡음 분산](https://www.mathworks.com/help/comm/ref/awgn.html)
- [MathWorks pwelch: 복소 양측 PSD와 주파수 범위](https://www.mathworks.com/help/signal/ref/pwelch.html)
