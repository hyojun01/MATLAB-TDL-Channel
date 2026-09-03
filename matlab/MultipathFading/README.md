# Multipath fading 학습 스크립트

[MathWorks Multipath Fading Channel 예제](https://www.mathworks.com/help/comm/ug/multipath-fading-channel.html)를 바탕으로 작성한 한국어 실험 노트입니다. 기존 `matlab/AWGN` 코드는 변경하지 않습니다.

## 실행

필요 환경: **MATLAB R2025a 이상, Communications Toolbox, Signal Processing Toolbox**.

`study_multipath_fading.m`을 MATLAB **Live Editor**로 열고 **Run All**을 누르세요. 설명, 수식, 표, 그래프가 각 절에 표시되는 텍스트 형식 Live Script입니다. 처음에는 위에서 아래로 실행하고, 이후 설정을 바꿔 다시 전체 실행하세요.

프로젝트 루트가 MATLAB Current Folder이면 다음 명령으로 열 수 있습니다.

```matlab
open('matlab/MultipathFading/study_multipath_fading.m')
```

코드만 실행할 때는 다음 명령도 가능합니다. 일반 `run`에서는 중간 그래프가 다음 절에서 교체될 수 있으므로 학습에는 Live Editor를 권합니다.

```matlab
run('matlab/MultipathFading/study_multipath_fading.m')
```

## 구성

| 절 | 실험 |
|---|---|
| 1 | 경로 전력 정규화, PDP, 평균 지연, RMS delay spread |
| 2–3 | Rayleigh/Rician 이득, 진폭 분포, K-factor |
| 4–5 | Doppler에 따른 시간 변화, Jakes PSD, 자기상관 |
| 6–8 | RRC QPSK, 분수 지연, 실제 채널 통과, 광대역/협대역 응답과 성상도 |
| 9 | 선택 가능한 공식 채널 뷰어, `reset`/`release`/프레임 상태 |
| 10 | 결과 구조체와 파라미터 변경 과제 |

## 공식 예제와의 관계

4개 경로 지연 `[0 5 10 15] us`, 평균 전력 `[0 -3 -6 -9] dB`, 두 샘플률 500/20 kHz, 최대 Doppler 200 Hz, Rician K=10, 직접파 Doppler 100 Hz를 출발점으로 삼습니다.

학습을 위해 아래 내용을 추가하거나 변경했습니다.

- 신호 대역폭과 샘플률을 혼동하지 않도록 심볼당 4샘플, roll-off 0.25인 RRC 필터를 명시합니다. 두 경우의 심볼률은 125/5 ksym/s, 전체 이론 대역폭은 156.25/6.25 kHz입니다.
- 통계 실험은 시간에 따라 변하는 채널을, ISI 비교 실험은 **Doppler=0인 고정 랜덤 채널**을 사용합니다. 같은 seed로 두 샘플률의 물리적 경로 이득을 맞춥니다.
- 송수신 필터 지연과 채널 구현 지연을 보정한 뒤 성상도를 그립니다. 물리적 경로 간 지연은 없애지 않습니다.
- 입력 파형을 경로별 분수 지연 FIR과 이득으로 재구성하여 채널 출력과 일치하는지 검사합니다.
- 모든 송신 심볼을 안다는 가정의 단일 복소 이득 보상으로 ISI와 일정한 위상/진폭 변화를 구분합니다. 실용적인 수신기나 적응 등화기는 아닙니다.

Rician `KFactor`의 스칼라 설정은 **첫 경로에만** 직접파를 넣습니다. 경로별 전력 정규화는 평균적인 조건이며 매 프레임 수신 전력을 강제로 맞추지 않습니다. 자세한 정의는 [RicianChannel 문서](https://www.mathworks.com/help/comm/ref/comm.ricianchannel-system-object.html)를 참고하세요.

## 바꿔 볼 설정

스크립트 첫 절의 `mfConfig`를 수정합니다. 각 절은 앞 절의 변수를 사용합니다.

| 설정 | 영향 |
|---|---|
| `pathDelaysSeconds`, `averagePathGainsDb` | 지연 분산, PDP, 주파수 선택성. 두 벡터의 길이를 맞추세요. |
| `maximumDopplerHz` | 통계 실험의 시간 변화와 PSD. 기본 통계 Fs에서 0 초과 2000 Hz 이하. |
| `kFactorLinear`, `losDopplerHz` | Rician 통계 실험. K는 dB가 아닌 선형 전력비. |
| `sampleRatesHz`, `samplesPerSymbol`, `rolloff` | QPSK 심볼률과 대역폭. |
| `statisticsDurationSeconds` | 통계 관측 시간. 기본 5초. |
| `staticChannelSeed` | ISI 실험의 고정 채널 실현. 오차 수치는 seed에 따라 달라집니다. |
| `showBuiltInViewer` | `true`이면 공식 내장 impulse/frequency response 뷰어도 실행. |

`fadingStudy`에 설정, 표, 시간 이득, Doppler PSD/자기상관, 송수신 파형, 성상도 심볼 및 채널 응답이 모입니다. 스크립트는 파일을 자동 저장하지 않고 전역 RNG, MATLAB 검색 경로 또는 기존 workspace를 초기화하지 않습니다.

AWGN, 경로 손실, shadowing, BER, 5G NR 표준 TDL, FPGA 구현은 범위에 포함하지 않습니다. EVM은 한 채널 실현의 진단값으로, 규격 측정값이나 여러 채널에 대한 평균 성능이 아닙니다.
