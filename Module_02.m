clear
close all
clc

% Install subfolders
addpath irs
addpath signals
addpath tools


% load irs\meas_2025_10_8_13_13_13.mat
% h = h_norm;
% filename = 'HATS_test.wav';
% audiowrite(filename,h,48E3);

%% USER PARAMETERS

% Sampling frequency
fsHz = 16E3;

% Source signal
fileName = 'speech@24kHz.wav';

% Impulse response
% roomName = 'Cortex_45deg.wav';
% roomName = 'Room_A_45deg.wav';          
roomName = 'Room_D_45deg.wav';      % reverberant room
% roomName = 'HATS_test.wav';


% Window length
winSec = 8E-3;

% Integration time constant for coherence estimation in seconds
% tauSec = 0.1;

% tau -> 0 -> exp(-) -> alpha = 0 -> less weight on prev frame (fast response)
% tau -> infty -> exp(0) -> alpha 1 -> more weight on prev frame (slow response)


%% CREATE BINAURAL SIGNAL


% Load source signal
s = readAudio(fileName,fsHz);


% Load impulse response
h = readIR(roomName,fsHz);


% sound(s, fsHz) % anechoic speech signal

hL = h(:,1);
hR = h(:,2);


% Fast convolution using the overlap-save method
xL = convolveFFT_OLS(s, hL);
xR = convolveFFT_OLS(s, hR);


% sound(xL, fsHz) % filtered signal

% plot original and filtered signals in time-domain
t = (0:length(s)-1) / fsHz;
tFilter = (0:length(xL)-1) /fsHz;
figure;
subplot(2,1,1)
plot(t, s)
xlabel('Time (sec)')
ylabel('Amplitude (a.u.)')
title('Original speech signal')
xlim([0 tFilter(end)])
grid on
subplot(2,1,2)
plot(tFilter, xR)
hold on
plot(tFilter, xL)
xlabel('Time (sec)')
ylabel('Amplitude (a.u.)')
title('Filtered speech signal')
xlim([0 tFilter(end)])
legend('Right channel', 'Left channel')
grid on


% Compute STFT of right and left ear room filtered signals
N = 2*round(winSec * fsHz / 2);     % Window size
R = N/4;                            % Step size
M = pow2(nextpow2(N));              % DFT size
w = cola(N,R,'hamming', 'wola');    % Window


[XL,tSTFT,fSTFT] = stft(xL,fsHz,w,R,M);
[XR,~,~] = stft(xR,fsHz,w,R,M);


dRdB = 60;


title1 = 'STFT of left channel';
title2 = 'STFT of right channel';

plotSTFT(tSTFT, fSTFT, XL, fsHz, false, dRdB, title1);
plotSTFT(tSTFT, fSTFT, XR, fsHz, false, dRdB, title2);


%% PERFORM COHERENCE-BASED DEREVERBERATION
tauSec = 0.1;


[sL_hat, sR_hat, G] = dereverb(xL, xR, fsHz, winSec, tauSec);




tDereverb = (0:length(sL_hat)-1) / fsHz;


% Plot dereverberated signals
figure;
subplot(3,1,1)
plot(t, s)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Anechoic Signal')
grid on
subplot(3,1,2)
plot(tFilter, xR)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Reverbed Left Channel Signal')
grid on
subplot(3,1,3)
plot(tDereverb, sR_hat)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Dereverberated Left Channel Signal')
grid on



%% Objective evaluation

% Binaural reference without reverb (convolve signal with BRIR from left
% and right channel)


% Create BRIR
brir = cat(2, hL, hR);

% Identify the direct-path component of the BRIR
hDirect = splitBRIR(brir, fsHz);

hDirect_L = hDirect(:,1);
hDirect_R = hDirect(:,2);

sL_ref = convolveFFT_OLS(s, hDirect_L);
sR_ref = convolveFFT_OLS(s, hDirect_R);





figure;
subplot(4,1,1)
plot(t, s)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Anechoic Signal')
grid on
subplot(4,1,2)
plot(tFilter, xR)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Reverbed XX Channel Signal')
grid on
subplot(4,1,3)
plot(tDereverb, sR_hat)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Dereverberated XX Channel Signal')
grid on
subplot(4,1,4)
plot(tDereverb, sR_ref)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Reference XX Channel Signal')
grid on



%% PLOT STFTs

% Compute STFT of right and left ear room filtered signals
N = 2*round(winSec * fsHz / 2);     % Window size
R = N/4;                            % Step size
M = pow2(nextpow2(N));              % DFT size
w = cola(N,R,'hamming', 'wola');    % Window


% sound(s, fsHz)
% sound(sL_ref, fsHz)
% sound(sL_hat, fsHz)
sound(xL, fsHz)


[SL,t_S,f_S] = stft(sL_ref,fsHz,w,R,M);      % STFT Dry Signal
[XL,t_XL,f_XL] = stft(xL,fsHz,w,R,M);       % STFT Reverberant Signal
[SL_hat,t_SL,f_SL] = stft(sL_hat,fsHz,w,R,M);       % STFT Dereverberant Signal




title1 = 'STFT of Dry Signal (left channel)';
title2 = 'STFT of Reverberant Signal (left channel)';
title3 = 'STFT of Dereverberant Signal (left channel)';

plotSTFT(tSTFT, fSTFT, SL, fsHz, false, dRdB, title1);
plotSTFT(tSTFT, fSTFT, XL, fsHz, false, dRdB, title2);
plotSTFT(tSTFT, fSTFT, SL_hat, fsHz, false, dRdB, title3);



%% Plot reverberant and enhanced signals together

tauSec = 0.1;

[sL_hat, sR_hat, G] = dereverb(xL, xR, fsHz, winSec, tauSec);

tDereverb = (0:length(sL_hat)-1) / fsHz;
tFilter = (0:length(xL)-1) /fsHz;


figure;
subplot(2,1,1)
plot(tFilter, xR)
hold on
plot(tDereverb, sR_hat)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Reverbed vs Dereverbed Signal (right channel)')
legend('Reverbered', 'Dereverbed')
grid on
subplot(2,1,2)
plot(tFilter, xL)
hold on
plot(tDereverb, sL_hat)
xlabel('Time (sec)')
xlim([0 2])
ylabel('Amplitude (a.u.)')
title('Reverbed vs Dereverbed Signal (left channel)')
legend('Reverbered', 'Dereverbed')
grid on


%% DRR improvement

xL = convolveFFT_OLS(s, hL);
xR = convolveFFT_OLS(s, hR);
[sL_hat, sR_hat, G] = dereverb(xL, xR, fsHz, winSec, tauSec);

sL = convolveFFT_OLS(s, hDirect_L);
sR = convolveFFT_OLS(s, hDirect_R);


DRR_pre = 10*log10((sum(sL.^2) + sum(sR.^2)) / (sum(xL-sL).^2 - sum(xR-sR).^2));
DRR_post = 10*log10((sum(sL.^2) + sum(sR.^2)) / (sum(sL_hat-sL).^2 - sum(sR_hat-sR).^2));
delta_DRR = DRR_post - DRR_pre;

disp(delta_DRR)



%% STFT of gain, G


tauSec = 0.05;

[sL_hat, sR_hat, G] = dereverb(xL, xR, fsHz, winSec, tauSec);

title4 = 'STFT of Gain function ';

plotSTFT(tSTFT, fSTFT, G, fsHz, false, 60, title4);


%%
GdB = 20*log10(G + eps);

figure;
imagesc(tSTFT, fSTFT*1e-3, GdB);
axis xy;
colormap(jet);
colorbar;
caxis([-40 0]);
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Estimated Gain (dB)');
