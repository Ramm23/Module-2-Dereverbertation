clear
close all
clc

% Install subfolders
addpath irs
addpath signals
addpath tools


%% USER PARAMETERS
% 
% 
% Sampling frequency
fsHz = 16E3;

% Source signal
fileName = 'speech@24kHz.wav';

% Impulse response
%roomName = 'Cortex_45deg.wav';
roomName = 'Room_A_45deg.wav';
% roomName = 'Room_D_45deg.wav';
% roomName = 'YourOwnHATSRecording.wav';

% Window length
% winSec = ?????

% Integration time constant for coherence estimation in seconds
% tauSec = ?????


%% CREATE BINAURAL SIGNAL
% 
% 
% Load source signal
s = readAudio(fileName,fsHz);

% Load impulse response
h = readIR(roomName,fsHz);

% Fast convolution using the overlap-save method


%% PERFORM COHERENCE-BASED DEREVERBERATION
% 
% 

%% ANALYSIS USING THE DRR
% Create the reference signal without reverberation 
% direct impulse response:
[h_direct, h_reverb] = splitBRIR(h,fsHz);

% getting binaural signals
sL = convolveFFT_OLS(s,h_direct(:,1),256, false);
sR = convolveFFT_OLS(s,h_direct(:,2),256, false);

%DRR
xL = s(:,1); %original dry signal
xR = s(:,2); %original dry signal
sL_hat = (:,1)%derevered signal
sR_hat = (:,2)%dereverbed signal
DRR_pre = (sum(sL.^2) + sum(sR.^2))/(sum(xL-sL).^2 - sum(xR-sR).^2);
DRR_post = (sum(sL.^2) + sum(sR.^2))/(sum(sL_hat-sL).^2 - sum(sR_hat-sR).^2);
delta_DRR = DRR_post - DRR_pre;

%%HOW TO PICK TAUS???


%%
s_direct = horzcat(sL_direct, sR_direct);
sound(s,fsHz)
