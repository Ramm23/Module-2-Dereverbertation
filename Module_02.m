clear; close all; clc;

% Install subfolders
addpath irs
addpath signals
addpath tools

%% USER PARAMETERS
fsHz    = 16e3;                 % Sampling frequency
fileName = 'speech@24kHz.wav';  % Source signal (will be resampled in readAudio)
roomName = 'Room_D_45deg.wav';  % Impulse response (BRIR) : 'Room_A_45deg.wav' 'Room_D_45deg.wav' 'Cortex_45deg.wav' 'HATS.wav'

% Window length (analysis window) and coherence smoothing time constant
winSec = 8e-3;      % per module default (8 ms)  
tauSec = 0.06;     % time constant for coherence smoothing (tune in experiments)

%% CREATE BINAURAL SIGNAL
% Load source signal (mono)
s = readAudio(fileName, fsHz);

% Load binaural impulse response (BRIR): left/right channels
% (make sure readIR returns [hL, hR] at fsHz)
[hFile, fsIR] = audioread(roomName);
if size(hFile,2) ~= 2, error('BRIR must be stereo'); end
if fsIR ~= fsHz, hFile = resample(hFile, fsHz, fsIR); end

% Optional: trim leading silence to the first peak (keeps headroom)
[~,iL] = max(abs(hFile(:,1))); [~,iR] = max(abs(hFile(:,2)));
lead = max(0, min(iL,iR)-1);
hFile = hFile(lead+1:end, :);

% **Normalize BRIR** (peak normalization is simple & safe here)
hFile = hFile ./ max(abs(hFile), [], 'all');

hL = hFile(:,1);
hR = hFile(:,2);


% Pad the shorter ear so both BRIRs have identical length
Lh = max(numel(hL), numel(hR));
if numel(hL) < Lh, hL(end+1:Lh) = 0; end
if numel(hR) < Lh, hR(end+1:Lh) = 0; end

xL = convolveFFT_OLS(s, hL);
xR = convolveFFT_OLS(s, hR);


%% STFT PARAMS 
% - N even
% - R = 25% of N (WOLA step)
% - M = next power of two >= N
N  = 2*floor((winSec*fsHz)/2);       % even length
R  = round(0.25 * N);                % 25% hop size  
M  = 2^nextpow2(N);                  % DFT size      

% Create analysis/synthesis window enforcing COLA (Hamming, WOLA)
% (uses provided cola.m)
w = cola(N, R, 'hamming', 'wola');  

%% STFT ANALYSIS (left/right)  [Eq. (3.1)]
XL = stft(xL, fsHz, w, R, M);
XR = stft(xR, fsHz, w, R, M);

%% SHORT-TERM INTERAURAL COHERENCE (IC)
% alpha from tau, hop, and Fs  [Eq. (3.12)]
alpha = exp(-R/(tauSec*fsHz));

% Estimate complex short-term coherence 
C = estCoherence(XL, XR, alpha);    

% Coherence-based gain 
G = abs(C).^2;

%% APPLY GAIN (magnitude domain)  
SL = (abs(XL) .* G) .* exp(1j*angle(XL));
SR = (abs(XR) .* G) .* exp(1j*angle(XR));

%% ISTFT SYNTHESIS (WOLA)  
sL = istft(SL, w, R, 'wola', numel(xL)); %dereverbed signal left
sR = istft(SR, w, R, 'wola', numel(xR)); %dereverbed signal right

%% LISTEN / SAVE / QUICK LOOKS

% --- Combine into stereo signals (L/R)
revLR = [xL(:), xR(:)];
revLR = revLR ./ max(abs(revLR), [], 'all') * 0.9;   % Reverberant stereo keep headroom (-1 dBFS-ish)   


derevLR = [sL(:), sR(:)];
derevLR = derevLR ./ max(abs(derevLR), [], 'all') * 0.9;    % Dereverberated stereo

% --- Clean up numerical issues (avoid NaNs/Infs)
revLR(~isfinite(revLR))   = 0;
derevLR(~isfinite(derevLR)) = 0;

% --- Normalize both files to -1 dBFS independently
targetAmp = 10^(-1/20);  % -1 dBFS ≈ 0.891
peakRev   = max(abs(revLR), [], 'all');
peakDerev = max(abs(derevLR), [], 'all');

if peakRev > 0
    revLR = revLR * (targetAmp / peakRev);
end
if peakDerev > 0
    derevLR = derevLR * (targetAmp / peakDerev);
end



% --- Save as 24-bit stereo WAVs
% Create output folder name
outDir = 'simulated_audio';

% Check if folder exists; if not, create it
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

[~, signal, ~] = fileparts(fileName);
[~, room, ~]   = fileparts(roomName);

% Build full paths for the output files
fn_reverb   = fullfile(outDir, strcat(signal, '_', room, '_reverb_LR.wav'));
fn_dereverb = fullfile(outDir, strcat(signal, '_', room, '_dereverb_LR.wav'));



% Save audio files
audiowrite(fn_reverb,   revLR,   fsHz, 'BitsPerSample', 24);
audiowrite(fn_dereverb, derevLR, fsHz, 'BitsPerSample', 24);

disp(['Saved files to folder: ' outDir]);


%% --- Plotting left ear as example

% Align lengths (just before plotting)
Lmin = min([numel(xL), numel(sL)]);
xLplt = xL(1:Lmin);
sLplt = sL(1:Lmin);

t = (0:Lmin-1)/fsHz;

% TIME DOMAIN OVERLAY
figure('Name','Time-domain (Left Ear)');
plot(t, xLplt, 'LineWidth', 0.8); hold on;
plot(t, sLplt, 'LineWidth', 0.8);
grid on; xlabel('Time [s]'); ylabel('Amplitude');
legend('Reverberant L','Dereverberated L','Location','best');
title('Reverberant vs Dereverberated (Left)');


% SPECTROGRAM
figure;
subplot(2,1,1);
spectrogram(xLplt, hamming(512), 384, 512, fsHz, 'yaxis');
title('Reverberant Left');
subplot(2,1,2);
spectrogram(sLplt, hamming(512), 384, 512, fsHz, 'yaxis');
title('Dereverberated Left');

%% ANALYSIS USING THE DRR
% Create the reference signal without reverberation 
% --- Build the DRR reference from the direct-path BRIR ---
[h_direct, ~] = splitBRIR(hFile, fsHz);           % use the same (trimmed, normalized) hFile

sL_dir = convolveFFT_OLS(s, h_direct(:,1));       % clean direct-path left
sR_dir = convolveFFT_OLS(s, h_direct(:,2));       % clean direct-path right

% Raw reverberant (from full BRIR) and raw dereverbed (from ISTFT)
xL_rev = xL;  xR_rev = xR;                        % reverberant
sL_der = sL;  sR_der = sR;                        % dereverbed

% Align lengths
Lmin = min([numel(sL_dir), numel(xL_rev), numel(sL_der)]);
sL_dir = sL_dir(1:Lmin); sR_dir = sR_dir(1:Lmin);
xL_rev = xL_rev(1:Lmin); xR_rev = xR_rev(1:Lmin);
sL_der = sL_der(1:Lmin); sR_der = sR_der(1:Lmin);

% DRR per assignment
num      = sum(sL_dir.^2) + sum(sR_dir.^2);
den_pre  = sum((xL_rev - sL_dir).^2) + sum((xR_rev - sR_dir).^2);
den_post = sum((sL_der - sL_dir).^2) + sum((sR_der - sR_dir).^2);

DRR_pre   = 10*log10(num / max(den_pre,  1e-20));
DRR_post  = 10*log10(num / max(den_post, 1e-20));
delta_DRR = DRR_post - DRR_pre


%GAIN PLOT
Fpos = (0:M/2)*(fsHz/M);                 % one-sided freqs, length M/2+1
Tstf = (0:size(G,2)-1)*(R/fsHz);         % your STFT frame times
Gpos_dB = 20*log10(max(G(1:numel(Fpos),:), 1e-4));

figure('Name','Gain (|C|^2)'); 
imagesc(Tstf, Fpos, Gpos_dB); axis xy; colorbar;
xlabel('Time [s]'); ylabel('Frequency [Hz]');
clim([-40 0]); title('Coherence-based Gain (dB)');
