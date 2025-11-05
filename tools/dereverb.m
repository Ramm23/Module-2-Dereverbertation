function [sL, sR, G] = dereverb(xL, xR, fs, winSec, tauSec)
%dereverb Coherence−based dereverberation algorithm
%
%USAGE
% [sL,sR,G] = dereverb(xL,xR,fs)
% [sL,sR,G] = dereverb(xL,xR,fs,winSec,tauSec)
%
%INPUT ARGUMENTS
% xL : left ear input signal [nSamples x 1]
% xR : right ear input signal [nSamples x 1]
% fs : sampling frequency in Hertz
% winSec : window size in seconds across which the short−term
% interaural coherence function is estimated
% (default, winSec = 8E−3)
% tauSec : time constant controlling the smoothing of the auto− and
% cross−power spectral densitiy estimates across time
% (default, tauSec = ??????? )
%
%OUTPUT ARGUMENTS
% sL : enhanced left ear signal [nSamples x 1]
% sR : enhanced right ear signal [nSamples x 1]
% G : time−varying gain function


if nargin < 4 || isempty(winSec), winSec = 8E-3; end
if nargin < 5 || isempty(tauSec), tauSec = 0.1; end
if nargin < 3 || isempty(fs), fs = 16E3; end


% Compute STFT for each channel
N = 2*round(winSec * fs / 2);       % Window size
R = N/4;                            % Step size
M = pow2(nextpow2(N));              % DFT size
w = cola(N,R,'hamming', 'wola');    % Window

[XL,~,~] = stft(xL,fs,w,R,M);
[XR,~,~] = stft(xR,fs,w,R,M);


% Compute smoothing coefficient alpha
alpha = exp(-R / (tauSec * fs));


% Compute interaural coherence function (IC)
C = estCoherence(XL, XR, alpha);

% Compute gain function from IC
G = abs(C).^2;


% Apply the coherence gain to the left and right signals
SL_hat = abs(XL) .* G .* exp(1j * angle(XL));
SR_hat = abs(XR) .* G .* exp(1j * angle(XR));



% Compute ISTFT
sL = istft(SL_hat, w, R, 'wola', numel(xL));
sR = istft(SR_hat, w, R, 'wola', numel(xR));