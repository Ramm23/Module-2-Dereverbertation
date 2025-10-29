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
roomName = 'Cortex_45deg.wav';
% roomName = 'Room_A_45deg.wav';
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


