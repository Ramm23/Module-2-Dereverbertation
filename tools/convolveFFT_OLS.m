function [y] = convolveFFT_OLS(x, h, N, bZeroPad)
% convolveFFT_OLS  Overlap-Save convolution using FFT
% 
%   y = convolveFFT_OLS(x, h, N, bZeroPad)
%
%   Inputs:
%       x         - Input signal (column or row)
%       h         - Impulse response (filter)
%       N         - Block length for FFT processing
%       bZeroPad  - 1 to return full convolution (Nx + M - 1)
%                   0 to return same length as input (Nx)
%
%   Output:
%       y         - Convolved output signal

    L = N;

    Nx = length(x);
    M = length(h);
    Length = Nx + M - 1;
    Blocklenght = L; 
    Blocknumber = ceil(Nx / Blocklenght);
    
    x = x'; %transpose x because we pass it in as a column
    x = [zeros(1,M-1) x]; %add zeros to beginning

    x = [x zeros(1,mod(length(x),(Blocknumber*L+M-1)))]; %add zeros to end

    y = zeros(1,length(x)-(M-1)); 

    hx = zeros(1,Blocklenght+M-1);
    hx(1:length(h)) = h;
    h = hx;

    H = fft(h);
    
    % segmenting
    for seg = 0:Blocknumber-1
        x_segmented = x(seg*L+1:(seg+1)*L+(M-1));
        X_segmented = fft(x_segmented);
        Y = X_segmented .* H;
        y_segmented = ifft(Y);
        y_segmented = y_segmented(M:end);
        y(seg*L+1:(seg+1)*L) = y_segmented;
    end

        % Trim or pad based on bZeroPad flag
    if bZeroPad
        % Full convolution length
        y = y(1 : Length);
    else
        % Same length as input (like conv(x,h,'same'))
        y = y(M : Length);
    end

    y = y'; %transpose y because we want to return it as a column
end