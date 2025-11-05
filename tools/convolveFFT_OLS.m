function y = convolveFFT_OLS(x,h,N,bZeroPad)
% convolveFFT_OLS  Linear convolution using Overlap-Save FFT method
%
% y = convolveFFT_OLS(x,h,N,bZeroPad)
%
% x        : input signal [Nx x 1]
% h        : impulse response [M x 1]
% N        : DFT size (power of 2). If empty, uses optimalN().
% bZeroPad : if true, zero-pad x with M-1 samples to get full conv length
%
% y        : linear convolution result (Nx or Nx+M-1 samples)

x = x(:);
h = h(:);
Nx = length(x);
M  = length(h);

% --- Choose DFT size ---
if nargin < 3 || isempty(N)
    [N, ~] = optimalN(Nx,M);
end
if N < M
    error('N must be >= M for overlap-save to work.');
end

% --- Parameters ---
L = N - (M-1);        % number of "new" samples per block
H = fft(h,N);         % N-point FFT of impulse response

% --- Optional zero-padding at end ---
if nargin < 4 || bZeroPad
    x = [x; zeros(M-1,1)];
    Nx = length(x);
end

% --- Pad first block with M-1 zeros (overlap-save) ---
x = [zeros(M-1,1); x];

% --- Preallocate result ---
y = zeros(Nx,1);

% --- Process blocks ---
numBlocks = ceil((length(x)- (M-1))/L);
outIdx = 1;

for b = 0:numBlocks-1
    idx = b*L + (1:N);        % take N samples with overlap
    if idx(end) > length(x)   % zero-pad last block if needed
        x_block = zeros(N,1);
        x_block(1:length(x)-b*L) = x(idx(idx<=length(x)));
    else
        x_block = x(idx);
    end

    % FFT-based circular convolution
    Y_block = ifft(fft(x_block).*H);

    % Discard first M-1 corrupted samples
    valid = real(Y_block(M:end));

    % Write L valid samples to output
    endIdx = min(outIdx+L-1,length(y));
    y(outIdx:endIdx) = valid(1:(endIdx-outIdx+1));
    outIdx = outIdx + L;
end
end
