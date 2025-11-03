function y = convolveFFT_OLS(x,h,nfft)
%convolveFFT_OLS  Linear convolution via overlap-save (streaming FFT).
%   y = convolveFFT_OLS(x,h)
%   y = convolveFFT_OLS(x,h,nfft)  % optional FFT length (power-of-2 >= M)
%
% Returns y with length N+M-1 (same as conv(x,h)).

    % --- shape & sizes
    x = x(:); h = h(:);
    N = numel(x); M = numel(h);
    if N==0 || M==0, y = []; return; end
    if M==1, y = x*h; return; end

    % --- FFT size (compact but safe default)
    if nargin < 3 || isempty(nfft)
        nfft = 2^nextpow2(max(1024, 2*M)); % >= M, power-of-2
    end
    if nfft < M
        nfft = 2^nextpow2(M);
    end
    L = nfft - (M - 1);                  % valid output per block

    % --- precompute filter FFT
    H = fft(h, nfft);

    % --- allocate output
    ylen = N + M - 1;
    y    = zeros(ylen,1);

    % --- input with headroom for tail; prepend M-1 zeros for OLS
    %     ensure enough total samples so the loop can read full nfft blocks
    numBlocks = ceil(ylen / L);
    needIn    = (numBlocks-1)*L + nfft;          % total samples needed incl. M-1 lead
    xpadlen   = max(0, needIn - ((M-1) + N));
    xz        = [zeros(M-1,1); x; zeros(xpadlen,1)];

    % --- process
    inIdx  = 1;    % index into xz for current block start
    outIdx = 1;    % index into y for writing valid samples

    for b = 1:numBlocks
        xb = xz(inIdx : inIdx + nfft - 1);       % M-1 overlap + L new
        Yb = ifft( fft(xb) .* H );               % circular conv
        valid = Yb(M:end);                       % drop first M-1 (overlap)
        k = min(L, ylen - outIdx + 1);           % last block may be short
        y(outIdx:outIdx+k-1) = valid(1:k);
        inIdx  = inIdx + L;
        outIdx = outIdx + k;
    end

    % (optional) strip tiny imag from roundoff
    if ~isreal(y), y = real(y); end
end
