function w = genWin(N,type,opt)
%genWin   Create generalized cosine windows
%
%USAGE
%   w = genWin(N)
%   w = genWin(N,type,opt)
%
%INPUT ARGUMENTS
%      N : window size in samples
%   type : string specifying window type (default, type = 'hamming')
%          'rect'      = rectangular window
%          'hann'      = hann window
%          'hamming'   = hamming window 
%          'blackmann' = blackmann window
%          'flattop'   = flattop window
%    opt : string specifying window property (default, opt = 'periodic')
%          'periodic'  = unique maximum at N/2+1
%          'symmetric' = symmetric start and end-points
%
%OUTPUT ARGUMENTS
%   w : window function [N x 1]
% 
%   genWin(...) plots the window function in a new figure.
%
%   See also cola, iscola, stft and istft.

%   Developed with Matlab 9.2.0.538062 (R2017a). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/09/07
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 1 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(opt);  opt  = 'periodic'; end
if nargin < 2 || isempty(type); type = 'hamming';  end
    

%% CREATE GENERALIZED COSINE WINDOW
% 
% 
% Check if window should be periodic or symmetric
if strcmpi(opt,'periodic')
    N = N + 1;
end

% Force windows to be symmetric
if rem(N,2)
    % Odd-length window
    M = (N + 1) / 2;
else
    % Even-length window
    M = N / 2;
end

% Select window type
switch lower(type)
    case 'rect'
        a = 1; 
        b = 0; 
        c = 0;
        d = 0;
        e = 0;
    case 'hann'
        a = 0.5; 
        b = 0.5; 
        c = 0;
        d = 0;
        e = 0;
    case 'hamming'
        a = 0.54; 
        b = 0.46; 
        c = 0;
        d = 0;
        e = 0;
    case 'blackmann'
        a = 0.42; 
        b = 0.5; 
        c = 0.08;
        d = 0;
        e = 0;
    case 'flattop'
        a = 0.21557895; 
        b = 0.41663158; 
        c = 0.277263158;
        d = 0.083578947;
        e = 0.006947368;
    otherwise
        error('Window type "%s" is not supported.',lower(type));
end

% Normalized sample indices
n = (0:M-1)'/(N-1);

% Generalized form of cosine windows
w = a - b * cos(2*pi*n) + c * cos(4*pi*n) - d * cos(6*pi*n) ...
    + e * cos(8*pi*n);

% Force first element to zero
if strcmpi(type,'blackmann')
   w(1) = 0;
end

% Force window to be symmetric
if rem(N,2)
    % Odd-length window
    w = [w; w(end-1:-1:1)];
else
    % Even-length window
    w = [w; w(end:-1:1)];
end

% Check if window should be periodic or symmetric
if strcmpi(opt,'periodic')
    % Window was created for N + 1 samples, remove last element
    w(end) = []; N = N - 1;
end


%% PLOT WINDOW FUNCTION
%
%
% If no output is specified
if nargout == 0
    figure;
    hold on;
    plot(0:N-1,w,'color',[0 0.3895 0.9712])
    legend(lower(type))
    xlabel('$n$','interpreter','latex')
    ylabel('$w[n]$','interpreter','latex')
    grid on;
    xlim([0 N-1])
end
