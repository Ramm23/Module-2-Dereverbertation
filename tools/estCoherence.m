function C = estCoherence(XL, XR, alpha)
%estCoherence   Estimate the short-term interaural coherence
%
%USAGE
% C = estCoherence(XL, XR, alpha)
%
%INPUT ARGUMENTS
%       XL: left ear STFT representation [M x L]
%       XR: right ear STFT representation [M x L]
%    alpha: smoothing coefficient which controls the averaging of the auto-
%    and cross power spectral densities across adjacent time frames
%
%OUTPUT ARGUMENTS
%       C: short-term interaural coherence [M x L]

% Compute auto- and cross-power 
phiLL = XL .* conj(XL);
phiRR = XR .* conj(XR);
phiLR = XL .* conj(XR);

% Perform smoothing
a = [1, -alpha];
b = [1, -alpha];

phiLL_smooth = filter(b,a,phiLL,[],2);
phiRR_smooth = filter(b, a, phiRR, [], 2);
phiLR_smooth = filter(b, a, phiLR, [], 2);

% Compute the short-term interaural coherence
C = phiLR_smooth ./ sqrt(phiLL_smooth.*phiRR_smooth);


end