function [hrf] = HRF_DoubleGamma(t,options)
%
%  [hrf] = HRF_DoubleGamma(t,options)
%
%% References
%
% Minako Uga, Ippeita Dan, Toshifumi Sano, Haruka Dan, Eiju Watanabe
% Optimizing the general linear model for functional near-infrared
% spectroscopy: an adaptive hemodynamic response function approach.
% Neurophotonics 1(1), 015004 (Jul–Sep 2014)
%
%% Parameters
%
% t - Double[].
%   Vector of time samples in [s]
%
% options - Struct with the following options
%   .tau_p - Scalar. Optional. Default is 6 [s].
%       first peak delay in [s]. 
%
%   .tau_d - Scalar. Optional. Default is 10 [s].
%       Delay of undershoot in [s]. 
%
%   .A - Scalar. Optional. Default is 6.
%       Amplitude ratio between the peaks
%
%
%% Output
%
% hrf - Double[]. Same length as t.
%   The hemodynamic response function as a function of time t.
%
%
%
% Copyright 2025
% @author Felipe Orihuela-Espina
%
% See also 
%


%% Log
%
% 3-Sep-2025: FOE
%   + Function created.
%

tau_p = 6;
tau_d = 10;
A = 6; 
if(exist('options','var'))
    %%Options provided
    if(isfield(options,'tau_p'))
        tau_p = options.tau_p;
    end
    if(isfield(options,'tau_d'))
        tau_d = options.tau_d;
    end
    if(isfield(options,'A'))
        A = options.A;
    end
end

firstTerm  = (t.^tau_p .* exp(-t))/factorial(tau_p);
secondTerm = (t.^(tau_p+tau_d) .* exp(-t))/(A*factorial(tau_p+tau_d));
hrf = (firstTerm - secondTerm)';

end