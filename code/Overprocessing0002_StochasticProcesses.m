function Overprocessing_StochasticProcesses
%Overprocessing_StochasticProcesses
%
%
% Dealing with stochastic processes.
%
% Sampling from stochastic process illustrated using the inverse
% transform method (see auxiliary function sampling).
%
%
% Copyright 2025
% @author Felipe Orihuela-Espina
%
% See also 
%


%% Log
%
% 5-Sep-2025: FOE
%   + File created.
%


opt.fontSize  = 18;
opt.lineWidth = 1.5;

%Sampling spaces
% -- For the time series
fs  = 10; %Sampling frequency [Hz]
t   = (0:(1/fs):30)';
% -- For the stochastic process
sfs = 100; %Sampling frequency [AU] 



nSamples = size(t,1);



%Observations
% CentralTendency / Signal
Q_mu = sin(2*pi*1*t+0);
% Uncertainty / Noise
qe       = (-10:(1/sfs):10)';
qe_mu    = 0;
qe_sigma = 1;
Q_e  = (1/(qe_sigma*sqrt(2*pi))) * exp((-(qe-qe_mu).^2)/(2*qe_sigma^2)); % Gaussian process
%figure, plot(qe,Q_e);

%Illustrate the sampling
% [outcome] = sampling(size(t),[qe Q_e]);
% histogram(outcome);
Q_s = sampling([nSamples,1],[qe Q_e]);




%Hypothesis
% CentralTendency / Signal
P_mu = sin(2*pi*1*t+0);
% Uncertainty / Noise
pe   = (0:(1/sfs):10)';
pe_lambda = 1;
%P_e  = (e_lambda.^(pe)./factorial(floor(pe))) * exp(-e_lambda); % Poisson process
P_e  = (pe_lambda.^(pe)./gamma(pe)) * exp(-pe_lambda); % (Smoothed) Poisson process (using the gamma approximation)
%figure, plot(pe,P_e);

P_s = sampling([nSamples,1],[pe P_e]);

Q = Q_mu + Q_s;
P = P_mu + P_s;


%Central tendency
tol  = max(size(Q_mu))*eps(norm(Q_mu)); %Matlab default tolerance in pinv
A_mu = P_mu*pinv(Q_mu,tol); %Pipeline for central tendency
%Uncertainty
tol = max(size(Q_s))*eps(norm(Q_s)); %Matlab default tolerance in pinv
A_s = P_s*pinv(Q_s,tol); %Pipeline for uncertainty
%Verifying the solution
disp(['Verification on the central tendency (0 means correct): ' ...
        num2str(any(any((P_mu-A_mu*Q_mu)>tol)))]);
disp(['Verification on the uncertainty (0 means correct): ' ...
        num2str(any(any((P_s-A_s*Q_s)>tol)))]);



%% Render
tt   = (1:nSamples)'.*(1/fs);
nBins = 40;

hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
hAxis(1) = subplot(3,3,[1 2]);
plot(tt,Q,'Color','r',...
           'LineStyle','-', 'LineWidth', opt.lineWidth);
title('Observations','FontSize',opt.fontSize);
hAxis(2) = subplot(3,3,3);
histogram(Q_s,nBins);
title('Gaussian uncertainty','FontSize',opt.fontSize);

hAxis(3) = subplot(3,3,[4 5]);
plot(tt,P,'Color','b',...
           'LineStyle','-', 'LineWidth', opt.lineWidth);
title('Hypothesis','FontSize',opt.fontSize);
hAxis(4) = subplot(3,3,6);
histogram(P_s,nBins);
title('Poisson uncertainty','FontSize',opt.fontSize);


hAxis(5) = subplot(3,3,[7:9]); hold on,
hLegend(:,1)=plot(tt,P,'Color','b',...
           'LineStyle','-', 'LineWidth', opt.lineWidth);
hLegend(:,2)=plot(tt,(A_mu*Q_mu + A_s*Q_s ),'Color','g',...
           'LineStyle','--', 'LineWidth', opt.lineWidth);
title('Processed data','FontSize',opt.fontSize);
legend(hLegend(1,:),{'Hypothesis','Processed data'},'FontSize',opt.fontSize);


    set(hAxis([1,3,5]),'XLim',[0 tt(end)]);
    set(hAxis,'YLimitMethod','padded');
    set(hAxis,'Box','on');
    set(hAxis,'XGrid','on','YGrid','on');
    set(hAxis,'FontSize',opt.fontSize);
    xlabel(hAxis(3),'Time [sec]','FontSize',opt.fontSize);
    ylabel(hAxis,'[A.U.]','FontSize',opt.fontSize);


mySaveFig(hFig,['..' filesep 'media' filesep ...
        'Overprocessing0002_StochasticProcesses']);
close(gcf);

end

%% AUXILIARY FUNCTIONS
function [outcome] = sampling(s,pdf)
%Simulates sampling from a probability distribution function (pdf)
%
% [outcome] = sampling(s,pdf)
%
% This functions approximates the inverse transform method.
% 
%   1) Estimate the inverse cumulative distribution
%   function (iCDF)
%   2) generate a random number between 0 and 1 (uniform), 
%   3) plug that number into the iCDF to get the sampled value.
%
%% Input parameters
%
% s - Int[]. The desired size of the outcome.
% pdf - double[nx2]. A probability distribution function
%   The first column provides the sampling location and the second
%   column provides the probability value. For instance, for a Gaussian
%   distribution;
%
%       sfs = 100; %Sampling frequency [AU]
%       e   = (-10:(1/sfs):10)';
%       e_mu    = 0;
%       e_sigma = 1;
%       P_e  = (1/(e_sigma*sqrt(2*pi))) * exp((-(e-e_mu).^2)/(2*e_sigma^2));
%       pdf = [e P_e];
%
%
%% Output
%
% outcome - Double[Sized s]. A collection of random observations following pdf
%
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
% 6-Sep-2025: FOE
%   + File created.
%


seed = 1;
rng(seed);


% 1) Estimate the inverse cumulative distribution function (iCDF)
if sum(pdf(:,2)) ~= 1
    %Normalize
    pdf(:,2) = pdf(:,2)/sum(pdf(:,2));
end


% Regularize. Reduce flat regions and improve numerical stability.
smoothing_factor = 5;
pdf(:,2) = smooth(pdf(:,2), smoothing_factor);  


%Compute the CDF
cdf(:,2) = cumsum(pdf(:,2));
cdf(:,1) = pdf(:,1);



% Remove duplicate CDF values - Not conceptually needed, but numerically
% needed to avoid the "Sample points must be unique." error during
% interpolation.
[unique_cdf, idx] = unique(cdf(:,2));
unique_x = cdf(idx,1);
cdf_forInterp = [unique_x unique_cdf];

% figure, hold on
% plot(cdf(:,1),cdf(:,2),'r-');
% plot(cdf_forInterp(:,1),cdf_forInterp(:,2),'b--');




% Define a grid of cumulative probabilities
N = 99; %Number of quantile points
u = linspace(0, 1, N);  
% Inverse CDF via interpolation
x_inv = interp1(cdf_forInterp(:,2), cdf_forInterp(:,1), u, 'linear', 'extrap');
    %x_inv(i) such that P(X ≤ x_inv(i)) ≈ u(i)

iCDF = [u; x_inv]';
% figure,
% plot(iCDF(:,1),iCDF(:,2));

% 2) generate a random number between 0 and 1 (uniform), 
tmp = rand(s);

% 3) plug that number into the iCDF to get the sampled value.
outcome = interp1(iCDF(:,1), iCDF(:,2), tmp, 'linear', 'extrap');



end
