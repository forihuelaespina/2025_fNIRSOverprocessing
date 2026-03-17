function Overprocessing0005_SegregationalAnalysis
%Overprocessing0005_SegregationalAnalysis
%
%
% Illustrates the overprocessing in the case of segregational analysis.
%
%% Remarks
%
% Uses ICNNA v1.4.1beta
%
%
%
%
% Copyright 2025-26
% @author Felipe Orihuela-Espina
%
% See also 
%


%% Log
%
% 8-Oct-2025: FOE
%   + File created.
%
%
% -- ICNNA v1.4.1beta
%
% 17-Mar-2026: FOE
%   + Updated for ICNNA v1.4.1beta.
%


opt.fontSize  = 18;
opt.lineWidth = 1.5;

%The following table summarises the simulated cases;
%
%
%  Original observation \ Analysis Outcome | Significant / Non-Signficant
%  ----------------------------------------+-----------------------------
%   Ch1) Random                            | No
%   Ch2) Random                            | Yes
%   Ch3) Active (no noise)                 | No
%   Ch4) Active (no noise)                 | Yes
%   Ch5) Active (plus noise)               | No
%   Ch6) Active (plus noise)               | Yes
%   Ch7) Not active (no noise)             | No
%   Ch8) Not active (no noise)             | Yes
%   Ch9) Not active (plus noise)           | No
%   Ch10) Not active (plus noise)          | Yes
%
% Note that lifting is needed for Ch 7 and 8 (here simulated by setting the
% no active as no contrast yet at amplitude 1).
%
% Noise is uniform and fixed to 20%.
%

seed = 1;
rng(seed);

nChannels = 10;

fs = 10; %Sampling frequency in [Hz]
t  = 0:(1/fs):120; %in [s]
nSamples = length(t);

hrfAmplitude = 1;
hrf = hrfAmplitude * HRF_DoubleGamma(t);

%Simulate 2x15sec blocks
boxcar = zeros(nSamples,1);
blocks = [30 15; ...
    75 15]; %[Onset duration] in secs
for iBlock = 1:size(blocks,1)
    theOnset  = blocks(iBlock,1) * fs;
    theOffset = (blocks(iBlock,1) + blocks(iBlock,2)) * fs ;
    boxcar(theOnset:theOffset) = 1;
end
q   = conv(boxcar,hrf,"full");
q(nSamples+1:end) = [];
tmpAmplitude = max(q);


%% Hypothesis
P = [0 1 0 1 0 1 0 1 0 1]'; 
% figure
% imagesc(P')

%% Observations.
Q(:,1:2)  = 0.2*tmpAmplitude*rand(nSamples, 2) - (1/2); %Random
Q(:,3:4)  = [q q]; %Active (no noise)
Q(:,5:6)  = [q q] + 0.2*tmpAmplitude*(rand(nSamples, 2) - (1/2)); %Active (plus noise)
Q(:,7:8)  = ones(nSamples, 2); %Not active (no noise) - Lifted
Q(:,9:10) = zeros(nSamples, 2) ...
    + 0.2*(tmpAmplitude*(rand(nSamples, 2) - (1/2))); %Not active (plus noise)

%% Pipeline and verification

extendedP = [tmpAmplitude*rand(size(Q)); tmpAmplitude*P'];  %The scaling prevents some rounding errors
                                   %when retriving the solution due to the
                                   %magnitude of the signal.
extendedQ = [Q; tmpAmplitude*rand(size(P'))];


tol = max(size(extendedQ))*eps(norm(extendedQ)); %Matlab default tolerance in pinv
A = extendedP*pinv(extendedQ,tol); %Pipeline
tmpExtendedResult = (A*extendedQ);
tmpResult = tmpExtendedResult(end-size(P,2)+1:end,:)';
tmpResult = round(tmpResult./max(tmpResult));


% Verifying the solution
disp(['Verification (0 means correct): ' ...
    num2str(any(any((P-tmpResult)>tol)))]);

%% Render
offsetFactor = 1.1;  %For plotting.

tt   = (1:nSamples)'.*(1/fs);
cmap = jet(nChannels);
legendStr(1,nChannels) = {''};
for iCh=1:nChannels
    legendStr(1,iCh) = {['Ch. ' num2str(iCh)]};
end

hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.88]);
hAxis(1) = subplot(3,1,1);
hLegend = plot(tt',(Q+offsetFactor*[1:nChannels])',...
    'LineStyle','-', 'LineWidth', opt.lineWidth);
title('Observations','FontSize',opt.fontSize);
legend(hLegend,legendStr,'FontSize',opt.fontSize);
xlabel(hAxis(1),'Time [sec]','FontSize',opt.fontSize);
ylabel(hAxis(1),'[A.U.]','FontSize',opt.fontSize);
set(hAxis(1),'XLim',[0 tt(end)]);

hAxis(2) = subplot(3,1,2);
hLegend = imagesc(P');
colormap('gray');
title('Hypothesis','FontSize',opt.fontSize);
set(hAxis(2),'XLim',[0.5 nChannels+0.5]);
set(hAxis(2),'YTick',[]);
xlabel(hAxis(2),'Channel','FontSize',opt.fontSize);

hAxis(3) = subplot(3,1,3); %hold on,
hLegend = imagesc(tmpResult');
colormap('gray');
title('Processed data','FontSize',opt.fontSize);
set(hAxis(3),'XLim',[0.5 nChannels+0.5]);
set(hAxis(3),'YTick',[]);
xlabel(hAxis(3),'Channel','FontSize',opt.fontSize);


set(hAxis(1),'YLimitMethod','padded');
set(hAxis,'Box','on');
set(hAxis(1),'XGrid','on','YGrid','on');
set(hAxis,'FontSize',opt.fontSize);

% mySaveFig(hFig,['..' filesep 'media' filesep ...
%     'Overprocessing0005_SegregationalAnalysis']);
% close(gcf);


end
