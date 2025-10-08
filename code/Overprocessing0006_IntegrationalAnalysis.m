function Overprocessing0006_IntegrationalAnalysis
%Overprocessing0006_IntegrationalAnalysis
%
%
% Illustrates the overprocessing in the case of integrational analysis.
%
% Copyright 2025
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


opt.fontSize  = 18;
opt.lineWidth = 1.5;

%The following table summarises the simulated cases;
%                                          
%  Original observation \ Analysis Outcome 
%  ----------------------------------------
%   Ch1) Random
%   Ch2) Random
%   Ch3) Active (no noise)
%   Ch4) Active (no noise)
%   Ch5) Active (plus noise)
%   Ch6) Active (plus noise)
%   Ch7) Not active (no noise)
%   Ch8) Not active (no noise)
%   Ch9) Not active (plus noise)
%   Ch10) Not active (plus noise)
%
% Note that lifting is needed for Ch 7 and 8 (here simulated by setting the
% no active as no contrast yet at amplitude 1).
%
% Noise is uniform and fixed to 20%.
%
%
% Two hypothesis are simulated;
%
%   1) Functional connectivity - Random symmetric matrix
%   2) Effective  connectivity - Random asymmetric matrix
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

for iCase = 1:2

    opt.case = iCase;

    %% Hypothesis
    switch (opt.case)
        case 1 %1) Functional connectivity - Random symmetric matrix
            upperTri = triu(randi([0,1], nChannels, nChannels), 1);
            P = eye(nChannels) + upperTri + upperTri';
        case 2 %2) Effective  connectivity - Random asymmetric matrix
            P = or(eye(nChannels), randi([0,1], nChannels, nChannels));
        otherwise
            error('Unexpected case.');
    end

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
    hAxis(1) = subplot(2,2,[1 2]);
    hLegend = plot(tt',(Q+offsetFactor*[1:nChannels])',...
        'LineStyle','-', 'LineWidth', opt.lineWidth);
    title('Observations','FontSize',opt.fontSize);
    legend(hLegend,legendStr,'FontSize',opt.fontSize);
    xlabel(hAxis(1),'Time [sec]','FontSize',opt.fontSize);
    ylabel(hAxis(1),'[A.U.]','FontSize',opt.fontSize);
    set(hAxis(1),'XLim',[0 tt(end)]);

    hAxis(2) = subplot(2,2,3);
    hLegend = imagesc(P');
    colormap('gray');
    title('Hypothesis','FontSize',opt.fontSize);
    set(hAxis(2),'XLim',[0.5 nChannels+0.5]);
    set(hAxis(2),'XTick',[1:10]);
    set(hAxis(2),'YLim',[0.5 nChannels+0.5]);
    set(hAxis(2),'YTick',[1:10]);
    xlabel(hAxis(2),'Channel','FontSize',opt.fontSize);
    ylabel(hAxis(2),'Channel','FontSize',opt.fontSize);
    set(hAxis(2),'XAxisLocation', 'top')
    axis equal
    axis tight

    hAxis(3) = subplot(2,2,4); %hold on,
    hLegend = imagesc(tmpResult');
    colormap('gray');
    title('Processed data','FontSize',opt.fontSize);
    set(hAxis(3),'XLim',[0.5 nChannels+0.5]);
    set(hAxis(3),'XTick',[1:10]);
    set(hAxis(3),'YLim',[0.5 nChannels+0.5]);
    set(hAxis(3),'YTick',[1:10]);
    xlabel(hAxis(3),'Channel','FontSize',opt.fontSize);
    ylabel(hAxis(3),'Channel','FontSize',opt.fontSize);
    set(hAxis(3),'XAxisLocation', 'top')
    axis equal
    axis tight


    set(hAxis(1),'YLimitMethod','padded');
    set(hAxis,'Box','on');
    set(hAxis(1),'XGrid','on','YGrid','on');
    set(hAxis,'FontSize',opt.fontSize);

    mySaveFig(hFig,['..' filesep 'media' filesep ...
        'Overprocessing0006_IntegrationalAnalysis_' ...
        num2str(opt.case,'%04d')]);
    close(gcf);

end
end
