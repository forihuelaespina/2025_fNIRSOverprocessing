function Overprocessing0003_ExperimentalData
%Overprocessing0003_ExperimentalData
%
%
% Reach a desired hypothesis from experimental data.
%
%% Data
%
% Public data available from: https://openneuro.org/datasets/ds005963/versions/1.0.0
%
%
% Rickson C. Mesquita (2025). FRESH Motor Dataset. OpenNeuro. [Dataset]
% doi: doi:10.18112/openneuro.ds005963.v1.0.0
%
% Appelhoff, S., Sanderson, M., Brooks, T., Vliet, M., Quentin, R.,
% Holdgraf, C., Chaumon, M., Mikulan, E., Tavabi, K., Höchenberger,
% R., Welke, D., Brunner, C., Rockhill, A., Larson, E., Gramfort,
% A. and Jas, M. (2019). MNE-BIDS: Organizing electrophysiological
% data into the BIDS format and facilitating their analysis.
% Journal of Open Source Software 4: (1896).
% https://doi.org/10.21105/joss.01896
%
%
%
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
% 3-Sep-2025: FOE
%   + File created.
%
% 3-Oct-2025: FOE
%   + Re-run to simplify the output figure
%
% -- ICNNA v1.4.1beta
%
% 17-Mar-2026: FOE
%   + Updated for ICNNA v1.4.1beta.
%


opt.fontSize  = 18;
opt.lineWidth = 1.5;

seed = 1;
rng(seed);


%% Preliminaries
iSubject  = 1;
iSession  = 1;
switch (iSession)
    case 1
        sessName = 'left2s';
    case 2
        sessName = 'left3s';
    case 3
        sessName = 'right2s';
    case 4
        sessName = 'right3s';
end

% Load the data file
srcFolder = ['..' filesep 'data' filesep];
filename  = ['sub-' num2str(iSubject,'%02d') ...
    '_ses-' sessName '_task-FRESHMOTOR_nirs.snirf'];


r = rawData_Snirf();
r = r.import([srcFolder 'sub-' num2str(iSubject,'%02d') filesep ...
    'ses-' sessName filesep ...
    'nirs' filesep filename]); %Load the .snirf file.
nimg = r.convert(); %Convert to Hb.

%The following info is available from companion .json file;
%
% sub-01_ses-left2s_task-FRESHMOTOR_nirs.json
%
%...but it is manually added here for simplicity.
%nMeasurements = 136;
SamplingFrequency = 8.928571428571429; %in [Hz]



%Simulated cases:
% 1 - Hypothesis; Random
% 2 - Hypothesis; All channels active
% 3 - Hypothesis; All channels inactive
for iCase = 1:3

    opt.case = iCase;


    nChannels = nimg.nChannels;
    nSamples  = nimg.nSamples;
    offsetFactor = 1.4;


    %% Generate the hypothesis
    % P (hypothesis) and Q (observations)

    %Build the active hypothesis for a single channel
    tmpCond = nimg.timeline.getConditions('motor');
    cevents = tmpCond.cevents;
    nEvents = size(cevents,1);
    boxcar  = zeros(nSamples,1);
    for iEv = 1:nEvents
        boxcar(cevents(iEv).onsets:cevents(iEv).onsets+cevents(iEv).durations) = ...
            cevents(iEv).amplitudes;
    end

    t   = (0:1/SamplingFrequency:30);
    hrf = HRF_DoubleGamma(t);
    p   = conv(boxcar,hrf,"full");
    p(nSamples+1:end) = [];

    % figure, hold on
    % plot(boxcar,'r-')
    % plot(p,'b-')

    %Generate a hypothesis for every channel
    switch (opt.case)
        case 1 %Land in random
            activeChannels = rand(1,nChannels)>0.5;
            P = zeros(nSamples,nChannels);
            P(:,activeChannels) = repmat(p,1,sum(activeChannels));

        case 2 %Land in all active
            P = repmat(p,1,nChannels);

        case 3 %Land in all inactive
            P = zeros(nSamples,nChannels);

        otherwise
            error('Unexpected case.');
    end

    % figure
    % plot(P+offsetFactor*[1:nChannels]);


    %% Extract the experimental observations
    Q = nimg.data(:,:,1); %Use HbO2 only for exemplary purposes.


    %% Set the pipeline
    tol = max(size(Q))*eps(norm(Q)); %Matlab default tolerance in pinv
    A = P*pinv(Q,tol); %Pipeline
    %Verifying the solution
    disp(['Verification (0 means correct): ' num2str(any(any((P-A*Q)>tol)))]);



    %% Render
    tt   = (1:nSamples)'.*(1/SamplingFrequency);

    %Pick only a subset of channels for visualization purposes.
    tmpVisibleChannelsMask = rand(1,nChannels)<0.25;
    nVisibleChannels       = sum(tmpVisibleChannelsMask);

    %cmap = jet(nChannels);
    %legendStr(1,nChannels) = {''};
    % for iCh=1:nChannels
    %     legendStr(1,iCh) = {['Ch. ' num2str(iCh)]};
    % end
  

    hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
    hAxis(1) = subplot(3,1,1);
    hold on
    plot(tt,Q(:,tmpVisibleChannelsMask)',...
        'LineStyle','-', 'LineWidth', opt.lineWidth);
    title('Observations','FontSize',opt.fontSize);

    hAxis(2) = subplot(3,1,2);
    plot(tt,(P(:,tmpVisibleChannelsMask)+offsetFactor*(1:nVisibleChannels))',...
        'LineStyle','-', 'LineWidth', opt.lineWidth);
    title('Hypothesis','FontSize',opt.fontSize);

    hAxis(3) = subplot(3,1,3); hold on,
    hLegend(:,1)=plot(tt,(P(:,tmpVisibleChannelsMask)+offsetFactor*(1:nVisibleChannels))',...
        'LineStyle','-', 'LineWidth', opt.lineWidth);
    tmpProc = A*Q;
    hLegend(:,2)=plot(tt,(tmpProc(:,tmpVisibleChannelsMask))+offsetFactor*(1:nVisibleChannels),...
        'Color','g',...
        'LineStyle','--', 'LineWidth', opt.lineWidth);
    title('Processed data','FontSize',opt.fontSize);


    set(hAxis,'XLim',[0 tt(end)]);
    set(hAxis,'YLimitMethod','padded');
    set(hAxis,'Box','on');
    set(hAxis,'XGrid','on','YGrid','on');
    set(hAxis,'FontSize',opt.fontSize);
    xlabel(hAxis(3),'Time [sec]','FontSize',opt.fontSize);
    ylabel(hAxis,'[A.U.]','FontSize',opt.fontSize);


    for iAx = 1:3
        axes(hAxis(iAx));
        tmpYLim = ylim();
        line(([[cevents.onsets]' [cevents.onsets]'].*(1/SamplingFrequency))',...
            repmat([tmpYLim(1) tmpYLim(2)],nEvents,1)',...
            'Color','k',...
            'LineStyle','-','LineWidth',opt.lineWidth);
    end
    axes(hAxis(3))
    legend(hLegend(1,:),{'Hypothesis','Processed data'},'FontSize',opt.fontSize);


    % mySaveFig(hFig,['..' filesep 'media' filesep ...
    %     'Overprocessing0003_ExperimentalData' ...
    %     num2str(opt.case,'%04d')]);
    % close(gcf);

    clear hLegend

end



end