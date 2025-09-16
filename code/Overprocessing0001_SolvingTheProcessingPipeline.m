function Overprocessing0001_SolvingTheProcessingPipeline
%Overprocessing0001_SolvingTheProcessingPipeline
%
%
% Reach a desired hypothesis despite random observations.
%
% Copyright 2025
% @author Felipe Orihuela-Espina
%
% See also 
%


%% Log
%
% 3-Sep-2025: FOE
%   + File created.
%
% 6-Sep-2025: FOE
%   + Added the case 5 - lifting example.
%   + Beautified the figures.
%
% 9-Sep-2025: FOE
%   + Added the case 6.
%


opt.fontSize  = 18;
opt.lineWidth = 1.5;

%Simulated cases:
% 1 - Hypothesis; All random (fixed seed)
% 2 - Hypothesis; all channels HRF
% 3 - Hypothesis; 1 HRF, 1 flat,  1 random (fixed seed)
% 4 - Hypothesis; "Rotated" bases i.e. P (hypothesis) and Q (observations)
%                   are two orthogonal signals each composed of a single
%                   frequency e.g. A=exp(jw_1t) and B=exp(jw_2t)
% 5 - Lifting example: Observations: Flat; Hypothesis; same as 3.
% 6 - Sending to infty: Observations: Random; Hypothesis; has some infinity.
for iCase = 6%1:6

    opt.case = iCase;

    nChannels = 3;

    fs = 10; %Sampling frequency in [Hz]
    t  = 0:(1/fs):30; %in [s]
    hrfAmplitude = 2.5;
    hrf = hrfAmplitude * HRF_DoubleGamma(t);


    %% Hypothesis
    switch (opt.case)
        case 1 %Random to random
            seed = 1;
            rng(seed);

            tmpSamples = 300;
            P = rand(tmpSamples, nChannels); %Hypothesis.

            offsetFactor = 1.1;  %For plotting. For random examples

        case 2 %Random to activity
            P  = repmat(hrf,2,nChannels); %Hypothesis.
            offsetFactor = 0.08; %For plotting. For double gamma HRF

        case 3 %Random to heterogeneous
            seed = 1;
            rng(seed);
            P  = repmat(hrf,2,1);
            tmpSamples = size(P,1);
            P = [P zeros(tmpSamples,1) rand(tmpSamples, 1)];
            %Hypothesis: HRF + flat + random.

            offsetFactor = 0.5; %For plotting. For double gamma HRF

        case 4
            nChannels = 1;
            w1 = 10/(2*pi);
            P  = imag(exp(1i * w1 * t))'; %Hypothesis.

            offsetFactor = 1; %For plotting. For double gamma HRF

        case 5  %Lifting and vectorization
            seed = 1;
            rng(seed);
            P  = repmat(hrf,2,1);
            tmpSamples = size(P,1);
            P = [P zeros(tmpSamples,1) rand(tmpSamples, 1)];
            %Hypothesis: HRF + flat + random.

            offsetFactor = 0.5; %For plotting. For double gamma HRF

        case 6 %Random to infinity
            seed = 1;
            rng(seed);

            p = tan(0:pi/32:4*pi)';
            p(p>10*16)    = Inf; %Force the infinity instead of MATLAB's rounding
            p(p<-1*10*16) = -Inf;
            tmpSamples = length(p);
            P = repmat(p,1,nChannels); %Hypothesis.

            offsetFactor = 1.1;  %For plotting. For random examples


        otherwise
            error('Unexpected case.');
    end

    nSamples = size(P,1);

    % figure
    % plot(P)

    %% Observations.
    switch (opt.case)
        case {1,2,3,6}
            Q = rand(nSamples, nChannels); 
        case 4
            w2 = 15/(2*pi);
            Q  = imag(exp(1i * w2 .* t))';
        case 5
            Q = zeros(nSamples, nChannels);
    end

    %% Pipeline and verification

    if opt.case == 5 %Lifting and vectorization
        %Lift and vectorized P and Q
        P_lifted = [P; zeros(1,nChannels)];
        P_lifted = reshape(P_lifted,numel(P_lifted),1);
        Q_lifted = [Q; ones(1,nChannels)];
        Q_lifted = reshape(Q_lifted,numel(Q_lifted),1);

        tol = max(size(Q_lifted))*eps(norm(Q_lifted)); %Matlab default tolerance in pinv
        A_lifted = P_lifted*pinv(Q_lifted,tol); %Pipeline
        
        tmpResult = reshape(A_lifted*Q_lifted,nSamples+1,nChannels);
        %Unlift result
        tmpResult(end,:) = [];

    else
        tol = max(size(Q))*eps(norm(Q)); %Matlab default tolerance in pinv
        A = P*pinv(Q,tol); %Pipeline
        tmpResult = (A*Q);
    end

    % Verifying the solution
    disp(['Verification (0 means correct): ' ...
            num2str(any(any((P-tmpResult)>tol)))]);

    %% Render
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

    hAxis(2) = subplot(3,1,2);
    hLegend = plot(tt',(P+offsetFactor*[1:nChannels])',...
            'LineStyle','-', 'LineWidth', opt.lineWidth);
    title('Hypothesis','FontSize',opt.fontSize);
    legend(hLegend,legendStr,'FontSize',opt.fontSize);

    hAxis(3) = subplot(3,1,3); hold on,
    hLegend(:,1)=plot(tt,P+offsetFactor*[1:nChannels],...
        'LineStyle','-', 'LineWidth', opt.lineWidth);
    hLegend(:,2)=plot(tt,tmpResult+offsetFactor*[1:nChannels],...
        'Color','g',...
        'LineStyle','--', 'LineWidth', opt.lineWidth);
    title('Processed data','FontSize',opt.fontSize);
    legend(hLegend(1,:),{'Hypothesis','Processed data'},...
            'FontSize',opt.fontSize);


    set(hAxis,'XLim',[0 tt(end)]);
    set(hAxis,'YLimitMethod','padded');
    set(hAxis,'Box','on');
    set(hAxis,'XGrid','on','YGrid','on');
    set(hAxis,'FontSize',opt.fontSize);
    xlabel(hAxis(3),'Time [sec]','FontSize',opt.fontSize);
    ylabel(hAxis,'[A.U.]','FontSize',opt.fontSize);

    % mySaveFig(hFig,['..' filesep 'media' filesep ...
    %     'Overprocessing0001_SolvingTheProcessingPipeline_Case' ...
    %     num2str(opt.case,'%04d')]);
    % close(gcf);


end

end
