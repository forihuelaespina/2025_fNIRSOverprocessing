function Overprocessing0007_OperationsAsPointShiftsForfNIRS
%Overprocessing0007_OperationsAsPointShiftsForfNIRS
%
%
% Illustrates the point of expressing operations as point shifts in
%the signals space but in a manner more closely relatable to fNIRS
%than the generic sinusoidals.
%
% In this example, we simulate two fNIRS channels as linear combinations
% of canonical basis components (HRF, drift, cardiac, Mayer-wave) plus noise.
%
% Each processing step is expressed as a linear operator (matrix) acting
% on the observed signals P:
%
% 1) High-pass filter (HP):
%    - Removes slow baseline drift component from each channel.
%    - Constructed by zeroing the drift coefficient in the true weights:
%      w_HP = w_true; w_HP(Drift) = 0
%    - Operator: A_HP = Q_HP * pinv(P)
%
% 2) Physiological regression (PR):
%    - Removes fast physiological fluctuations (cardiac + Mayer-wave).
%    - Constructed by zeroing cardiac and LFO coefficients in the true weights:
%      w_PR = w_true; w_PR([Cardiac, LFO]) = 0
%    - Operator: A_PR = Q_PR * pinv(P)
%
% 3) Overprocessing (OV):
%    - Projects the observed signals to the target hypothesis (Q).
%    - Operator: A_over = Q * pinv(P)
%
% The resulting processed signals:
%    X_HP   = A_HP * P
%    X_PR   = A_PR * P
%    X_over = A_over * P
%
% These operators and signals can be visualized in:
%    - Time-domain plots
%    - 2D signal-space plots (showing shifts along basis components)
%    - 3D signal-space plots (visualizing trajectories in HRF–Drift–Cardiac space)
%
% This framework illustrates the conceptual point: each processing
% operation corresponds to a **shift along a meaningful component** in the
% signal space.
%
%
%% Remarks
%
% Uses ICNNA v1.4.1beta
%
%
%
%
%
% Copyright 2026
% @author Felipe Orihuela-Espina
%
% See also
%

%% Log
%
% 10-Mar-2026: FOE
%   + File created. In response to reviewers of JAIHC.
%
% -- ICNNA v1.4.1beta
%
% 17-Mar-2026: FOE
%   + Updated for ICNNA v1.4.1beta.
%

opt.fontSize  = 20;
opt.lineWidth = 1.5;
opt.markerSize = 14;


seed = 1;
rng(seed);


%% Reviewer comment context
%
% In the original submission, Figures 4-8 were based on a generic example
% on sinusoidals. The reviewer indicated:
%
% "I agree with the authors' intention to employ highly general
% examples to introduce their conceptual framework..."
%
% So I give it a thought and what follows is a didactic fNIRS-specific
% example to replace the sinusoid-only illustration from Figures 4–8.
%
% In principle, fNIRS signals are well-described as the sum of:
%
%   + an HRF (slow canonical response ~0.01–0.1 Hz)
%   + slow baseline drift (~0.005–0.02 Hz)
%   + cardiac pulsation (~1 Hz)
%   + Mayer waves (~0.1 Hz)
%   + additive noise
%
% These components behave like basis vectors in the signal space.
% Shifting along one component corresponds to a meaningful processing
% step (high-pass, regression, physiological noise removal).
%

%% 1) Define basis signals (HRF, drift, cardiac, LFO)

% Time base
Fs = 10;              % 10 Hz typical fNIRS sampling
T  = 40; %120;             % In seconds
t  = (0:1/Fs:T-1/Fs)';

% HRF
opt2.tau_p = 6;
opt2.tau_d = 10;
opt2.A = 6;
hrf = HRF_DoubleGamma(t,opt2);
hrf = hrf(:); % ensure column vector
hrf = hrf / max(hrf);

% Linear drift
drift = linspace(0,1,length(t))';
drift = drift(:); % ensure column vector
drift = drift - mean(drift);

% Cardiac component (~1 Hz)
card = sin(2*pi*1*t);
card = card(:); % ensure column vector
card = card / max(abs(card));

% Mayer-wave (~0.1 Hz)
lfo = sin(2*pi*0.1*t);
lfo = lfo(:); % ensure column vector
lfo = lfo / max(abs(lfo));

% Pack basis into matrix
B = [hrf, drift, card, lfo];

%% Convenience legends

legOrig = 'Original';
legHP   = 'High-pass (drift removed)';
legPR   = 'Physio-regression (cardiac+LFO removed)';
legOV   = 'Overprocessed (truth inverted)';
legTgt  = 'Target';

%% 2) Construct synthetic fNIRS signals

%% 2a) Construct two synthetic fNIRS channel observations

% True underlying coefficients for the responsive channel (Ch1)
w_true_ch1 = [0.60; -0.20; 0.10; 0.05];

% Non-responsive channel (Ch2)
w_true_ch2 = [0.00; -0.15; 0.08; 0.04];

% Construct signals
x1 = B * w_true_ch1 + 0.05 * randn(size(t));
x2 = B * w_true_ch2 + 0.05 * randn(size(t));

% Pack observations
P = [x1, x2];

%% Plot observations

offset = 1*max(abs(P(:)));

hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.88]);
hold on
plot(t,P(:,1),'b','LineWidth',opt.lineWidth);
plot(t,P(:,2)+offset,'r','LineWidth',opt.lineWidth)

title('Observations','FontSize',opt.fontSize)
xlabel('Time [s]','FontSize',opt.fontSize)
ylabel('HbO_2 [A.U.]','FontSize',opt.fontSize)

legend('Channel 1 (responsive)','Channel 2 (non-responsive)')

box on, grid on
set(gca,'FontSize',opt.fontSize);

mySaveFig(hFig,['..' filesep 'media' filesep ...
    'Overprocessing0007_OperationsAsPointShifts_Observations']);
%close(gcf);




%% 2b) Construct target fNIRS signals (hypothesis)

% Target for Channel 1 -> flat
% Target for Channel 2 -> clean HRF

w_target_ch1 = [0;0;0;0];
w_target_ch2 = [1;0;0;0];

q1 = B * w_target_ch1;
q2 = B * w_target_ch2;

Q = [q1,q2];



%% 3) Define processing operations as operators

% --- High-pass (remove drift)

w_HP_ch1 = w_true_ch1; w_HP_ch1(2)=0;
w_HP_ch2 = w_true_ch2; w_HP_ch2(2)=0;

qHP1 = B*w_HP_ch1;
qHP2 = B*w_HP_ch2;

Q_HP = [qHP1 qHP2];
A_HP = Q_HP * pinv(P);

% --- Physio regression (remove cardiac + LFO)

w_PR_ch1 = w_true_ch1; w_PR_ch1(3:4)=0;
w_PR_ch2 = w_true_ch2; w_PR_ch2(3:4)=0;

qPR1 = B*w_PR_ch1;
qPR2 = B*w_PR_ch2;

Q_PR = [qPR1 qPR2];
A_PR = Q_PR * pinv(P);

% --- Overprocessing (project to hypothesis)

A_over = Q * pinv(P);

%% 4) Apply operators

X_HP   = A_HP * P;
X_PR   = A_PR * P;
X_over = A_over * P;

%% 5) Plot time-domain transformations

offset = 1*max(abs(P(:)));

hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
tiledlayout(3,1)

nexttile
plot(t,P(:,1),'b','LineWidth',opt.lineWidth); hold on
plot(t,P(:,2)+offset,'r','LineWidth',opt.lineWidth)

plot(t,X_HP(:,1),'b--','LineWidth',opt.lineWidth)
plot(t,X_HP(:,2)+offset,'r--','LineWidth',opt.lineWidth);

title('High-pass filter (drift removed)','FontSize',opt.fontSize)

legend({'Ch1 original','Ch2 original','Ch1 processed','Ch2 processed'})
box on, grid on
set(gca,'FontSize',opt.fontSize);



nexttile
plot(t,P(:,1),'b','LineWidth',opt.lineWidth); hold on
plot(t,P(:,2)+offset,'r','LineWidth',opt.lineWidth)

plot(t,X_PR(:,1),'b--','LineWidth',opt.lineWidth)
plot(t,X_PR(:,2)+offset,'r--','LineWidth',opt.lineWidth);

title('Physiological regression','FontSize',opt.fontSize)

legend({'Ch1 original','Ch2 original','Ch1 processed','Ch2 processed'})
box on, grid on
set(gca,'FontSize',opt.fontSize);



nexttile
plot(t,P(:,1),'b','LineWidth',opt.lineWidth); hold on
plot(t,P(:,2)+offset,'r','LineWidth',opt.lineWidth)

plot(t,X_over(:,1),'b--','LineWidth',opt.lineWidth)
plot(t,X_over(:,2)+offset,'r--','LineWidth',opt.lineWidth)

title('Overprocessing','FontSize',opt.fontSize)

legend({'Ch1 original','Ch2 original','Ch1 processed','Ch2 processed'})
box on, grid on
set(gca,'FontSize',opt.fontSize);


mySaveFig(hFig,['..' filesep 'media' filesep ...
    'Overprocessing0007_OperationsAsPointShifts_HighPass']);
%close(gcf);




%% 6) Signal space representations

W_orig = pinv(B) * P;
W_HP   = pinv(B) * X_HP;
W_PR   = pinv(B) * X_PR;
W_OV   = pinv(B) * X_over;

% Coordinates
wHRF   = 1;
wDrift = 2;
wCard  = 3;
wLFO   = 4;

% 2D signal space plots
hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
tiledlayout(1,3)

nexttile
plotSignalSpace2D(gca,W_orig,W_HP,W_PR,W_OV,wDrift,wHRF,...
    'HRF vs Drift');

nexttile
plotSignalSpace2D(gca,W_orig,W_HP,W_PR,W_OV,wCard,wHRF,...
    'HRF vs Cardiac');

nexttile
plotSignalSpace2D(gca,W_orig,W_HP,W_PR,W_OV,wDrift,wCard,...
    'Drift vs Cardiac');

mySaveFig(hFig,['..' filesep 'media' filesep ...
    'Overprocessing0007_OperationsAsPointShifts_2DSignalSpace']);
%close(gcf);


%% 3D signal space plot
opt.markerSize = 150;

hFig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
hold on

% Channel markers
%mk1 = 'o'; % channel 1
%mk2 = 's'; % channel 2
mk = {'o','s','<','>','^','v','h','d'}; %Markers

% Original
scatter3(W_orig(1,1),W_orig(2,1),W_orig(3,1),...
        opt.markerSize,'b',mk{1},'filled');
scatter3(W_orig(1,2),W_orig(2,2),W_orig(3,2),...
        opt.markerSize,'r',mk{2},'filled');

% High-pass
scatter3(W_HP(1,1),W_HP(2,1),W_HP(3,1),...
        opt.markerSize,'b',mk{3},'filled');
scatter3(W_HP(1,2),W_HP(2,2),W_HP(3,2),...
        opt.markerSize,'r',mk{4},'filled');

% Physio regression
scatter3(W_PR(1,1),W_PR(2,1),W_PR(3,1),...
        opt.markerSize,'b',mk{5},'filled');
scatter3(W_PR(1,2),W_PR(2,2),W_PR(3,2),...
        opt.markerSize,'r',mk{6},'filled');

% Overprocessed
scatter3(W_OV(1,1),W_OV(2,1),W_OV(3,1),...
        opt.markerSize,'b',mk{7},'filled');
scatter3(W_OV(1,2),W_OV(2,2),W_OV(3,2),...
        opt.markerSize,'r',mk{8},'filled');

% Channel labels near original points
text(W_orig(1,1),W_orig(2,1),W_orig(3,1),'  Ch1',...
        'FontSize',opt.fontSize,'Color','b');
text(W_orig(1,2),W_orig(2,2),W_orig(3,2),'  Ch2',...
        'FontSize',opt.fontSize,'Color','r');

xlabel('HRF coefficient','FontSize',opt.fontSize)
ylabel('Drift coefficient','FontSize',opt.fontSize)
zlabel('Cardiac coefficient','FontSize',opt.fontSize)

title('Signal space (HRF–Drift–Cardiac)','FontSize',opt.fontSize)
%legend('Original','High-pass','Physiological regression','Overprocessed')
legend({'Orig Ch1','Orig Ch2',...
        'HP Ch1','HP Ch2',...
        'PR Ch1','PR Ch2',...
        'Overproc Ch1','Overproc Ch2'})

box on, grid on

view(35,25)
axis vis3d
camproj perspective

mySaveFig(hFig,['..' filesep 'media' filesep ...
    'Overprocessing0007_OperationsAsPointShifts_3DSignalSpace']);
%close(gcf);


end


%% AUXILIARY FUNCTIONS

function plotSignalSpace2D(ax,W_orig,W_HP,W_PR,W_OV,ix,iy,ttl)


opt.fontSize  = 20;
opt.lineWidth = 1.5;
opt.markerSize = 14;

% Channel markers
mk1 = 'o'; % channel 1
mk2 = 's'; % channel 2


hold on

% ----- Original
h1 = plot(ax,W_orig(ix,1),W_orig(iy,1),['k' mk1],...
    'MarkerFaceColor','k','MarkerSize',opt.markerSize);
h2 = plot(ax,W_orig(ix,2),W_orig(iy,2),['k' mk2],...
    'MarkerFaceColor','k','MarkerSize',opt.markerSize);

% ----- High-pass
h3 = plot(ax,W_HP(ix,1),W_HP(iy,1),['r' mk1],...
    'MarkerFaceColor','r','MarkerSize',opt.markerSize);
h4 = plot(ax,W_HP(ix,2),W_HP(iy,2),['r' mk2],...
    'MarkerFaceColor','r','MarkerSize',opt.markerSize);

% ----- Physio regression
h5 = plot(ax,W_PR(ix,1),W_PR(iy,1),['b' mk1],...
    'MarkerFaceColor','b','MarkerSize',opt.markerSize);
h6 = plot(ax,W_PR(ix,2),W_PR(iy,2),['b' mk2],...
    'MarkerFaceColor','b','MarkerSize',opt.markerSize);

% ----- Overprocessing
h7 = plot(ax,W_OV(ix,1),W_OV(iy,1),['g' mk1],...
    'MarkerFaceColor','g','MarkerSize',opt.markerSize);
h8 = plot(ax,W_OV(ix,2),W_OV(iy,2),['g' mk2],...
    'MarkerFaceColor','g','MarkerSize',opt.markerSize);

% ----- Processing paths (dashed lines) with start/end markers
plot(ax,[W_orig(ix,1) W_HP(ix,1)], [W_orig(iy,1) W_HP(iy,1)],'r-',...
    'LineWidth',opt.lineWidth)
plot(ax,[W_orig(ix,2) W_HP(ix,2)], [W_orig(iy,2) W_HP(iy,2)],'r--',...
    'LineWidth',opt.lineWidth)

plot(ax,[W_orig(ix,1) W_PR(ix,1)], [W_orig(iy,1) W_PR(iy,1)],'b-',...
    'LineWidth',opt.lineWidth)
plot(ax,[W_orig(ix,2) W_PR(ix,2)], [W_orig(iy,2) W_PR(iy,2)],'b--',...
    'LineWidth',opt.lineWidth)

plot(ax,[W_orig(ix,1) W_OV(ix,1)], [W_orig(iy,1) W_OV(iy,1)],'g-',...
    'LineWidth',opt.lineWidth)
plot(ax,[W_orig(ix,2) W_OV(ix,2)], [W_orig(iy,2) W_OV(iy,2)],'g--',...
    'LineWidth',opt.lineWidth)


xlabel(ax,['w_' num2str(ix)],'FontSize',opt.fontSize)
ylabel(ax,['w_' num2str(iy)],'FontSize',opt.fontSize)

title(ax,['Signal space: ' ttl],'FontSize',opt.fontSize)

% legend(ax,[h1 h2 h3 h4 h5 h6 h7 h8],...
%     {'Ch1 orig','Ch2 orig',...
%     'Ch1 HP','Ch2 HP',...
%     'Ch1 PR','Ch2 PR',...
%     'Ch1 Overproc.','Ch2 Overproc.'},...
%     'Location','best');
legend(ax,...
    {'Ch1 orig','Ch2 orig',...
     'Ch1 HP','Ch2 HP',...
     'Ch1 PR','Ch2 PR',...
     'Ch1 Overproc.','Ch2 Overproc.'},...
    'Location','best');

xlim(ax,padlim([W_orig(ix,:) W_HP(ix,:) W_PR(ix,:) W_OV(ix,:)]));
ylim(ax,padlim([W_orig(iy,:) W_HP(iy,:) W_PR(iy,:) W_OV(iy,:)]));
box on, grid on,
set(gca,'FontSize',opt.fontSize);

axis equal

end



% function drawArrow(hAxis,a,b,ix,iy,col)
% %Helper function
% h = quiver(hAxis,a(ix),a(iy), ...
%        b(ix)-a(ix), ...
%        b(iy)-a(iy), ...
%        0,...
%        'Color',col,...
%        'LineWidth',2,...
%        'MaxHeadSize',0.8,...
%        'AutoScale','off');
% 
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% end


function L = padlim(v)
% Pad limits a bit for nicer framing when plotting
v = v(~isnan(v) & ~isinf(v));
if isempty(v)
    L = [-1 1];
    return;
end
mn = min(v);
mx = max(v);
pad = 0.08*(mx-mn + eps);
if (pad == 0)
    pad = 0.1;
end
L = [mn-pad, mx+pad];
end