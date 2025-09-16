function Overprocessing0004_RemasterOHBM2010
%Overprocessing0004_RemasterOHBM2010
%
% Remasters the OHBM 2010 Fig 2.
%
% Copyright 2025
% @author Felipe Orihuela-Espina
%
% See also 
%


%% Log
%
% 8-Sep-2025: FOE
%   + File created.
%


srcDir = ['..' filesep 'media' filesep];
opt.fontSize  = 18;
opt.lineWidth = 1.5;


%Figure 2 - Subfigure 1
%Experiment: VM; Session = Moto; Stimulus = A
filename = 'OHBM2010_0003_NumActivationChannels_ExpVM_SessMoto_StimA';
hFig = openfig([srcDir filename '.fig']);
hAxis = findall(hFig,'type','axes');
%Beware! The first axis is the suptitle.
nAxis = numel(hAxis);

set(hFig,'Units','normalized','Position',[0.05 0.05 0.9 0.88]);
set(hAxis(1),'FontSize',opt.fontSize+2);
set(hAxis(2:end),'FontSize',opt.fontSize);
xlabel(hAxis(2:end),'Task [s]');
ylabel(hAxis(2:end),'Baseline [s]');

tmpPos = get(hAxis(2),'XTick');
tmpStr = get(hAxis(2),'XTickLabels');
tmpPos = tmpPos(1:4:end);
tmpStr = tmpStr(1:4:end);
set(hAxis(2:end),'XTick',tmpPos);
set(hAxis(2:end),'XTickLabels',tmpStr);

tmpPos = get(hAxis(2),'YTick');
tmpStr = get(hAxis(2),'YTickLabels');
tmpPos = tmpPos(1:3:end);
tmpStr = tmpStr(1:3:end);
set(hAxis(2:end),'YTick',tmpPos);
set(hAxis(2:end),'YTickLabels',tmpStr);

for iAxis = 2:nAxis
    hAxis(iAxis).XAxis.Label.Rotation = 20;
    hAxis(iAxis).YAxis.Label.Rotation = -35;
end


set(gcf,"PaperPositionMode","auto")
mySaveFig(hFig,[srcDir filename '_Remastered']);
close(gcf);




%Figure 2 - Subfigure 2
%Experiment: DP; Session = SM; Stimulus = Shadow
filename = 'OHBM2010_0003_NumActivationChannels_ExpDP_SessSM_StimShadow';

hFig = openfig([srcDir filename '.fig']);
hAxis = findall(hFig,'type','axes');
%Beware! The first axis is the suptitle.
nAxis = numel(hAxis);

set(hFig,'Units','normalized','Position',[0.05 0.05 0.9 0.88]);
set(hAxis(1),'FontSize',opt.fontSize+2);
set(hAxis(2:end),'FontSize',opt.fontSize);
xlabel(hAxis(2:end),'Task [s]');
ylabel(hAxis(2:end),'Baseline [s]');

tmpPos = get(hAxis(2),'XTick');
tmpStr = get(hAxis(2),'XTickLabels');
tmpPos = tmpPos(1:4:end);
tmpStr = tmpStr(1:4:end);
set(hAxis(2:end),'XTick',tmpPos);
set(hAxis(2:end),'XTickLabels',tmpStr);

tmpPos = get(hAxis(2),'YTick');
tmpStr = get(hAxis(2),'YTickLabels');
tmpPos = tmpPos(1:3:end);
tmpStr = tmpStr(1:3:end);
set(hAxis(2:end),'YTick',tmpPos);
set(hAxis(2:end),'YTickLabels',tmpStr);

for iAxis = 2:nAxis
    hAxis(iAxis).XAxis.Label.Rotation = 20;
    hAxis(iAxis).YAxis.Label.Rotation = -35;
end


set(gcf,"PaperPositionMode","auto")
mySaveFig(hFig,[srcDir filename '_Remastered']);
close(gcf);





end
