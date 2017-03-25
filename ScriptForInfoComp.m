fileDateName = '2014_09_30 plate exp1_tracking_export.mat';
DetermineMutualInformationSelimkhanovMethodIterations_newparfor(fileDateName)
fileDateName = '2017_01_30 plate exp2_tracking_export.mat';
DetermineMutualInformationSelimkhanovMethodIterations_newparfor(fileDateName)
fileDateName = '2017_03_13 plate exp4_tracking_export.mat';
DetermineMutualInformationSelimkhanovMethodIterations_newparfor(fileDateName)
close all
fileDateName = '2014_09_30 plate exp1_tracking_export.mat';
dispName = 'exogenous cmv';
plotInformation(fileDateName,dispName)
fileDateName = '2017_01_30 plate exp2_tracking_export.mat';
dispName = 'endogenous';
plotInformation(fileDateName,dispName)
fileDateName = '2017_03_13 plate exp4_tracking_export.mat';
dispName = 'cmv endogenous';
plotInformation(fileDateName,dispName)

f=gcf;
for h=f.Children'
    if strcmp(h.Type,'legend')
        h.delete;
    end
end
f=gcf;
for h=f.Children'
l = legend(h,'show');
l.EdgeColor = 'w';
h.FontSize = 12;
h.FontName = 'helvetica';
h.Color = [0.97 0.97 0.97];
h.GridLineStyle = '--';
h.GridColor = 'k';
h.Box = 'off';
h.XGrid = 'on';
h.YGrid = 'on';
h.LineWidth = 1.5;
h.XColor = 'k';
h.YColor = 'k';
h.XLim = [0 80];
end
f.Color = [1 1 1];