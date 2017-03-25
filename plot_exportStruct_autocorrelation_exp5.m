function plot_exportStruct_autocorrelation_exp5

% close all
%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
%     [~,b ] = regexp(mdir,'/');
    [~,b] = regexp(mdir,'Tracking\w*/');
        if isempty(b)
%             [~,b] = regexp(mdir,'\');
            [~,b] = regexp(mdir,'Tracking\w*\');
        end
    parentdir = mdir(1:b);
    loaddir = strcat(parentdir,'Export');
    exportdir = strcat(parentdir,'LookingAtData');
cd(loaddir);

[~,b ] = regexp(mdir,'/');
    if isempty(b)
        [~,b] = regexp(mdir,'\');
    end
mfiledir =mdir(1:b(end));

exportdirz = exportdir;
%load the exported tracking structure
% FileName = uigetfile('*export.mat');%choose file to load
FileName = '2017_03_13 plate exp5_tracking_export.mat';
cd(loaddir)
load(FileName)


%load metadata associated with the experiment (requires manual input if
%there is ambiguity
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*metaData.mat');
cd(loaddir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    filename = uigetfile();
    metaData = load(filename);
end

timeVec = metaData.timeVec;

%load information regarding doses and scenes and tgfbeta addition
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*DoseAndScene*');
cd(loaddir)
filelist = dir(datequery);
    if isempty(filelist)
       dosestruct = makeDoseStruct; %run function to make doseStruct 
    else
        dosestructstruct = load(char(filelist.name));
        dosestruct = dosestructstruct.dosestruct;
    end
    
    %exportstruct
    %datastruct
    %dosestruct
    
    
    
    
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored 
colormapChoice = 'lines';
darkenFactor = 1.5;
    
%determine the scenes present in the experiment   
scenestr = 'scene';
sceneListArray = vertcat({exportStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosesstruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);


    dose = dosestruct(indicestwo).dose;
    frame = dosestruct(indicestwo).tgfFrame;
    
    dosestr = dosestruct(indicestwo).dosestr;
    framestr = dosestruct(indicestwo).tgfFramestr;
    
    
    [exportStruct(indices).dose] = deal(dose);
    [exportStruct(indices).frame] = deal(frame);
    [exportStruct(indices).dosestr] = deal(dosestr);
    [exportStruct(indices).framestr] = deal(framestr);
end

doseListArray = vertcat({exportStruct.dosestr});
doseList = unique(doseListArray);
    

    
%determine details needed for plotting such as when Tgfbeta is added, etc
%medianSmadbkg
stimulationFrame = exportStruct(1).frame;
smadTracesString = 'medianNucEGFP'; %value to plot
% smadTracesString = 'medianSmadbkg';
reporterTracesString = 'medianNucRFP';
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;




%establish the color map for plotting
coloringArray = vertcat({exportStruct.(coloringChoice)});
coloringList = unique(coloringArray);
indices = true(1,length(exportStruct));
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
    figure(1)
    cmap = colormap(parula(length(coloringList).*2));
    cmap = colormap(colormapChoice)./darkenFactor;
    close 1
    cmap = cmap; %darken the cmap



%assign a color array using the created colormap based on the choices above
colormapMatrix = zeros(length(coloringArrayTrunc),size(cmap,2));
for i=1:length(coloringArrayTrunc)
   cA = coloringArrayTrunc{i};
   idx = strcmp(uniqueColoring,cA);
   colormapMatrix(i,:) = cmap(idx,:);
%    colormapArray{i} = colorNames{idx};
   colormapArray{i} = cmap(idx,:);
end

%need to determine the number of scenes present and choose the time vector
%depending on the scene from which it was imaged
%THIS WORKS FOR NOW BUT NEEDS TO BE CHANGED
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);
    coloringArray = vertcat({exportStruct.(coloringChoice)});
    coloringList = unique(coloringArray);
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
for i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idx = strcmp(sceneListArray,sceneChoice);
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);
    
   stimulationFrame = dosestruct(idxtwo).tgfFrame;
   timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame);  
end
% setTequalZeroToStimulation = timeVector(stimulationFrame);
% xtickTimeVector = timeVector - setTequalZeroToStimulation;



indices = true(1,length(exportStruct));
%function to exract the cell traces, normalized and not
[smadCellTracesNorm,smadCellTraces] = extractTraces(exportStruct,indices,smadTracesString,finalFrame,stimulationFrame);
[reporterCellTracesNorm,reporterCellTraces] = extractTraces(exportStruct,indices,reporterTracesString,finalFrame,stimulationFrame);



plottingMat = smadCellTraces;
% plottingMat = reporterCellTraces;
% plottingMat = smadCellTracesNorm;
plottingMatNorm = smadCellTracesNorm;


    % conditionsArray = {'(s01|s02|s03|s04)','(s05|s06|s07|s08)','(s09|s10|s11|s12)','(s18|s19|s20|s21)'};
    conditionsArray = {'(s01|s02|s03|s04)'};
    doseList = {'0.07','0.0'}; %unstimulated vs stimulated cells exp5
%         doseList = {'0','2.4'}; %tgf stimulated cells exp4
    for condidx  = 1:length(doseList)

       
        f = figure(212);
        plottingMat = smadCellTraces; channelName = 'endo NG-Smad3';
        [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
%             sst =24;
%             ssp =26;
%             sortedpmat = sortedpmat(:,[1:sst ssp:size(sortedpmat,2)]);
%             timeMatrixForPlot = timeMatrixForPlot(:,[1:sst ssp:size(timeMatrixForPlot,2)]);
            tmat = timeMatrixForPlot./60;
            colormapMatrix = colormap(parula(size(sortedpmat,1)));
            subplot(2,3,1+(3*(condidx-1)));
                p=plot(tmat',sortedpmat','LineWidth',2);
    %             p=plot(sortedpmat','LineWidth',2);
                set(p, {'color'}, num2cell(colormapMatrix,2));
                title(channelName)
                xlim([(min(tmat(:))) max(tmat(:))])
                ylabel('NG-S3 nuc fluorescence');
                xlabel('hours')
        
        plottingMat = smadCellTracesNorm; channelName = 'fold change';
        [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
%             sortedpmat = sortedpmat(:,[1:sst ssp:size(sortedpmat,2)]);
%             timeMatrixForPlot = timeMatrixForPlot(:,[1:sst ssp:size(timeMatrixForPlot,2)]);
            colormapMatrix = colormap(parula(size(sortedpmat,1)));
            subplot(2,3,2+(3*(condidx-1)));
                p=plot(timeMatrixForPlot'./60,sortedpmat','LineWidth',2);
                set(p, {'color'}, num2cell(colormapMatrix,2));
                title(channelName)
                nvalue = size(sortedpmat,1);
                disp(nvalue)
                xlim([(min(tmat(:))) max(tmat(:))])
                ylabel('NG-S3 fold-change');
                xlabel('hours')
                
        plottingMat = reporterCellTracesNorm; channelName = 'reporter fold change';
        [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
%             sortedpmat = sortedpmat(:,[1:sst ssp:size(sortedpmat,2)]);
%             timeMatrixForPlot = timeMatrixForPlot(:,[1:sst ssp:size(timeMatrixForPlot,2)]);
            colormapMatrix = colormap(parula(size(sortedpmat,1)));
            subplot(2,3,3+(3*(condidx-1)));
                p=plot(timeMatrixForPlot'./60,sortedpmat','LineWidth',2); hold on
                set(p, {'color'}, num2cell(colormapMatrix,2));
                title(channelName)
                nvalue = size(sortedpmat,1);
                disp(nvalue)
                xlim([(min(tmat(:))) max(tmat(:))])
                ylim([0 10])
                ylabel('caga12 mCherry fold-change');
                xlabel('hours')
        
%         subplot(1,2,1);plot(sortedpmatNorm');
%         ylim([0 10])
%         xlim([0 100])
%         title('smad')
%         subplot(1,2,2);plot(reportersortedpmatNorm');
%         ylim([0 10])
%         xlim([0 100])
%         title('reporter')




        plottingMat = smadCellTraces; channelName = 'endo NG-Smad3';
        [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
%             only perform auto correlation calculation with same cells
%             throughout the timelapse imaging
            invec = isnan(sortedpmat);
            sinvec = sum(~invec,1);   
            numindx = find(sinvec(:,stimulationFrame:end)<max(sinvec).*0.70,1,'first')+stimulationFrame;
            spvec = sortedpmat(:,numindx);
            spidx = isnan(spvec);
            
            endframe =numindx;
            startframe = stimulationFrame;
            startframe = stimulationFrame+15;%jump ahead 1 hr
%             startframe = stimulationFrame+30;%jump ahead 2 hrs
            sortedpmat = sortedpmat(~spidx,startframe:numindx);
            tmatforplot = timeMatrixForPlot(~spidx,startframe:numindx);

            Ti = tmatforplot;
            Fi = sortedpmat;
            sortFrame = 1;
[timeVec,alphaAuto,alphaAutoErrUp,alphaAutoErrDown,Xi,Ri,Ti,Frankedi] = calcAutoCorrAustin(Ti,Fi,sortFrame);        

        ff = figure(0111)
        subplot(2,2,1+(2*(condidx-1)));
            tmat = tmatforplot./60;
            pmat = Frankedi;
                p = plot(tmat',pmat','LineWidth',2);
                    colormapMatrix = colormap(parula(size(pmat,1)));
                    set(p, {'color'}, num2cell(colormapMatrix,2));
                    title(channelName)
                    ylabel('NG-S3 nuc fluorescence');
                    xlabel('hours')
                    nvalue = size(pmat,1);
                    nvaluestr = strcat('n = ',num2str(nvalue));
                    t = text(0,0,nvaluestr);
                        t.Units = 'normalized';
                        t.HorizontalAlignment = 'right';
                        t.Position = [0.95 0.9];
                    dstr = strcat(doseList{condidx});
                    t = text(0,0,dstr);
                        t.Units = 'normalized';
                        t.HorizontalAlignment = 'right';
                        t.Position = [0.95 0.85];
                        
                    xlim([(min(tmat(:))) max(tmat(:))])
        subplot(2,2,2+(2*(condidx-1)));
            tmat = Ti./60;
            pmat = Ri;
                p = plot(tmat',pmat','LineWidth',2);
                    colormapMatrix = colormap(parula(size(pmat,1)));
                    set(p, {'color'}, num2cell(colormapMatrix,2));
                    title('Ranks')
                    axis('tight')
                    xlim([(min(tmat(:))) max(tmat(:))])
                    ylabel('ranks');
                    xlabel('hours')
        figure(555)
            tmat = timeVec./60;
            pmat = alphaAuto;
                pp(condidx) = plot(tmat',pmat','LineWidth',2);hold on
                    title('AutoCorrelation')
                    xlim([(min(tmat(:))) max(tmat(:))])
                    pp(condidx).DisplayName = strcat(channelName,'-',dstr);
%                     lgd = legend(doseList');
                    ylim([0 1])
                    ylabel('auto correlation')
                    xlabel('hours')



        f.Position = [100 100 1800 500];
        ff.Position = [100 100 1000 500];

    end
    
    for h = f.Children'
        h.FontSize = 10;
        h.FontName = 'helvetica';
    end
end


function channelinputs =channelregexpmaker(channelstoinput)
    channelinputs = '(';
    for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
        if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i});
        elseif i < length(channelstoinput)
            channelinputs = strcat(channelinputs,'|',channelstoinput{i});
        else
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
        end
    end
end


function [cellTracesNorm,cellTraces] = extractTraces(exportStruct,indices,xTracesString,finalFrame,stimulationFrame)
% extract the cell traces for the desired number of frames
cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
cellTraces = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]

% normalize by basal values
basalLength = 3;
if (stimulationFrame-basalLength)<1
    basalLength=0;
end

basalVector = nanmedian(cellTraces(:,stimulationFrame-basalLength:stimulationFrame),2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');
end


function [timeVec,alphaAuto,alphaAutoErrUp,alphaAutoErrDown,Xi,Ri,Ti,Frankedi] = calcAutoCorrAustin(Ti,Fi,sortFrame)
%determine auto correlation based on method in Sigal et al 2006 paper
%method to calculate auto correlation is from austin et al 2006 Gene
%network shaping of inherent noise spectra (Science). Gives identical
%results to Sigal et al
alphaAutoErrUp=[];
alphaAutoErrDown=[];
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
[Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame);



N = traceLength;
M = numberOfCells;
alphaAuto = nan(1,N);
for tau = 1:N
%adjust tau
    j=tau-1;
%determine denominator of autocorrelation
    numerato = nan(1,M);
    denominato = nan(1,M);
    for m = 1:M
        
        %determine numerator;
            numer = nan(1,N-j);
            for n = 1:N-j
                numer(1,n) = Xi(m,n).*Xi(m,n+j);
            end
            numerato(m) = nanmean(numer);
        
        %determine denominator
            denom = nan(1,N);
            for n = 1:N
                denom(1,n) = Xi(m,n).^2;
            end
            denominato(m) = nanmean(denom);
        
    end
    %compute numerator and denominator
    numerator = nanmean(numerato);
    denominator = nanmean(denominato);
    
%determine autocorrelation
    alphaAuto(tau) = numerator./denominator;
end

    timeVec = Ti(1,:);
    
end

function [timeVec,alphaAuto,alphaAutoErrUp,alphaAutoErrDown,Xi,Ri,Ti,Frankedi] = calcAutoCorrAustinTZero(Ti,Fi,sortFrame)
%determine auto correlation based on method in milo et al 2006 paper
alphaAutoErrUp=[];
alphaAutoErrDown=[];
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
[Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame);



N = traceLength;
M = numberOfCells;
alphaAuto = nan(1,N);
for tau = 1:N
%adjust tau
    j=tau-1;
%determine denominator of autocorrelation
    numerato = nan(1,M);
    denominato = nan(1,M);
    for m = 1:M
        
        %determine numerator;
            numer = nan(1,N-j);
            for n = 1:1
                numer(1,n) = Xi(m,n).*Xi(m,n+j);
            end
            numerato(m) = nanmean(numer);
        
        %determine denominator
            denom = nan(1,N);
            for n = 1:1
                denom(1,n) = Xi(m,n).^2;
            end
            denominato(m) = nanmean(denom);
        
    end
    %compute numerator and denominator
    numerator = nanmean(numerato);
    denominator = nanmean(denominato);
    
%determine autocorrelation
    alphaAuto(tau) = numerator./denominator;
end

    timeVec = Ti(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting
    figure(995)
    colormapMatrix = colormap(parula(size(Xi,1)));
    p = plot(Ti'./60,Xi','LineWidth',1.5,'Color',[0.5 0.1 0.1]);hold on
    set(p, {'color'}, num2cell(colormapMatrix,2));
    xlabel('hours');
            ylabel('total nuclear fluorescence (au)')
            title('Level of endogenous nuclear NG-Smad3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame)
removeNaN=1;
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
sortDim = 1; %or 1?

%first sort your poulation based on TgfFrame
    [~,I] = sort(Fi(:,sortFrame),sortDim);

%next build your sorted fluorescence matrix
    Frankedi = zeros(size(Fi));
    for i = 1:numberOfCells
        Frankedi(i,:) = Fi(I(i),:);
    end

%next determine the rank of each trace at a given point
    Ri = zeros(size(Frankedi));
    for j = 1:traceLength
        FrankediFrame = Frankedi(:,j);
        [~,I] = sort(Frankedi(:,j),sortDim); %I(1) is the index of the highest value (not the rank)
        rankI = zeros(size(FrankediFrame));
            rankI(I) = 1:length(I);
            if removeNaN==1
            rankI(I(isnan(FrankediFrame(I)))) = NaN; %remove all NaN values;
            end
            Ri(:,j) = rankI;
    end


%next determine the relative ranks of cells
    Xi = zeros(size(Ri));
    for i = 1:numberOfCells
        for j = 1:traceLength
            Xi(:,j) = Ri(:,j) - nanmean(Ri(:,j));
        end
    end

end


function [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix)
        doseArray = [exportStruct.dose];
        dlog = (doseArray == str2double(doseList{condidx}));
        idx = dlog;
        plottingMatforSort = plottingMat(idx,:);
        plotmatframesorting = plottingMatforSort(:,stimulationFrame);
        [~,indsort] = sort(plotmatframesorting);
        sortedpmat = zeros(size(plottingMatforSort));
        for inin = 1:length(indsort)
            sortedpmat(inin,:) = plottingMatforSort(indsort(inin),:);
        end
        timeMatrixForPlot = timeMatrix(idx,:);
end