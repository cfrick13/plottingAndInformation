function plotInformation(fileDateName,dispName)
datename = fileDateName(1:21);
%determine the parent folder and directories
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b); %specifies folder in which all analysis is being done
    exportdir = strcat(parentdir,'Export'); %specifies where data is exported
    cd(exportdir);
    loaddir = strcat(parentdir,'Export'); %specifies where data is exported

    [~,b] = regexp(mdir,'/');
            if isempty(b)
                [~,b] = regexp(mdir,'\');
            end

    mfiledir =mdir(1:b(end)); %specifies location of matlab function file    

%load the file
    cd(mfiledir)
    flist = dir(strcat('fchighINFORMATIONzSELIMKHANOV*',datename,'*.mat'));
    fname = char(flist.name);
    load(fname);
    
    %load metadata associated with the experiment (requires manual input if there is ambiguity)
    FileName = fileDateName;
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
%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
    timeVec = metaData.timeVec;
    tvec = timeVec(1,:);
    selectedFeatures = [stimulationFrame-2:stimulationFrame+30];
    tm = round(tvec(selectedFeatures)-tvec(stimulationFrame));
    
    %plot for figure
%     infoInUnitsOfBitsMatrix(iter,cyclenumber,samplingcycle,selectedFeatures)
    info = squeeze(infoInUnitsOfBitsMatrix(:,:,1,:));
    for d3 = 1:size(info,3)
        i3 = squeeze(info(:,:,d3));
        mi3 = nanmean(i3,1); %mean of all iterations
        %determine input distribution at which info is max
        [imx,idx] = max(mi3,[],2);
        ebarstd(d3) = nanstd(i3(:,idx));
        imax(d3) = imx;
    end
    numberOfCells = length([PRofScellarray{:,1}]);
    nstr = ['n = ' num2str(numberOfCells)];

    
    if strcmp(datename,'2014_09_30 plate exp1')
%         lc = [0 0.4980 0];
        lc = [0 0.4470 0.7410];
    elseif strcmp(datename,'2017_01_30 plate exp2')
        lc = [0.6353 0.0784 0.1843];
    elseif strcmp(datename,'2017_03_13 plate exp4')
%         lc = [0 0.4470 0.7410];
        lc = [0 0.4980 0];
    end
    f=figure(44);
    subplot(2,1,1);
    lsm = length(selectedFeatures);
    tidx = 1:lsm;
    idx = 1:lsm;
    e = errorbar(tm(tidx),imax(idx),ebarstd(idx));hold on
    e.DisplayName = [dispName '(' nstr ')'];
    e.Color = lc;
%     plot(1:length(imax),imax);
    ylim([0 1.4]);
    title('Smad3 abundance channel capacity');
    xlabel('minutes after Tgfbeta addition');
    e.LineWidth = 1.5;
    xlim([min(tm) max(tm)]);
    legend('show')
    ylabel('channel capacity (bits)')
    
    subplot(2,1,2);
    lsm = length(selectedFeatures);
    tidx = 1:lsm;
    idx = [1:lsm]+lsm-1;
    e = errorbar(tm(tidx),imax(idx),ebarstd(idx));hold on
    e.DisplayName = [dispName '(' nstr ')'];
    e.LineWidth = 1.5;
    e.Color = lc;
%     plot(1:length(imax),imax);
    ylim([0 1.4]);
    xlim([min(tm) max(tm)]);
    title('Smad3 fold-change channel capacity');
    ylabel('channel capacity (bits)');
    xlabel('minutes after Tgfbeta addition');
    legend('show')
    f.Position = [680 85 785 893];

    
    %first conclusion is that iterations are not necessary if full sample size
    infoTwo = squeeze(infoInUnitsOfBitsMatrix(:,1,1,:));
    
    %iterations are necessary if sample size is reduced
    infoTwo = squeeze(infoInUnitsOfBitsMatrix(:,1,end,:));
    
    
%     figure(33)
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     plot((1./percentageofsamplesize),info);hold on
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     plot((1./percentageofsamplesize),info);
% 
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')
% 
%     ylim([0 2])
%     xlim([0 2])
%     xlabel('1/samplesize')
%     ylabel('estimated mutual information')
%     title('jacknife sampling without replacement to determine bias due to sample size')
end