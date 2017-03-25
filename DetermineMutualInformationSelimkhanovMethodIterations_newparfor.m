%DetermineMutualInformation(SigFeaturesArray,k,selectedFeatures,dimensions)
% see Estimating Mutual Information by local gaussian approximation by
% Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger (2008)
%   "SigFeaturesArray" is a cell array where each cell is a different condition and
%   the items contained within the cells are matrices with the experimental
%   measurement being made (such as a measure of signaling response) in
%   (dim=2). (dim=1) corresponds to the cell (or obseration# where the
%   measurement was made)
%   
%  SigFeaturesArray = {1,i} where n = conditions (e.g. doses)
%  SigFeaturesArray{1,i} = [n,m] where each m is a response type (measure) for n
%  number of cells (observation)
%     
%     the matrix within SigFeaturesArray is assemble as shown below
%  [measure1|observation1  measure2|observation1 ..... measureN|observation1;...
%          :                         :                           :          
%          :                         :                           :
%  [measure1|observationM  measure2|observationM ..... measureN|observationM]                                     :


%   "k" is the number of nearest neighbors
%   "selectedFeatures" defines the indices of the specific signaling features within "SigFeaturesArray" to be
%   quantified. Should be a vector of scalars.
%   "dimensions" is the dimensions to determine the information at (1d ([x]), 2d ([x,y]), or 3d ([x,y,z]) timepoint
%   combinations. Should be a scalar : 1, 2 or 3.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('C:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled\')
% load('preinfo-2016-02-29.mat')
% load('preinfo-2016-03-03.mat')
% % cd('C:\Users\zeiss\Documents\MATLAB\')
% load('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/FIGURES/2016_02_29/Figure3/INFORMATIONz01-Mar-2016-1d.mat')
% 
% cycle=1;
% % fnames = sort(fieldnames(SCALAR))
% for i = 1:length(SCALAR)
% sc = SCALAR{i};
% if size(sc,2)>50
% INFOyo{cycle} = sc(1:end-1,:);
% cycle=cycle+1;
% end
% end
% selectedFeaturesoned = [[1:21] [65:78] [93:106]];
% Chris_Info_CalculationsUpdatedALLDMACnew(INFOyo, 3, selectedFeaturesoned,1)

% load('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/FIGURES/2016_02_29/Figure3/INFORMATIONz01-Mar-2016-1d.mat')
% cycle=1;
% % fnames = sort(fieldnames(SCALAR))
% for i = 1:length(SCALAR)
% sc = SCALAR{i};
% if size(sc,2)>50
% INFOyo{cycle} = sc(1:end-1,:);
% cycle=cycle+1;
% end
% end
% selectedFeaturesoned = [[1:21] [65:78] [93:106]];
% SigFeaturesArray = SCALAR;
% selectedFeatures = [65:78];
% k=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = DetermineMutualInformationSelimkhanovMethodIterations_newparfor(fileDateName)
k=3;
dimensions=1;
datename = fileDateName(1:21);
%determine the location of the matlab function and establish export 
%directory in relation to that filepath
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b); %specifies folder in which all analysis is being done
    exportdir = strcat(parentdir,'Export'); %specifies where data is exported
    cd(exportdir);

    [~,b] = regexp(mdir,'/');
            if isempty(b)
                [~,b] = regexp(mdir,'\');
            end

    mfiledir =mdir(1:b(end)); %specifies location of matlab function file

%loads the specific file you want

    [SigFeaturesArray,dosestruct] = convert_exportStruct_doses_toInformation(fileDateName);
    stimulationFrame = dosestruct(1).tgfFrame;
    pmatdim2 = dosestruct(1).pmatdim2;
    numTpoints = 8;
%determine which frames to get info from
    featureVec = round(linspace(stimulationFrame,stimulationFrame+30,numTpoints)-2);

    %BEST
%     selectedFeatures = [stimulationFrame-2:stimulationFrame+30];
%     selectedFeatures = [selectedFeatures selectedFeatures+pmatdim2];

    %NEXT BEST
    selectedFeatures = [stimulationFrame-1:stimulationFrame+11];
    selectedFeatures = [selectedFeatures selectedFeatures+pmatdim2];
    
    dimensions=1;

    nInputSignalConditionsInStruct = length(SigFeaturesArray);
    featuresVector = selectedFeatures;
    assembledFeatures = chooseQFfunct(SigFeaturesArray,featuresVector,dimensions); %assembles an array of 1d ([x]), 2d ([x,y]), or 3d ([x,y,z]) timepoint combinations depedning on the dimension input
    lengthOfFeatures = length(assembledFeatures);    


%define the smoothness and the means (mus) for the signaling input
%probability distributions
    smthness=2000;
    number_of_PofS=50; %number of probability distributions to test
    lowsigma=0.0001;
    highsigma=1;
    InputSignals = generateInputDistributions(nInputSignalConditionsInStruct,smthness,number_of_PofS,lowsigma,highsigma); %generates a continuum of unimodal, bimodal, trimodal inputs
    %%%%%%%%%%%
    InputSignals(1,:) = ones(1,size(InputSignals,2))./size(InputSignals,2); %set the first input to be a uniform input distribution. 
    InputSignals(1,:) = InputSignals(50,:); %set the first input to be a uniform input distribution. 

    %%%%%%%%%%%


% drawnow
disp(datename)
    iterations = 30;
    percentageofsamplesize = 0.6 : 0.1 : 0.7;
    infoInUnitsOfBitsMatrix = zeros(iterations,size(InputSignals,1),length(percentageofsamplesize),lengthOfFeatures);
    for iter = 1:iterations
        %Define Data
        DataCell = cell(nInputSignalConditionsInStruct,length(assembledFeatures)); %initialize
        PRofScellarray = makePRofScellarray(DataCell,SigFeaturesArray,assembledFeatures);  %and then determine PRofScellarray
        numberOfSignals = size(PRofScellarray,1);


        for samplingcycle = 1:length(percentageofsamplesize)
            disp([iter samplingcycle]); %Display cycle to track progress
            samplesizepercentage = percentageofsamplesize(samplingcycle);
            m = numberOfSignals;
            responsesGivenSvectorarray = cell(m,lengthOfFeatures);
                for ff=1:lengthOfFeatures
                    for i = 1:m 
                        ResponsesGivenS = PRofScellarray{i,ff};
                        responsesGivenSvector =ResponsesGivenS';
                        if sum(sum(isnan(responsesGivenSvector)))>0
                            idx = reshape(isnan(ResponsesGivenS),size(ResponsesGivenS));
                            [~,col,~] = find(idx==1);
                            [~,dim2] = size(ResponsesGivenS);
                            nline = 1:dim2;
                            nline(col)=[];
                            responsesGivenSvector = ResponsesGivenS(:,nline)';  
                        end %reshape responseGivenSvector


                        sampledResponsesGivenSvector = datasample(responsesGivenSvector,round(length(responsesGivenSvector).*samplesizepercentage),'Replace',false);
            %             bootresponsesGivenSvector = bootstrp(1,@(x) x,responsesGivenSvector); %cannot perform bootstrap because you end up with duplicate values and an information of -INF
                        responsesGivenSvectorarray{i,ff} =sampledResponsesGivenSvector;

                    end
                end
                iyo = zeros(size(InputSignals,1),lengthOfFeatures);
                parfor cyclenumber = 1:size(InputSignals,1)
                    q = InputSignals(cyclenumber,:); %vector where each number is a probability of a give signal (condition) 
                    

                    %conditional entropy, H(R|S)
                    entropydiffRS = calculateEntropyRS(responsesGivenSvectorarray,numberOfSignals,q,lengthOfFeatures,k,assembledFeatures,dimensions);
                    %unconditional entropy, H(R)
                    entropydiffR = calculateEntropyR(responsesGivenSvectorarray,numberOfSignals,q,lengthOfFeatures,k,assembledFeatures,dimensions);

                    iyo(cyclenumber,:) = entropydiffR - entropydiffRS;
                end
                infoInUnitsOfBitsMatrix(iter,:,samplingcycle,:) = iyo;
                
        end %samplecycle loop
        stophere=1;
    end
    stophere=1;


    
    cd(mfiledir)
    save(strcat('fchighINFORMATIONzSELIMKHANOV',datename,'-',num2str(dimensions),'-',num2str(iterations),'.mat'),'-v7.3');
    
    
    
    %plot for figure
    figure(33)
    info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
    plot((1./percentageofsamplesize),info);hold on
    info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
    plot((1./percentageofsamplesize),info);

    info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
    infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
    errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')
    info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
    infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
    errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')

    ylim([0 2])
    xlim([0 2])
    xlabel('1/samplesize')
    ylabel('estimated mutual information')
    title('jacknife sampling without replacement to determine bias due to sample size')



    % save the data from the information run
    % cd('/Users/frick/Documents/MATLAB')
%     cd(exportdir)
%     save(strcat('INFORMATIONzSELIMKHANOV',date,'-',num2str(dimensions),'-',num2str(iterations),'d.mat'));
end



function DataCell = makePRofScellarray(DataCell,SigFeaturesArray,assembledFeatures)
%construct a cellarray of PRofS where dim1 is the signal number and dim2 is
%the feature of interest to analyze (such as cellular response)
for i = 1:size(DataCell,1)
    for j=1:length(assembledFeatures)
    DataI = SigFeaturesArray{i};
    CH = DataI(assembledFeatures{j},:);
    DataCell{i,j} = CH;
    end
end
end

function entropydiffRS = calculateEntropyRS(responsesGivenSvectorarray,numberOfSignals,q,lengthOfFeatures,k,assembledFeatures,dimensions)


m = numberOfSignals;    
Epsarray = cell(m,lengthOfFeatures);      
for ff=1:lengthOfFeatures
    for i = 1:m 
        responseSGivenSvector = responsesGivenSvectorarray{i,ff};
        [~, Eps] = knnsearch(responseSGivenSvector,responseSGivenSvector,'K', k,'Distance','euclidean');
        Eps = Eps(:,k);%Eps is the distance x 2
        Epsarray{i,ff} = Eps;
    end
end


entropydiffRS = zeros(1,lengthOfFeatures);
for ff=1:lengthOfFeatures

        featuresSelected = assembledFeatures{ff};
        if dimensions >1 
            dmnsn = determineDimension(featuresSelected,responseSGivenSvector);
        else %else statement is for speed
            dmnsn=1;
        end

        entropydiffRSfori = zeros(1,m);
        for i = 1:m 
            responseSGivenSvector = responsesGivenSvectorarray{i,ff};

            Nnumi = length(responseSGivenSvector);
            Eps = Epsarray{i,ff};

            log2val = zeros(1,Nnumi);
            for j=1:Nnumi
                Epsj  = Eps(j); %twice the distance
                Vd = (pi.^(dmnsn./2))./gamma((dmnsn./2)+1);
                log2arg = k./(Nnumi.*Vd.*(Epsj.^dmnsn));
                log2val(j) = log2(log2arg);
            end
            
            entropydiffRSfori(i)  = (q(i)./Nnumi).*(sum(log2val)); %in equation 2.12 ni = Nx (from 2.1) for a given signal qi
        end

    entropydiffRS(ff) = -sum(entropydiffRSfori);
end

end


function entropydiffR = calculateEntropyR(responsesGivenSvectorarray,numberOfSignals,q,lengthOfFeatures,k,assembledFeatures,dimensions)

m = numberOfSignals;

    
for ff=1:lengthOfFeatures
    for i = 1:m 
        responseSGivenSvector_i = responsesGivenSvectorarray{i,ff};
        for w=1:m
            responsesGivenSvector_w = responsesGivenSvectorarray{w,ff};    
            [~, Epsj] = knnsearch(responsesGivenSvector_w,responseSGivenSvector_i,'K', k,'Distance','euclidean');
            Epsj = Epsj(:,k);%Eps is the distance x 2
            Epsarray{i,w,ff} = Epsj;
        end
    end
end
        
for ff=1:lengthOfFeatures
    entropydiffRfori = zeros(1,m);
    for i = 1:m 
        
        featuresSelected = assembledFeatures{ff};
        
        if dimensions >1  %determine dimension
            dmnsn = determineDimension(featuresSelected,responseSGivenSvector_i);
        else %else statement is for speed
            dmnsn=1;
        end



        Vd = (pi.^(dmnsn./2))./gamma((dmnsn./2)+1);     
        responseSGivenSvector_i = responsesGivenSvectorarray{i,ff};
        Nnumi = length(responseSGivenSvector_i);
        
            log2val = zeros(1,Nnumi);
            for j = 1:Nnumi
                log2arg = zeros(1,w);
                for w=1:m
                    responsesGivenSvector_w = responsesGivenSvectorarray{w,ff};
                    Nnumw = length(responsesGivenSvector_w);
                    Epsvector = Epsarray{i,w,ff};
                    Epsj = Epsvector(j);
                    log2arg(w) = q(w).*(k./(Nnumw.*Vd.*(Epsj.^dmnsn)));
                end

                log2argsum = sum(log2arg); %sum from w to m
                log2val(j) = log2(log2argsum);
            end
        
        log2valsum = sum(log2val); %sum from j to Nnumi
        entropydiffRfori(i) = (q(i)./(Nnumi)).*log2valsum;
    end
        
    entropydiffR(ff) = -sum(entropydiffRfori); %sum from i to m
end

    stophere=1;
end





















function SamSigs = generateInputDistributions(scalarLength,smthness,number_of_PofS,lowsigma,highsigma)
Asigma = linspace(lowsigma,highsigma,number_of_PofS./2);
CDsigma = sort(Asigma,'descend');
Asigma = horzcat(Asigma,CDsigma); %lowsigma->highsigma->lowsigma with length of number_of_PofS
CDsigma = horzcat(CDsigma,ones(size(Asigma)).*lowsigma); %highsigma -> lowsigma, for half of length of number_of_PofS and then lowsigma->lowsigma for rest of lenght of number_of_PofS
%by altering the sigma we can get wide or narrow distributions so that when
%we run the script below, the dominant distributions will either be
%centered at 1 (unimodal), 0 and 2 (bimodal) or 0,1,and,2 (trimodal).
SamSigs = zeros(number_of_PofS,scalarLength);
    for i=1:length(Asigma)
        A = normrnd(1,Asigma(i),smthness,1);  %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=1 and sigma = Asigma(i)
        C = normrnd(0,CDsigma(i),smthness,1); %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=0 and sigma = CDsigma(i)
        D = normrnd(2,CDsigma(i),smthness,1); %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=2 and sigma = CDsigma(i)
        E = vertcat(A,C,D); %assembles A,C,D into a long vector
        h=histogram(E); 
        h.Normalization = 'probability';
        h.BinEdges = linspace(0,2,scalarLength+1);
        h.Normalization = 'probability';
        SamSigs(i,:) = h.Values; %take the values from the histogram 
    end
f=h.Parent;
ff=f.Parent;
ff.delete;
imagesc(SamSigs) %displays image of the assembled input distributions. These are P(S). where S is the signal (or condition)
end

function dmnsn = determineDimension(featuresSelected,responsesGivenSvector)
[~,epsdim] = knnsearch(featuresSelected',featuresSelected','K',2);    
    if epsdim == 0  
        dmnsn = 1;
    else
        epsdim = epsdim(:,2);
        dimd = sum(epsdim==0);
        if dimd ==3
            dmnsn = 1;
        elseif dimd ==2
            dmnsn = 2;
        else 
            dmnsn = size(responsesGivenSvector,2); %dimension is either 1D or 2D
        end
    end
end

function chosenScalars = chooseQFfunct(SCALAR,featuresVector,dimensions)
cycle = 1;
if isempty(featuresVector)   
    featuresVector = 1:size(SCALAR{1},1);   

    if dimensions ==1
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
        chosenScalars{cycle} = [i];
        cycle=cycle+1;
        end  
    end
    
    if dimensions ==2
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
            chosenScalars{cycle} = [i j];
            cycle=cycle+1;
            end
        end
    end

    if dimensions ==3
    chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
                for kk =featuresVector
                chosenScalars{cycle} = [i j kk];
                cycle=cycle+1;
                end
            end
        end 
    end
      

else

    if dimensions ==1
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
        chosenScalars{cycle} = [i];
        cycle=cycle+1;
        end
    end


    if dimensions ==2
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
            chosenScalars{cycle} = [i j];
            cycle=cycle+1;
            end
        end
    end

    if dimensions ==3
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
                for kk =featuresVector
                chosenScalars{cycle} = [i j kk];
                cycle=cycle+1;
                end
            end
        end
    end


end
end