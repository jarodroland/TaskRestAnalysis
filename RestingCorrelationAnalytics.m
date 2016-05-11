% RestingCorrelationAnalytics.m

% ptNumber = 145;
subjectID = num2str(ptNumber);
dataDir = ['E:\Data\ECoG Task-Rest\' subjectID '\'];
restDataFile = [dataDir 'Rest\RestingCorrelations.mat'];

load(restDataFile);

numFreqBands = length(blpCorrelations);
numChannels = size(blpCorrelations(1).restingCorrelationMatrix, 2);
if(numChannels > 64)
    numChannels = 64;
end

for freqBand = 1:numFreqBands

    corrMatrix = blpCorrelations(freqBand).restingCorrelationMatrix;
    corrMatrix = corrMatrix(1:numChannels, 1:numChannels);

    % all pair-wise electrode-electrode correlation values
    indices = find( triu(corrMatrix, 1) );
    corrValues = corrMatrix(indices );
    corrValuesSorted = sort(corrValues);
    figAllChannels = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
    hold on;
    plot(corrValuesSorted, '.k')
    axis('fill')
    plot([0 length(corrValuesSorted)], [0 0], '-k');
    xlim([0 length(corrValuesSorted)]);
    set(gca,'units','normalized','position',[0.02, 0.05, 0.96, 0.9])
    title(['Subject ' subjectID ', Freq ' num2str(blpCorrelations(freqBand).bandPassFrequencies(1)) '-' num2str(blpCorrelations(freqBand).bandPassFrequencies(2)) 'Hz: All pair-wise electrode cross-correlation values sorted']);
    
    % per channel correlation values
    chanCorrValuesMax = zeros(1, numChannels);
    chanCorrValuesMin = zeros(1, numChannels);
    chanCorrValuesMean = zeros(1, numChannels);
    chanCorrValuesStd = zeros(1, numChannels);
    for chan = 1:numChannels
        chanCorrValues = corrMatrix(chan, :);
        chanCorrValues(chan) = [];
        chanCorrValues(isnan(chanCorrValues)) = 0;
        chanCorrValuesMax(chan) = max(chanCorrValues);
        chanCorrValuesMin(chan) = min(chanCorrValues);
        chanCorrValuesMean(chan) = mean(chanCorrValues);
        chanCorrValuesStd(chan) = std(chanCorrValues);
    end
    
    figPerChannel = figure('units','normalized','outerposition',[0 0.2 1 0.3]);
    hold on;
    errorbar(chanCorrValuesMean, chanCorrValuesStd, 'xk');
    plot(chanCorrValuesMax, '^k');
    plot(chanCorrValuesMin, 'vk');
    xlim([0 numChannels+1]);
    set(gca,'units','normalized','position',[0.02, 0.1, 0.96, 0.8])
    title(['Subject ' subjectID ', Freq ' num2str(blpCorrelations(freqBand).bandPassFrequencies(1)) '-' num2str(blpCorrelations(freqBand).bandPassFrequencies(2)) 'Hz: Mean, STD, and Min/Max cross-correlations per electrode']);
    
    % save figure files
    fileOut = sprintf('%sFigures\\Cross Correlation Values\\CrossCorr-%s-AllCh-Freq%03i-%03iHz', dataDir, subjectID, blpCorrelations(freqBand).bandPassFrequencies(1), blpCorrelations(freqBand).bandPassFrequencies(2));
%     print(figAllChannels, [fileOut '.svg'],  '-dsvg');     % save figure as SVG
    print(figAllChannels, [fileOut '.png'],  '-dpng', '-r300');      % save figure as PNG
    savefig(figAllChannels, [fileOut '.fig']);              % save figure as FIG
    close(figAllChannels);
    fileOut = sprintf('%sFigures\\Cross Correlation Values\\CrossCorr-%s-PerCh-Freq%03i-%03iHz', dataDir, subjectID, blpCorrelations(freqBand).bandPassFrequencies(1), blpCorrelations(freqBand).bandPassFrequencies(2));
    print(figPerChannel, [fileOut '.png'],  '-dpng', '-r300');
    savefig(figPerChannel, [fileOut '.fig']);
    close(figPerChannel);
end

