
ptNumber = 139;
subjectID = num2str(ptNumber);
dataDir = ['E:\Data\ECoG Task-Rest\' subjectID '\'];
restDataFile = [dataDir 'Rest\RestingCorrelations.mat'];
taskDataFile = [dataDir 'Task\TaskActivations.mat'];

load(restDataFile);
load(taskDataFile);

% numChannels = size(blpCorrelations(freqBand).restingCorrelationMatrix, 1);
numChannels = size(channels, 2);
maxCorrelation = zeros(numChannels, 1);
maxActivation = zeros(numChannels, 1);

timeStart = preMoveLength;
timeStop = postMoveLength + 1.0 * samplingRate;

for freqBand = 1:size(blpCorrelations, 2)
    % find the closest frequency in the freqList to the blpCorrelation pair
    [~, freqIdx1] = min(abs(freqList - blpCorrelations(freqBand).bandPassFrequencies(1) ));
    [~, freqIdx2] = min(abs(freqList - blpCorrelations(freqBand).bandPassFrequencies(2) ));
    
    % find the max coreelation and activation scores for all channels
    for channelIdx = 1:numChannels
        channel = channels(channelIdx);
        crossChannelIdxList = 1:numChannels;
        crossChannelIdxList(crossChannelIdxList == channelIdx) = [];      % remove channel from channel list to discard auto-correlation value
        maxCorrelation(channelIdx) = max(blpCorrelations(freqBand).restingCorrelationMatrix(channelIdx, crossChannelIdxList));
        maxActivation(channelIdx) = max( mean( mean( diffSpectra(freqIdx1:freqIdx2, channel, timeStart:timeStop), 3), 1) );
        if(max( -mean( mean( diffSpectra(freqIdx1:freqIdx2, channel, timeStart:timeStop), 3), 1) ))
            maxActivation(channelIdx) = -max( -mean( mean( diffSpectra(freqIdx1:freqIdx2, channel, timeStart:timeStop), 3), 1) );
        end
    end
    
    % creat figure
    figScatterPlot = figure();
    plot(maxActivation, maxCorrelation, 'o');
    hold on;
    
    % calculate and plot linear regression line
    X = [ones(length(maxActivation),1) maxActivation];
    b = X \ maxCorrelation;
    plot(maxActivation, X*b, '-k')

    % calculate correlation coefficient
    corrCoef = corrcoef(maxActivation, maxCorrelation);
    Rsquare = corrCoef(1, 2) ^ 2;
    
    % clean up figure
    title(['Movement - timeFreqBand ' num2str(blpCorrelations(freqBand).bandPassFrequencies(1)) ' - ' num2str(blpCorrelations(freqBand).bandPassFrequencies(2)) ' - R2=' num2str(Rsquare)]);
    xlabel('Max Task Activation');
    ylabel('Max Rest Correlation');
    ylim(gca, [0 1.0]);
    
    % save PNG figure
    fileOut = sprintf('%sFigures\\Scatter Max Task v Rest - %s - FreqBand %03i-%03iHz.png', dataDir, subjectID, blpCorrelations(freqBand).bandPassFrequencies(1), blpCorrelations(freqBand).bandPassFrequencies(2));
    print(figScatterPlot, fileOut,  '-dpng');    % save figure
    close(figScatterPlot);
end