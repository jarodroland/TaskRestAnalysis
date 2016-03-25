% function MonkeyTaskRestAnalysis(ptNumber)

%% Setup
flags.saveFigures = false; 
flags.plotNoisyChannel = false;

% define data files
subjectID = 'Monkey K';
dataDir = ['E:\Data\Moran Lab ECoG Data\' subjectID '\'];
metaDataFile = [dataDir 'Events_K20120322.mat'];
taskDataFile = [dataDir 'raw_ecog_64chans_K20120322.mat'];
restingDataFile = [dataDir 'baseline_10days.mat'];


numChannels = 64;	% MAGICNUMBER: only use first 64 channels in each task file


% load monkey data if not already loaded
if(exist('signal', 'var') ~= 1)
    metaData = load(metaDataFile, 'event', 'start_time', 'stop_time');
    taskData = load(taskDataFile, 'chans', 'time', 'Fs', 'bad_channels');
    signal = taskData.chans;
    samplingRate = taskData.Fs;
    noisyChannels = taskData.bad_channels;
end

%DEBUG
signal = signal(:, 1:numChannels);  % truncture signal to numChannels size

% define variables
signalLen = size(signal, 1);
channels = 1:min(numChannels, size(signal, 2));
numChannels = length(channels);
numEvents = size(metaData.event, 2);

% create a simulus code sequence from the event codes and time stamps
stimulusCode = zeros(signalLen, 1);
for event = 1:numEvents
    stimulusCode(round( metaData.start_time(event) * samplingRate ) : min(round( metaData.stop_time(event) * samplingRate ), signalLen)) = metaData.event(event);
end


% clear('metaData');

% freqList = 2.^(1:.25:8);                    % exponential frequency distribution
freqList = 2.^(1:.125:8);                   % exponential frequency distribution

% % Nick's frequency and span lists
% freqList = [1, 1.3748, 1.8593, 2.4754, 3.2467, 4.198, 5.3545, 6.7417, 8.3839, 10.3043, 12.5237, 15.0601, 17.9281, 21.1386, 24.6979, 28.6084, 32.8676, 37.4689, 42.4012, 47.6496, 71.3898, 77.8856, 84.5496, 91.3519, 98.2627, 105.2527, 133.4449, 140.4231, 147.3279, 154.1397, 160.8411, 167.4164, 192.2115, 197.9879, 203.5827, 208.9922];
% freqSpan = [1, 1.0635, 1.131, 1.2028, 1.2791, 1.3603, 1.4467, 1.5385, 1.6362, 1.7401, 1.8505, 1.968, 2.0929, 2.2258, 2.3671, 2.5173, 2.6771, 2.8471, 3.0278, 3.22, 4.1188, 4.3803, 4.6584, 4.9541, 5.2686, 5.603, 7.167, 7.622, 8.1058, 8.6204, 9.1676, 9.7496, 12.4711, 13.2627, 14.1046, 15];
% frequencyIndexToPlot = [10 24];       % index of frequencies (10 & 100Hz) to use in Time-Power plots, used for Nick's list



%% Common average re-referencing 
% % group our channels (by amplifier) to common-average reference within each group
carGroups = channels;
carGroups(carGroups > numChannels) = [];                        % limit to numChannels
noisyBooleanIndices = ismember(carGroups, noisyChannels);       % index of noisy channels
carGroupsSansNoisy = carGroups;                       % make a copy to edit
carGroupsSansNoisy(noisyBooleanIndices) = [];                   % remove noisy channels from CAR groups

% normalize signal (Z score)
signalReRef = zeros(size(signal));
meanSignalCARGroups = mean(signal(:, carGroupsSansNoisy), 2);     % don't include noisy electrodes in the mean
signalReRef(:) = bsxfun(@minus, signal, meanSignalCARGroups);

%% Pre-Process Signal

% arctan spike filter (filter any spikes that jump above 5 standard deviations)
signal5Std = repmat(5 * std(signalReRef), size(signalReRef, 1), 1);
signalReRef = signal5Std .* atan(signalReRef ./ signal5Std);

% notch filter 60 & 120 Hz mains hum
mainsNotchFilter60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', samplingRate);
mainsNotchFilter120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 119, 'HalfPowerFrequency2', 121, 'DesignMethod', 'butter', 'SampleRate', samplingRate);
signalReRefNotched = filtfilt(mainsNotchFilter60, signalReRef);
signalReRefNotched = filtfilt(mainsNotchFilter120, signalReRefNotched);

% normalize signal by subtracting mean and dividing standard-deviation or inter-quartile range
signalReRefNotchedMean = median(signalReRefNotched, 1);
% taskSignalReRefNotchedIQR = iqr(taskSignalReRefNotched{i}, 1);
signalReRefNotchedStd = std(signalReRefNotched, 0, 1);
signalNorm = bsxfun(@rdivide, bsxfun(@minus, signalReRefNotched,  signalReRefNotchedMean), signalReRefNotchedStd);

clear('signalReRefNotchedStd', 'signalReRefNotchedMean', 'signalReRefNotched', 'signalReRef', 'signal5Std')

% subChannels = [1:4; 5:8];%; 9:12; 13:16; 17:20; 21:24; 25:28; 29:32; 33:36; 37:40; 41:43; 44:48; 49:52; 53:56; 57:60; 61:64];
subChanSize = 4;    % number of channels to process at a time %NOTE: must be divsor of numChannels
subChannels = reshape(1:numChannels, subChanSize, numChannels/subChanSize)';
numSubChan = size(subChannels, 1);

%% Process Task Data

% define our epochs
holdIndex = find(metaData.event == 2);      %MAGICNUMBER: 2 defines center hold period
moveIndex = find(metaData.event == 3);      %MAGICNUMBER: 3 defines movement period
holdStartIndex = round(metaData.start_time(holdIndex) * samplingRate);
holdStopIndex = round(metaData.stop_time(holdIndex) * samplingRate);
moveStartIndex = round(metaData.start_time(moveIndex) * samplingRate);
moveStopIndex = round(metaData.stop_time(moveIndex) * samplingRate);

itiLen = round(0.5 * samplingRate);        % define ITI as 0.5s before movement start (during hold period)
trialLen = round(2.0 * samplingRate);      % define trial as 2.0s after movement start

itiMeanSpectra = zeros(length(freqList), numChannels, itiLen);
trialMeanSpectra = zeros(length(freqList), numChannels, trialLen);

% process spectra per channel due to large signal size
for subChanListIdx = 1:numSubChan
    subChanList = subChannels(subChanListIdx, :);
    fprintf('Processing channels %n to %n of %n\n', subChanList(1), subChanList(end), numChannels);
    
    % calculate spectra for the whole file signal
    clear('taskSpectra');
    taskSpectra = gabor_cov_fitted(signalNorm(:, subChanList), freqList, samplingRate, 3);                             % Zac's Gabor
    %     taskSpectra = permute(abs(gabor_response(taskSignalNorm, freqList, freqSpan, samplingRate)), [3 2 1]);      % Nick's Gabor (CPU)
    %     taspSpectraGPU = permute(abs(gpu_gabor_response(taskSignalNorm, freqList, freqSpan, samplingRate)), [3 2 1]);
    %     s = spectrogram(taskSignalNorm(:, 1), 10, [], freqList, samplingRate, 'PSD');

    % average across all ITI and Trial epochs
    numTrials = length(moveStartIndex);
    for i = 1:numTrials
        itiMeanSpectra(:, subChanList, :)   = itiMeanSpectra(:, subChanList, :)   + taskSpectra(:, :, moveStartIndex(i)-itiLen:moveStartIndex(i)-1);
        trialMeanSpectra(:, subChanList, :) = trialMeanSpectra(:, subChanList, :) + taskSpectra(:, :, moveStartIndex(i)       :moveStartIndex(i)+trialLen-1);
    end
    itiMeanSpectra(:, subChanList, :) = itiMeanSpectra(:, subChanList, :) ./ numTrials;
    trialMeanSpectra(:, subChanList, :) = trialMeanSpectra(:, subChanList, :) ./ numTrials;
end

% % average spectra across movement periods
% numMove = length(moveStartIndex);
% moveStartStopIndex = zeros(1, numMove * 2);
% moveStartStopIndex(1:2:numMove*2) = moveStartIndex;
% moveStartStopIndex(2:2:numMove*2) = moveStopIndex;
% itiMeanSpectra = mean(taskSpectra(:, :, moveStartStopIndex), 3);
% 
% % average spectra across movement periods
% numHold = length(holdStartIndex);
% holdStartStopIndex = zeros(1, numHold * 2);
% holdStartStopIndex(1:2:numHold*2) = holdStartIndex;
% holdStartStopIndex(2:2:numHold*2) = holdStopIndex;
% trialMeanSpectra = mean(taskSpectra(:, :, holdStartStopIndex), 3);


%% Plot spectrogram
preMoveLength = itiLen;

% concatenate the iti and trial spectra
catMeanSpectra = cat(3, itiMeanSpectra, trialMeanSpectra);

itiTimeMeanSpectra = mean(itiMeanSpectra, 3);   % time averaged ITI spectra
diffSpectra = bsxfun(@rdivide, bsxfun(@minus, catMeanSpectra, itiTimeMeanSpectra), itiTimeMeanSpectra);   % subtract and divide the time averaged ITI spectra
cLimMax = max(max(max(abs(diffSpectra))));

% create color map with a grey mid region
cMapLen = 65; 
cMapGreyWidth = 4; 
cMapGreyIntensity = 0.95; % 1.0 for white
greyZeroJet = jet(cMapLen);
greyZeroJet(ceil(cMapLen/2)-cMapGreyWidth:floor(cMapLen/2)+cMapGreyWidth, :) = repmat([cMapGreyIntensity, cMapGreyIntensity, cMapGreyIntensity], cMapGreyWidth*2+mod(cMapLen+1,2), 1);

figSpectrogram = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for channel = channels
    % don't plot noisy channels
    if(~flags.plotNoisyChannel && ~isempty(find(noisyChannels == channel, 1)))
        continue;
    end
    
    subplot(8, 8, channel);
    colormap(greyZeroJet);
    
    imagesc(squeeze(diffSpectra(:, channel, :)));
%     imagesc(squeeze(moveMeanSpectra(:, channel, :)));
    
    % plot closest frequency reference lines
    hold on;
    freqRefs = [10, 100];   % (Hz) list reference frequencies to plot horizontal lines at
    for freq = freqRefs
        [~, freqIndex] = min(abs(freqList - freq));
        plot([0, size(diffSpectra, 3)], [freqIndex freqIndex], 'b:');
    end
        
    axis = gca;
    freqLabelIndex = round(linspace(1, length(freqList), 5));
    set(axis, 'YDir', 'normal')
    caxis([-cLimMax cLimMax]);
    
    % label axis on left and bottom most only
    if mod(channel, 8) == 1
        axis.YTick = freqLabelIndex;
        axis.YTickLabel = round(freqList(freqLabelIndex));
    else
        axis.YTick = freqLabelIndex;
        axis.YTickLabel = [];
    end
    if channel > 56
%         axis.XTick = 1:round(samplingRate):size(trialMeanSpectra, 3);
        axis.XTick = [1, preMoveLength, preMoveLength+trialLen/2, preMoveLength+trialLen];
        axis.XTickLabel = round((axis.XTick-preMoveLength)/samplingRate, 1);
    else
        axis.XTick = 1:round(samplingRate):size(trialMeanSpectra, 3);
        axis.XTickLabel = [];
    end
    
    % plot move onset marker
    hold on;
    plot([preMoveLength preMoveLength], [0 length(freqList)], 'b');
    
%     % colorbar on right most only
%     if mod(channel, 8) == 0
%         colorbar()
%     end;
end

if(flags.saveFigures)
    fileOut = [dataDir 'Figures\Spectrogram - ' subjectID ' - Move vs Hold - Zac''s Gabor.fig'];
    savefig(figSpectrogram, fileOut);
    print(figSpectrogram, fileOut, '-dpng');
    clear('fileOut');
    close(figSpectrogram);
end

%% Plot Time-Power figures for a specific frequencies
frequencyIndexToPlot = [20 46];       % index of frequencies (10 & 100Hz) to use in Time-Power plots, used for Carl's list (used w/ Zac's Gabor) %TODO replace with a lookup fxn to the nearest frequency bin

figTimePowerPlot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for channel = channels
    % don't plot noisy channels
    if(~flags.plotNoisyChannel && ~isempty(find(noisyChannels == channel, 1)))
        continue;
    end
    
    subplot(8, 8, channel);
    plot(squeeze(diffSpectra(frequencyIndexToPlot(2), channel, :)), 'r');
    hold on;
    plot(squeeze(diffSpectra(frequencyIndexToPlot(1), channel, :)), 'b');

    axis = gca;
%     % label axis on left and bottom most only
%     if channel >= 56
%         axis.XTick = 1:round(samplingRate):size(trialMeanSpectra, 3);
%         axis.XTickLabel = round((axis.XTick-preMoveLength)/samplingRate, 1);
%     else
%         axis.XTickLabel = [];
%     end
    
    % plot move onset marker
    
    ylim(axis, [-0.5 0.5]);
    hold on;
    plot([preMoveLength preMoveLength], [-0.5 0.5], 'k');
end

if(flags.saveFigures)
    fileOut = [dataDir 'Figures\TimePowerPlot - ' subjectID ' - ' num2str(freqList(frequencyIndexToPlot(1))) 'Hz & ' num2str(freqList(frequencyIndexToPlot(2))) 'Hz - Move vs Hold - Zac''s Gabor.fig'];
    savefig(figTimePowerPlot, fileOut);
    print(figTimePowerPlot, fileOut, '-dpng');s
    clear('fileOut');
    close(figTimePowerPlot);
end

% end % function TaskRestAnalsysi()