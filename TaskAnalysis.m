% function TaskRestAnalysis(ptNumber)

%% Setup
% define data files
% % Pt 139 (PT55) data files
ptNumber = 139;
% dataDir = 'E:\Data\ECoG Task-Rest\139\';
% metaDataFile = [dataDir 'Task\PT55_ReachingTask_DataStructure.mat'];
% kineticsDataFile = [dataDir 'Task\PT55_ReachingTask3D_Contra_Kinematics_All.mat'];
% trialsDataFile = [dataDir 'Task\PT55_ReachingTask3D_Contra_QuestionableTrials_Move_and_Spectra_All.mat'];
% restingDataFile = [dataDir 'Rest\121004-B5-7-INV_7E_EDF_seg2.mat'];

% % Pt 145 (PT__) data files
% dataDir = 'E:\Data\ECoG Task-Rest\145\';
% metaDataFile = [dataDir 'Task\145_ReachingTask_DataStructure.mat'];
% kineticsDataFile = [dataDir 'Task\145_ReachingTask3D_Contra_Kinematics_All.mat'];
% trialsDataFile = [dataDir 'Task\145_ReachingTask3D_Contra_QuestionableTrials_Move_and_Spectra_All.mat'];
% restingDataFile = [dataDir 'Rest\.mat'];

subjectID = num2str(ptNumber);
dataDir = ['E:\Data\ECoG Task-Rest\' subjectID '\'];
metaDataFile = [dataDir 'Task\' subjectID '_ReachingTask_DataStructure.mat'];
kineticsDataFile = [dataDir 'Task\' subjectID '_ReachingTask3D_Contra_Kinematics_All.mat'];
trialsDataFile = [dataDir 'Task\' subjectID '_ReachingTask3D_Contra_QuestionableTrials_Move_and_Spectra_All.mat'];
% restingDataFile = [dataDir 'Rest\ .mat'];
restingDataFile = [dataDir 'Rest\121004-B5-7-INV_7E_EDF_seg2.mat'];
outDataFile = [dataDir 'Task\TaskActivations.mat'];


numChannels = 64;	% MAGICNUMBER: only use first 64 channels in each task file


% load meta data
metaData = load(metaDataFile, 'DataStructure');
metaData = metaData.DataStructure;
numChannels = min(numChannels, max(metaData.Channels));
channels = metaData.Channels(metaData.Channels <= numChannels);     % limit to numChannels
carGroups = metaData.CARgroups;

% define variables
flags.saveFigures = false; 
saveData = struct();
samplingRate = metaData.samplefreq;
numCARGroups = size(carGroups, 2);
freqList = 2.^(1:.125:8);                   % exponential frequency distribution
% freqList = 2.^(1:.25:8);                    % exponential frequency distribution

% % Nick's frequency and span lists
% freqList = [1, 1.3748, 1.8593, 2.4754, 3.2467, 4.198, 5.3545, 6.7417, 8.3839, 10.3043, 12.5237, 15.0601, 17.9281, 21.1386, 24.6979, 28.6084, 32.8676, 37.4689, 42.4012, 47.6496, 71.3898, 77.8856, 84.5496, 91.3519, 98.2627, 105.2527, 133.4449, 140.4231, 147.3279, 154.1397, 160.8411, 167.4164, 192.2115, 197.9879, 203.5827, 208.9922];
% freqSpan = [1, 1.0635, 1.131, 1.2028, 1.2791, 1.3603, 1.4467, 1.5385, 1.6362, 1.7401, 1.8505, 1.968, 2.0929, 2.2258, 2.3671, 2.5173, 2.6771, 2.8471, 3.0278, 3.22, 4.1188, 4.3803, 4.6584, 4.9541, 5.2686, 5.603, 7.167, 7.622, 8.1058, 8.6204, 9.1676, 9.7496, 12.4711, 13.2627, 14.1046, 15];

% index of frequencies to use in Time-Power plots       %TODO replace with a lookup fxn to the nearest frequency bin
% frequencyIndexToPlot = [10 24];       % frequency index at 10 & 100Hz for Nick's list
frequencyIndexToPlot = [10 46];       % frequency index at 10 & 100Hz for Carl's list (used w/ Zac's Gabor)


%% Common average re-referencing by channel groups
% group our channels (by amplifier) to common-average reference within each group
carGroupsSansNoisy = cell(numCARGroups, 1);
for carGroup = 1:numCARGroups      % average and re-reference the signal by group
    carGroups{carGroup}(carGroups{carGroup} > numChannels) = [];                        % limit to numChannels
    noisyBooleanIndices = ismember(carGroups{carGroup}, metaData.NoisyChannels);        % index of noisy channels
    carGroupsSansNoisy{carGroup} = carGroups{carGroup};                               	% make a copy to edit
    carGroupsSansNoisy{carGroup}(noisyBooleanIndices) = [];                             % remove noisy channels from CAR groups
end

% %% Process Resting Signal
% restingSamplingRate = 512;      %MAGICNUMBER: 512Hz sampling rate for clinical ECoG system
% restingSignal = load(restingDataFile, 'signals');
% restingSignal = double(restingSignal.signals);
% signalLen = size(restingSignal, 1);
% 
% % normalize signal by CAR groups
% meanRestingSignalCARGroups = cell(numCARGroups, 1);
% restingSignalReRef = zeros(size(restingSignal));
% for carGroup = 1:numCARGroups      % average and re-reference the signal by group
%     meanRestingSignalCARGroups{carGroup} = mean(restingSignal(:, carGroupsSansNoisy{carGroup}), 2);     % don't include noisy electrodes in the mean
% %     restingSignalReRef(:, carGroups{carGroup}) = restingSignal(:, carGroups{carGroup}) - repmat( meanRestingSignalCARGroups{carGroup}, size(carGroups{carGroup}) );
%     restingSignalReRef(:, carGroups{carGroup}) = bsxfun(@minus, restingSignal(:, carGroups{carGroup}), meanRestingSignalCARGroups{carGroup});
% end
% 
% % arctan spike filter (filter any spikes that jump above 5 standard deviations)
% restingSignal5Std = repmat(5 * std(restingSignalReRef), size(restingSignalReRef, 1), 1);
% restingSignalReRef = restingSignal5Std .* atan(restingSignalReRef ./ restingSignal5Std);
% 
% % notch filter 60 & 120 Hz mains hum
% mainsNotchFilter60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', restingSamplingRate);
% mainsNotchFilter120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 119, 'HalfPowerFrequency2', 121, 'DesignMethod', 'butter', 'SampleRate', restingSamplingRate);
% restingSignalReRefNotched = filtfilt(mainsNotchFilter60, restingSignalReRef);
% restingSignalReRefNotched = filtfilt(mainsNotchFilter120, restingSignalReRefNotched);
% 
% % Compute BLP (band-limited power), i.e. slow rhythms of power envelope
% % band-pass filter
% bandPassFrequencyList = [8, 10; 35, 50; 70, 100];
% for bandPassFrequencies = bandPassFrequencyList'
%     bandPassFilter = designfilt('bandpassiir', 'FilterOrder', 4, 'HalfPowerFrequency1', bandPassFrequencies(1), 'HalfPowerFrequency2', bandPassFrequencies(2), 'SampleRate', restingSamplingRate);
%     restingSignalReRefNotchedBandPass = filtfilt(bandPassFilter, restingSignalReRefNotched);
% 
%     % compute envelope of low pass filtered signal
%     % bandPassEnvelope = restingSignalReRefNotchedBandPass .^ 2;
%     restingBandPassEnvelope = abs(hilbert(restingSignalReRefNotchedBandPass));
% 
%     % low pass filter the envelope
%     BLPFrequencies = [0.1 1];
%     % lowPassFilter = designfilt('lowpassiir', 'FilterOrder', 4, 'PassbandFrequency', BLPFrequencies(2), 'SampleRate', samplingRate);
%     lowPassFilter = designfilt('bandpassiir', 'FilterOrder', 4, 'HalfPowerFrequency1', BLPFrequencies(1), 'HalfPowerFrequency2', BLPFrequencies(2), 'SampleRate', restingSamplingRate);
%     restingBandPassEnvelopeFiltered = filtfilt(lowPassFilter, restingBandPassEnvelope);
% 
%     % compute correlation matrix 
%     restingCorrelationMatrix = corrcoef(restingBandPassEnvelope(:, channels));     % limit only to useful channels (electrodes)
% 
%     % free up some memory
%     % clear('restingSignal', 'restingSignalReRef', 'restingSignalReRefNotched', 'restingSignalReRefNotchedBandPass', 'restingBandPassEnvelope');
% 
%     % plot the correlation matrix
%     channelToPlotList = [8 24];
%     for channelToPlot = channelToPlotList
%         figCorrelationMatrix = figure();
%         colorMapWhiteMiddle = jet();
%         colorMapWhiteMiddle(30:35, :) = repmat([1, 1, 1], 6, 1);    % white out the middle 6 indices
% 
% %         rValid = restingCorrelationMatrix(channelToPlot, :);        % use max relative to individual channel
%         rValid = restingCorrelationMatrix(channelToPlotList, :);    % use max relative to only channels being plotted
% %         rValid = restingCorrelationMatrix;                          % use max of all channels
%         rValid = rValid(rValid~=1);                             % exclude auto-correlation from max correlation
%         cLimMin = -max(abs(rValid));
%         cLimMax = max(abs(rValid));
%         % cLimMin = min(abs(restingCorrelationMatrix(channelToPlot, :)));
%         % cLimMax = max(abs(restingCorrelationMatrix(channelToPlot, :)));
%         imagesc(flipud(rot90(reshape([squeeze(restingCorrelationMatrix(channelToPlot, :)), 0, 0], 8, 8))), [cLimMin cLimMax]) % flip and rotate to get orientation of 1-8 left-to-righ on top row
%         colormap(colorMapWhiteMiddle);
%         colorbar()
%         fileOut = sprintf('%sFigures\\Correlation Matrix - %s - FreqBand %03i-%03iHz - Seed Ch %02i.png', dataDir, subjectID, bandPassFrequencies(1), bandPassFrequencies(2), channelToPlot);
%         print(figCorrelationMatrix, fileOut,  '-dpng');    % save figure
%         close(figCorrelationMatrix);
% 
%         % % plot a sample channel 
%         % figure();
%         % hold on;
%         % plot(restingSignalReRefNotchedBandPass((40000:60000), channelToPlot), 'r');
%         % plot(restingBandPassEnvelope((40000:60000), channelToPlot), 'b');
%         % plot(restingBandPassEnvelopeFiltered((40000:60000), channelToPlot), 'k');
%         % hold off;
%     end
% end
% 
% return

%% Process Task Data
% loop through data files and get terminal folder and filename for path
taskDataFiles = cell([size(metaData.ContraFile, 2) 1]);
for i = 1:size(metaData.ContraFile, 2)
    [path, name, ~] = fileparts(metaData.ContraFile(i).name);
    [~, string] = regexp(path, '\\(\w+)$', 'match', 'lineanchors'); % get terminal folder name
    taskDataFiles{i} = [dataDir 'Task\' path(string+1:end) '\' name '.mat'];
end

% define our epoch lengths
% itiLength = 2 * samplingRate;       % ITI lengths are usually 3 seconds; we will use the first _ seconds
itiLength = 0.8 * samplingRate;     % use the HoldA period starting 0.2s after until end of 1s
% moveLength = 1 * samplingRate;    % 80% of trials are < 1.0 sec, see histogram((moveEnd - moveOnset)/samplingRate)
preMoveLength = 3 * samplingRate;
postMoveLength = 4 * samplingRate;
moveLength = preMoveLength + postMoveLength;

% pre-allocate memory
numTaskFiles = size(taskDataFiles, 1);
% taskSignal = cell(numTaskFiles, 1);
meanTaskSignalCARGroups = cell(numTaskFiles, numCARGroups);
taskSignalReRef = cell(numTaskFiles, 1);
taskSignalReRefNotched = cell(numTaskFiles, 1);
taskLengthConcat = 0;
taskFileStart = cell(numTaskFiles, 1);
itiMeanSpectra = zeros(numTaskFiles, length(freqList), numChannels, itiLength);     % TODO: numChannels hard-coded b/c we won't know until we're inside the loop
trialMeanSpectra = zeros(numTaskFiles, length(freqList), numChannels, moveLength);

% load kinetics and trial label data
kineticsData = load(kineticsDataFile, 'HoldAStarts_used', 'MoveOnsets_used', 'HoldBStarts_used');   % 'ITIStarts_used', 
moveOnset = kineticsData.MoveOnsets_used;
% moveEnd = kineticsData.HoldBStarts_used;    % defines move end by start of Hold B section (which is defined as when subject actually reached the target)
% itiStarts = kineticsData.ITIStarts_used;
itiStarts = kineticsData.HoldAStarts_used + 0.2*samplingRate;       % use 0.2s after start of HoldA period for ITI 
trialsData = load(trialsDataFile, 'GoodTrials');
goodTrials = trialsData.GoodTrials;
clear('kineticsData', 'trialsData');

% loop through all task files and process signal individually
for i = 1:numTaskFiles
    fprintf('Processing file #%d of %d...\n', i, numTaskFiles);
    taskData = load(taskDataFiles{i}, 'signal', 'states', 'params');
    taskSignal = taskData.signal(:, 1:numChannels);
    taskLength = size(taskSignal, 1);
    taskFileStart{i} = taskLengthConcat + 1;
    taskLengthConcat = taskLengthConcat + taskLength;
    
    % common average re-references by groups
    taskSignalReRef{i} = zeros(size(taskSignal));
    for carGroup = 1:numCARGroups      % average and re-reference the signal by group
        meanTaskSignalCARGroups{i, carGroup} = mean(taskSignal(:, carGroupsSansNoisy{carGroup}), 2);
%         taskSignalReRef{i}(:, carGroups{carGroup}) = taskSignal{i}(:, carGroups{carGroup}) - repmat( meanTaskSignalCARGroups{i, carGroup}, size(carGroups{carGroup}) );
        taskSignalReRef{i}(:, carGroups{carGroup}) = bsxfun(@minus, taskSignal(:, carGroups{carGroup}), meanTaskSignalCARGroups{i, carGroup});
    end
    
    % arctan spike filter (filter any spikes that jump above 5 standard deviations)
    taskSignal5Std = repmat(5 * std(taskSignalReRef{i}), size(taskSignalReRef{i}, 1), 1);
    taskSignalReRef{i} = taskSignal5Std .* atan(taskSignalReRef{i} ./ taskSignal5Std);

    % notch filter 60 & 120 Hz mains hum
    mainsNotchFilter60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', samplingRate);
    mainsNotchFilter120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 119, 'HalfPowerFrequency2', 121, 'DesignMethod', 'butter', 'SampleRate', samplingRate);
    taskSignalReRefNotched{i} = filtfilt(mainsNotchFilter60, taskSignalReRef{i});
    taskSignalReRefNotched{i} = filtfilt(mainsNotchFilter120, taskSignalReRefNotched{i});
    
    % normalize signal by subtracting mean and dividing inter-quartile range
    taskSignalReRefNotchedMean = median(taskSignalReRefNotched{i}, 1);
%     taskSignalReRefNotchedIQR = iqr(taskSignalReRefNotched{i}, 1);
%     taskSignalNorm = (taskSignalReRefNotched{i} - repmat(taskSignalReRefNotchedMean, taskLength, 1)) ./ repmat(taskSignalReRefNotchedIQR, taskLength, 1);   % normalize by interquartile range
    taskSignalReRefNotchedStd = std(taskSignalReRefNotched{i}, 0, 1);
%     taskSignalNorm = (taskSignalReRefNotched{i} - repmat(taskSignalReRefNotchedMean, taskLength, 1)) ./ repmat(taskSignalReRefNotchedStd, taskLength, 1);   % normalize by standard deviation (z-score)
    taskSignalNorm = bsxfun(@rdivide, bsxfun(@minus, taskSignalReRefNotched{i},  taskSignalReRefNotchedMean), taskSignalReRefNotchedStd);
    
    % calculate spectra for the whole file signal
    taskSpectra = gabor_cov_fitted(taskSignalNorm, freqList, samplingRate, 3);                                  % Zac's Gabor
%     taskSpectra = permute(abs(gabor_response(taskSignalNorm, freqList, freqSpan, samplingRate)), [3 2 1]);      % Nick's Gabor (CPU)
%     taspSpectraGPU = permute(abs(gpu_gabor_response(taskSignalNorm, freqList, freqSpan, samplingRate)), [3 2 1]);
%     s = spectrogram(taskSignalNorm(:, 1), 10, [], freqList, samplingRate, 'PSD');

    % find the index to the first move trial
    firstTrialIndex = find(moveOnset >= taskFileStart{i}, 1);
    
    % segment kinetic data by task file
    taskITIStarts = itiStarts(itiStarts >= taskFileStart{i} & itiStarts <= taskLengthConcat);
    taskMoveOnset = moveOnset(moveOnset >= taskFileStart{i} & moveOnset <= taskLengthConcat);
    taskGoodTrials = goodTrials(moveOnset(goodTrials) >= taskFileStart{i} & moveOnset(goodTrials) <= taskLengthConcat);
    
    % references ITI and Move onsets and trial index relative to start of task (ie. correct for concat)
    taskITIStarts = taskITIStarts - taskFileStart{i} + 1;
    taskMoveOnset = taskMoveOnset - taskFileStart{i} + 1;
    taskGoodTrials = taskGoodTrials - firstTrialIndex + 1;
    
    % assert: number of good trials does not exceed total trials
    numTrials = size(taskGoodTrials, 2);
    if(numTrials > size(taskMoveOnset, 2))
        error('Assertion error: Number of good trials is greater than number of total trials for this task file (run).');
    end

    % segment signal into trial epochs by ITI and movement onset
    itiSpectra = zeros(numTrials, size(taskSpectra, 1), size(taskSpectra, 2), itiLength);
    trialSpectra = zeros(numTrials, size(taskSpectra, 1), size(taskSpectra, 2), moveLength);
    for trial = 1:numTrials-1
        itiSpectra(trial, :, :, :)   = taskSpectra(:, :, taskITIStarts(taskGoodTrials(trial)):taskITIStarts(taskGoodTrials(trial))+itiLength-1);
        trialSpectra(trial, :, :, :) = taskSpectra(:, :, taskMoveOnset(taskGoodTrials(trial))-preMoveLength:taskMoveOnset(taskGoodTrials(trial))+postMoveLength-1);
    end

    % average spectra across trials
    itiMeanSpectra(i, :, :, :) = mean(itiSpectra, 1);
    trialMeanSpectra(i, :, :, :) = mean(trialSpectra, 1);

end % for i = 1:numTaskFiles

% free up some memory
clear('taskData', 'taskSignal', 'taskSignalReRef', 'meanTaskSignalCARGroups');

% average spectra across runs (ie. each task files)
itiMeanSpectra = squeeze(mean(itiMeanSpectra, 1));
trialMeanSpectra = squeeze(mean(trialMeanSpectra, 1));

% % isolate spectra for a particular trial (ie. each task file)
% itiMeanSpectra = squeeze(itiMeanSpectra(1, :, :, :));
% trialMeanSpectra = squeeze(trialMeanSpectra(1, :, :, :));


%% Plot spectrogram
% high pass the difference spectra
itiTimeMeanSpectra = mean(itiMeanSpectra, 3);
% itiTimeMeanSpectra = repmat(itiTimeMeanSpectra, 1, 1, size(trialMeanSpectra, 3));   % time averaged ITI spectra
% diffSpectra = (trialMeanSpectra - itiTimeMeanSpectra) ./ itiTimeMeanSpectra;        % subtract and divide the time averaged ITI spectra
diffSpectra = bsxfun(@rdivide, bsxfun(@minus, trialMeanSpectra, itiTimeMeanSpectra), itiTimeMeanSpectra);   % subtract and divide the time averaged ITI spectra
cLimMax = max(max(max(abs(diffSpectra))));

% create color map with a grey mid region
cMapLen = 65; 
cMapGreyWidth = 4; 
cMapGreyIntensity = 0.95; % 1.0 for white
greyZeroJet = jet(cMapLen);
greyZeroJet(ceil(cMapLen/2)-cMapGreyWidth:floor(cMapLen/2)+cMapGreyWidth, :) = repmat([cMapGreyIntensity, cMapGreyIntensity, cMapGreyIntensity], cMapGreyWidth*2+mod(cMapLen+1,2), 1);

figSpectrogram = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for channel = channels
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
    if channel >= 56
        axis.XTick = 1:round(samplingRate):size(trialMeanSpectra, 3);
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
    fileOut = [dataDir 'Figures\Spectrogram - ' subjectID ' - Move vs HoldA - Zac''s Gabor.fig'];
    savefig(figSpectrogram, fileOut);
    print(figSpectrogram, fileOut, '-dpng');
    clear('fileOut');
    close(figSpectrogram);
end

%% Plot Time-Power figures for a specific frequencies
figTimePowerPlot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for channel = channels
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
    fileOut = [dataDir 'Figures\TimePowerPlot - ' subjectID ' - ' num2str(freqList(frequencyIndexToPlot(1))) 'Hz & ' num2str(freqList(frequencyIndexToPlot(2))) 'Hz - Move vs HoldA - Zac''s Gabor.fig'];
    savefig(figTimePowerPlot, fileOut);
    print(figTimePowerPlot, fileOut, '-dpng');s
    clear('fileOut');
    close(figTimePowerPlot);
end

% end % function TaskRestAnalsysi()