% function RestingCorrelations(ptNumber)

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
outDataFile = [dataDir 'Rest\RestingCorrelations.mat'];


numChannels = 64;               % MAGICNUMBER: only use first 64 channels in each task file
channelToPlotList = [8 24];     % MAGICNUMBER: electrodes 8 & 24 are in motor for Pt 139


% load meta data
metaData = load(metaDataFile, 'DataStructure');
metaData = metaData.DataStructure;
numChannels = min(numChannels, max(metaData.Channels));
channels = metaData.Channels(metaData.Channels <= numChannels);     % limit to numChannels
carGroups = metaData.CARgroups;

% define variables
saveData = struct();
samplingRate = metaData.samplefreq;
numCARGroups = size(carGroups, 2);
freqList = 2.^(1:.125:8);                   % exponential frequency distribution
% freqList = 2.^(1:.25:8);                    % exponential frequency distribution

saveData.subjectID = subjectID;
saveData.restingDataFile = restingDataFile;


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

%% Process Resting Signal
restingSamplingRate = 512;      %MAGICNUMBER: 512Hz sampling rate for clinical ECoG system
restingSignal = load(restingDataFile, 'signals');
restingSignal = double(restingSignal.signals);
signalLen = size(restingSignal, 1);

% normalize signal by CAR groups
meanRestingSignalCARGroups = cell(numCARGroups, 1);
restingSignalReRef = zeros(size(restingSignal));
for carGroup = 1:numCARGroups      % average and re-reference the signal by group
    meanRestingSignalCARGroups{carGroup} = mean(restingSignal(:, carGroupsSansNoisy{carGroup}), 2);     % don't include noisy electrodes in the mean
%     restingSignalReRef(:, carGroups{carGroup}) = restingSignal(:, carGroups{carGroup}) - repmat( meanRestingSignalCARGroups{carGroup}, size(carGroups{carGroup}) );
    restingSignalReRef(:, carGroups{carGroup}) = bsxfun(@minus, restingSignal(:, carGroups{carGroup}), meanRestingSignalCARGroups{carGroup});
end

% arctan spike filter (filter any spikes that jump above 5 standard deviations)
restingSignal5Std = repmat(5 * std(restingSignalReRef), size(restingSignalReRef, 1), 1);
restingSignalReRef = restingSignal5Std .* atan(restingSignalReRef ./ restingSignal5Std);

% notch filter 60 & 120 Hz mains hum
mainsNotchFilter60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', restingSamplingRate);
mainsNotchFilter120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 119, 'HalfPowerFrequency2', 121, 'DesignMethod', 'butter', 'SampleRate', restingSamplingRate);
restingSignalReRefNotched = filtfilt(mainsNotchFilter60, restingSignalReRef);
restingSignalReRefNotched = filtfilt(mainsNotchFilter120, restingSignalReRefNotched);

% Compute BLP (band-limited power), i.e. slow rhythms of power envelope
saveData.blpCorrelations = struct.empty;
bandPassFrequencyList = [8, 10; 8, 13; 13, 30; 35, 50; 70, 100; 100, 140];
for bandPassFrequencies = bandPassFrequencyList'
    bandPassFilter = designfilt('bandpassiir', 'FilterOrder', 4, 'HalfPowerFrequency1', bandPassFrequencies(1), 'HalfPowerFrequency2', bandPassFrequencies(2), 'SampleRate', restingSamplingRate);
    restingSignalReRefNotchedBandPass = filtfilt(bandPassFilter, restingSignalReRefNotched);

    % compute envelope of low pass filtered signal
    % bandPassEnvelope = restingSignalReRefNotchedBandPass .^ 2;
    restingBandPassEnvelope = abs(hilbert(restingSignalReRefNotchedBandPass));

    % low pass filter the envelope
    BLPFrequencies = [0.1 1];
    % lowPassFilter = designfilt('lowpassiir', 'FilterOrder', 4, 'PassbandFrequency', BLPFrequencies(2), 'SampleRate', samplingRate);
    lowPassFilter = designfilt('bandpassiir', 'FilterOrder', 4, 'HalfPowerFrequency1', BLPFrequencies(1), 'HalfPowerFrequency2', BLPFrequencies(2), 'SampleRate', restingSamplingRate);
    restingBandPassEnvelopeFiltered = filtfilt(lowPassFilter, restingBandPassEnvelope);

    % compute correlation matrix 
    restingCorrelationMatrix = corrcoef(restingBandPassEnvelopeFiltered(:, channels));     % limit only to useful channels (electrodes)

    % save the data to disc
    saveData.blpCorrelations(end+1).bandPassFrequencies = bandPassFrequencies;
    saveData.blpCorrelations(end).restingCorrelationMatrix = restingCorrelationMatrix;

    % free up some memory
    % clear('restingSignal', 'restingSignalReRef', 'restingSignalReRefNotched', 'restingSignalReRefNotchedBandPass', 'restingBandPassEnvelope');

    % plot the correlation matrix
    channelToPlotList = [8 24];
    for channelToPlot = channelToPlotList
        figCorrelationMatrix = figure();
        colorMapWhiteMiddle = jet();
        colorMapWhiteMiddle(30:35, :) = repmat([1, 1, 1], 6, 1);    % white out the middle 6 indices

%         rValid = restingCorrelationMatrix(channelToPlot, :);        % use max relative to individual channel
        rValid = restingCorrelationMatrix(channelToPlotList, :);    % use max relative to only channels being plotted
%         rValid = restingCorrelationMatrix;                          % use max of all channels
        rValid = rValid(rValid~=1);                             % exclude auto-correlation from max correlation
        cLimMax = max(abs(rValid));
        cLimMin = -cLimMax;
        % cLimMin = min(abs(restingCorrelationMatrix(channelToPlot, :)));
        % cLimMax = max(abs(restingCorrelationMatrix(channelToPlot, :)));
        imagesc(flipud(rot90(reshape([squeeze(restingCorrelationMatrix(channelToPlot, :)), 0, 0], 8, 8))), [cLimMin cLimMax]) % flip and rotate to get orientation of 1-8 left-to-righ on top row
        colormap(colorMapWhiteMiddle);
        colorbar()
        fileOut = sprintf('%sFigures\\Correlation Matrix - %s - FreqBand %03i-%03iHz - Seed Ch %02i.png', dataDir, subjectID, bandPassFrequencies(1), bandPassFrequencies(2), channelToPlot);
        print(figCorrelationMatrix, fileOut,  '-dpng');    % save figure
        close(figCorrelationMatrix);

        % % plot a sample channel 
        % figure();
        % hold on;
        % plot(restingSignalReRefNotchedBandPass((40000:60000), channelToPlot), 'r');
        % plot(restingBandPassEnvelope((40000:60000), channelToPlot), 'b');
        % plot(restingBandPassEnvelopeFiltered((40000:60000), channelToPlot), 'k');
        % hold off;
    end
end

save(outDataFile, '-struct', 'saveData');

