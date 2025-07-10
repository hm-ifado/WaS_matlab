function [eegMerged, streamsMerged] = WaS4_merge(subject, dataPath)

folders = dir(fullfile(dataPath, ['*WaS*' sprintf('%03d', subject) '*']));

if isempty(folders)
    warning(['no data found for subject ' sprintf('%03d', subject)]);
    return
end

subjectFolder = folders(1).name;

%% gather all required data

%% load eeg data from bdf file
files = dir(fullfile(dataPath, subjectFolder,  '*.bdf'));
ampTmp = {};

if isempty(files)
    warning('No EEG files found, stopping...')
    return
end
for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    ampTmp{f} = pop_biosig(filepath);
end

% merge files if > 1
amp = ampTmp{1};
cntIdx = find(strcmp({amp.chanlocs.labels}, 'CNT'));
ampCNTOffsets = [0];
for f = 2:length(ampTmp)
    newTime = amp.times(end) + ampTmp{f}.times + (1000/amp.srate);
    amp.times = [amp.times newTime];

    ampTmp{f}.data(cntIdx, :) = ampTmp{f}.data(cntIdx, :) + amp.data(cntIdx,end); % adjust counter value to continue from the last file
    ampCNTOffsets = [ampCNTOffsets amp.data(cntIdx, end)]; % offsets to adjust counter values in LSL data later (last counter value before new recording)

    amp.data = [amp.data, ampTmp{f}.data];
    amp.pnts = amp.pnts + ampTmp{f}.pnts;
end
disp('--- done reading EEG ---')
% amp <- eeg data


%% load faros data from edf file
files = dir(fullfile(dataPath, subjectFolder, '*.EDF'));
farosData = {};

if isempty(files)
    warning('No Faros files found, stopping...')
    return
end

% read files
for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    farosImport = edfread(filepath);
    farosInfo = edfinfo(filepath);
    ecgSampRate = 1/farosInfo.NumSamples(1);

    newECGData = cell2mat(farosImport.ECG);
    newACCData = [cell2mat(farosImport.Accelerometer_X), cell2mat(farosImport.Accelerometer_Y), cell2mat(farosImport.Accelerometer_Z)];
    newECGTime = (0:size(newECGData, 1)-1) * (1/farosInfo.NumSamples(1));
    newACCTime = (0:size(newACCData, 1)-1) * (1/farosInfo.NumSamples(2));

    newFarosData = [newECGTime', newECGData, zeros(length(newECGTime),size(newACCData,2))];
    for a = 1:size(newACCData,2)
        newFarosData(:, 2+a) = interp1(newACCTime, newACCData(:, a), newECGTime, 'makima');
    end

    farosData{f} = newFarosData;
end
disp('--- done reading ECG ---')
% farosData <- faros timestamps, eeg, acc X/Y/Z


%% load xsens data from mat file
% run convertXsens first to generate .mat
xsensData = [];
try
    files = dir(fullfile(dataPath, subjectFolder, '*xsens*mat'));
    filepath = fullfile(files(1).folder, files(1).name);
    disp(['reading ' filepath '...'])
    xsensTmp = load(filepath);
    xsensData = xsensTmp.xsensData;
    disp('--- done reading Xsens ---')
catch
    warning('No xsens.mat file found, continuing without Xsens data');
end
% xsensData <- xsens center of mass, angles etc


%% load lsl data from xdf file
files = dir(fullfile(dataPath, subjectFolder, '*.xdf'));
streamsTmp = {};

if isempty(files)
    warning('No LSL data found, stopping...');
    return
end

for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    streamsTmp{f} = load_xdf(filepath);
end

% if files > 1 look for timestamp reset (labrecorder pc restart)
if length(streamsTmp) > 1
    for f = 2:length(streamsTmp)
        reset = false;

        idx1 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f-1});
        idx2 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f});

        if streamsTmp{f}{idx2}.time_stamps(1) < streamsTmp{f-1}{idx1}.time_stamps(end)
            reset = true;
        end

        % fix timestamps
        if reset
            warning('Detected timestamp reset in LSL data, attempting to fix...');
            % use eeg counter if possible, ecg otherwise
            if length(ampTmp) == 1
                idxLast = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f-1});
                idxCurrent = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f});

                lastEnd = amp.time(amp.data(cntIdx,:) == streamsTmp{f-1}{idxLast}.time_series(end));
                currentFirst = amp.time(amp.data(cntIdx,:) == streamsTmp{f}{idxCurrent}.time_series(1));

                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(end) + (currentFirst - lastEnd)/1000 - streamsTmp{f}{idxCurrent}.time_stamps(1);
            elseif length(farosData) == 1
                idxLast = cellfun(@(x) contains(x.info.name, 'faros_ecg'), streamsTmp{f-1});
                idxCurrent = cellfun(@(x) contains(x.info.name, 'faros_ecg'), streamsTmp{f});

                newLSLTimeLast = streamsTmp{f-1}{idxLast}.time_stamps(1):ecgSampRate:streamsTmp{f-1}{idxLast}.time_stamps(end);
                newLSLECGLast = interp1(streamsTmp{f-1}{idxLast}.time_stamps, double(streamsTmp{f-1}{idxLast}.time_series), newLSLTimeLast, 'pchip');
                [r,lags] = xcorr(farosData{1}(:,2), newLSLECGLast);
                [~, maxIdx] = max(r);
                offsetLast = lags(maxIdx);

                newLSLTimeCurrent = streamsTmp{f}{idxCurrent}.time_stamps(1):ecgSampRate:streamsTmp{f}{idxCurrent}.time_stamps(end);
                newLSLECGCurrent = interp1(streamsTmp{f}{idxCurrent}.time_stamps, double(streamsTmp{f}{idxCurrent}.time_series), newLSLTimeCurrent, 'pchip');
                [r,lags] = xcorr(farosData{1}(:,2), newLSLECGCurrent);
                [~, maxIdx] = max(r);
                offsetCurrent = lags(maxIdx);

                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(1) + (offsetCurrent - offsetLast)*ecgSampRate - streamsTmp{f}{idxCurrent}.time_stamps(1);
            else
                warning("Could not fix LSL timestamps, stopping...");
                return
            end

            % use offset to adjust times for all parts of LSL data
            for nf = f:length(streamsTmp)
                for d = 1:length(streamsTmp{nf})
                    if ~isempty(streamsTmp{nf}{d}.time_stamps)
                        streamsTmp{nf}{d}.time_stamps = streamsTmp{nf}{d}.time_stamps + lslTimeOffset;
                    end
                end
            end
        end
    end

end

% merge files if > 1
streams = streamsTmp{1};
for f = 2:length(streamsTmp)
    nextStream = streamsTmp{f};
    for n = 1:length(nextStream)
        for id = 1:length(streams)
            if strcmp(streams{id}.info.name, nextStream{n}.info.name)
                if strcmp(streams{id}.info.type, "EEG") && nextStream{n}.time_series(1) < streams{id}.time_series(end)
                    % adjust EEG counter values
                    offsetIdx = find(ampCNTOffsets + nextStream{n}.time_series(1) > streams{id}.time_series(end), 1, 'first');
                    nextStream{n}.time_series = nextStream{n}.time_series + ampCNTOffsets(offsetIdx);
                end
                streams{id}.time_stamps = [streams{id}.time_stamps nextStream{n}.time_stamps];
                streams{id}.time_series = [streams{id}.time_series nextStream{n}.time_series];
            end
        end
    end
end
% streamNames = {};
% for id = 1:length(streams)
%     streamNames{id} = streams{id}.info.name;
% end
disp('--- done reading LSL ---')
% streams <- lsl recorder data


%% load gaze data
plDirs = dir([fullfile(dataPath, subjectFolder, '*pupilcloud_outside*')]);
plDirs = {plDirs.name};

if isempty(plDirs)
    warning('No eyetracker folders found, stopping...');
    return
end

plGazeData = [];
plFixData = [];
plEventData = [];
for f = 1:length(plDirs)
    % import gaze and fixation data
    disp(['reading gaze data from ' fullfile(dataPath, subjectFolder, plDirs{f}) '...'])
    try
        newGaze = importGaze(fullfile(dataPath, subjectFolder, plDirs{f}, 'time_aligned_gaze.csv'));
        newFix = importFixations(fullfile(dataPath, subjectFolder, plDirs{f}, 'fixations.csv'));
        newEvent = readtable(fullfile(dataPath, subjectFolder, plDirs{f}, 'events.csv'));
    catch ex
        warning(['missing files in folder:' newline ex.message]);
        return
    end

    plGazeData = vertcat(plGazeData, newGaze);
    plFixData = vertcat(plFixData, newFix);
    plEventData = vertcat(plEventData, newEvent);
end
disp('--- done reading gaze ---')
% plGazeData, plFixData, plEventData <- pupil labs gaze, fixations, events


%% replace lsl streams with offline recordings


%% eeg
streamsMerged = streams;
% eegIdx = ~cellfun(@isempty,regexp(streamNames,'Counter'));
% eegIdx = contains(streamNames, 'Counter', 'IgnoreCase', true);
eegIdx = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsMerged);

% determine beginning and end of recording
% look for counter gaps in LSL stream and determine b & e for each
% segment to preserve timestamps
gapIdx = find(diff(streams{eegIdx}.time_series) > 1);
gapIdx = [0 gapIdx length(streams{eegIdx}.time_series)];
beg = [];
ending = [];
for g = 1:length(gapIdx)-1
    nextBeg = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(gapIdx(g)+1));
    nextEnding = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(gapIdx(g+1)));

    beg = [beg nextBeg];
    ending = [ending nextEnding];
end

% replace counters with eeg data
newEEG = [];
for b = 1:length(beg)
    newEEG = [newEEG amp.data(1:end-1, beg(b):ending(b))];
end
streamsMerged{eegIdx}.time_series = newEEG;
disp('--- replaced EEG data ---');


%% ecg
% farosECGIdx = contains(streamNames, 'faros_ecg', 'IgnoreCase', true);
% farosACCIdx = contains(streamNames, 'faros_acc', 'IgnoreCase', true);
farosECGIdx = cellfun(@(x) contains(x.info.name, 'faros_ecg', 'IgnoreCase', true), streamsMerged);
farosACCIdx = cellfun(@(x) contains(x.info.name, 'faros_acc', 'IgnoreCase', true), streamsMerged);

% resample lsl data to faros sampling rate
for id = 1:2
    if id == 1
        farosIdx = farosECGIdx;
    elseif id == 2
        farosIdx = farosACCIdx;
    end

    lslFaros = streams{farosIdx};
    newLSLTime = lslFaros.time_stamps(1):ecgSampRate:lslFaros.time_stamps(end);
    newLSLData = [];
    for ch = 1:size(lslFaros.time_series, 1)
        newChData = interp1(lslFaros.time_stamps, double(lslFaros.time_series(ch,:)), newLSLTime, 'pchip');
        newLSLData = [newLSLData; newChData];
    end
    streamsMerged{farosIdx}.time_stamps = newLSLTime;
    streamsMerged{farosIdx}.time_series = newLSLData;
end
fDataMerged = [];
for f = 1:length(farosData)
    fData = farosData{f};

    % align data from faros and lsl through cross correlation
    [r,lags] = xcorr(streamsMerged{farosECGIdx}.time_series, fData(:,2));
    [~, maxIdx] = max(r(lags >= 0));
    maxIdx = maxIdx + find(lags == 0);
    offset = lags(maxIdx);
    offsetTime = offset*ecgSampRate;

    fData(:,1) = fData(:,1) + streamsMerged{farosECGIdx}.time_stamps(1) + offsetTime;

    fDataMerged = [fDataMerged; fData];
end

% replace ecg
idx = find(fDataMerged(:,1) >= streamsMerged{farosECGIdx}.time_stamps(1) & fDataMerged(:,1) <= streamsMerged{farosECGIdx}.time_stamps(end));
streamsMerged{farosECGIdx}.time_stamps = fDataMerged(idx,1)';
streamsMerged{farosECGIdx}.time_series = fDataMerged(idx,2)';

% replace acc
idx = find(fDataMerged(:,1) >= streamsMerged{farosACCIdx}.time_stamps(1) & fDataMerged(:,1) <= streamsMerged{farosACCIdx}.time_stamps(end));
streamsMerged{farosACCIdx}.time_stamps = fDataMerged(idx,1)';
streamsMerged{farosACCIdx}.time_series = fDataMerged(idx,3:5)';
disp('--- replaced ECG data ---');


%% gaze
% plGazeIdx = contains(streamNames, 'pupil_labs_Gaze', 'IgnoreCase', true);
% plEventIdx = contains(streamNames, 'pupil_labs_Event', 'IgnoreCase', true);
plGazeIdx = cellfun(@(x) contains(x.info.name, 'pupil_labs_Gaze', 'IgnoreCase', true), streamsMerged);
plEventIdx = cellfun(@(x) contains(x.info.name, 'pupil_labs_Event', 'IgnoreCase', true), streamsMerged);

% add lsl times to fixations
plFixData.lsl_times = zeros(size(plFixData,1),1);

for i = 1:size(plFixData,1)
    fixNum = plFixData.fixationId(i);
    fixId = find(plGazeData.fixationId == fixNum);
    fixId = fixId(1); % use start of fixation
    plFixData.lsl_times(i) = plGazeData.lsl_times(fixId);
end

% add lsl times to events
plEventData.lsl_times = zeros(size(plEventData,1),1);

% add known event times from LSL data
for e = 1:length(streamsMerged{plEventIdx}.time_stamps)
    idx = find(strcmp(plEventData.name, streamsMerged{plEventIdx}.time_series(e)));
    if ~isempty(idx)
        plEventData.lsl_times(idx) = streamsMerged{plEventIdx}.time_stamps(e);
    end
end

% calculate times for remaining events from nearest sync point
for e = 1:size(plEventData,1)
    if plEventData.lsl_times(e) == 0
        idx = find(plEventData.lsl_times ~= 0);
        [~,refIdx] = min(abs(idx - e));
        plEventData.lsl_times(e) = plEventData.lsl_times(idx(refIdx)) + (plEventData.timestamp_ns_(e) - plEventData.timestamp_ns_(idx(refIdx)))/10^9;
    end
end

% replace lsl data with pupil cloud data

% gaze:
gazeIdx = find(plGazeData.lsl_times >= streamsMerged{plGazeIdx}.time_stamps(1) & plGazeData.lsl_times <= streamsMerged{plGazeIdx}.time_stamps(end));
gazeDataNew = plGazeData(gazeIdx,:);
streamsMerged{plGazeIdx}.time_stamps = gazeDataNew.lsl_times';
streamsMerged{plGazeIdx}.time_series = [gazeDataNew.gazeXpx, gazeDataNew.gazeYpx, gazeDataNew.fixationId, gazeDataNew.blinkId, gazeDataNew.azimuthdeg, gazeDataNew.elevationdeg]';

streamsMerged{plGazeIdx}.info.desc.channels.channel{3}.label = 'fixationId';
streamsMerged{plGazeIdx}.info.desc.channels.channel{3}.eye = 'both';
streamsMerged{plGazeIdx}.info.desc.channels.channel{4}.label = 'blinkId';
streamsMerged{plGazeIdx}.info.desc.channels.channel{4}.eye = 'both';
streamsMerged{plGazeIdx}.info.desc.channels.channel{5}.label = 'azimuthdeg';
streamsMerged{plGazeIdx}.info.desc.channels.channel{5}.eye = 'both';
streamsMerged{plGazeIdx}.info.desc.channels.channel{6}.label = 'elevationdeg';
streamsMerged{plGazeIdx}.info.desc.channels.channel{6}.eye = 'both';

% fixations (new channel):
fixIdx = find(plFixData.lsl_times >= streamsMerged{plGazeIdx}.time_stamps(1) & plFixData.lsl_times <= streamsMerged{plGazeIdx}.time_stamps(end));
fixDataNew = plFixData(fixIdx,:);
plFixIdx = length(streamsMerged) + 1;
streamsMerged{plFixIdx}.time_stamps = fixDataNew.lsl_times';
streamsMerged{plFixIdx}.time_series = [fixDataNew.fixationId, fixDataNew.durationms, fixDataNew.fixationXpx, fixDataNew.fixationYpx, fixDataNew.azimuthdeg, fixDataNew.elevationdeg]';

streamsMerged{plFixIdx}.info.name = 'pupil_labs_Fixations';
streamsMerged{plFixIdx}.info.type = 'fixations';
streamsMerged{plFixIdx}.info.channel_count = '6';
streamsMerged{plFixIdx}.info.desc.channels.channel{1}.label = 'fixationId';
streamsMerged{plFixIdx}.info.desc.channels.channel{1}.eye = 'both';
streamsMerged{plFixIdx}.info.desc.channels.channel{2}.label = 'durationMs';
streamsMerged{plFixIdx}.info.desc.channels.channel{2}.eye = 'both';
streamsMerged{plFixIdx}.info.desc.channels.channel{3}.label = 'fixationX';
streamsMerged{plFixIdx}.info.desc.channels.channel{3}.eye = 'both';
streamsMerged{plFixIdx}.info.desc.channels.channel{4}.label = 'fixationY';
streamsMerged{plFixIdx}.info.desc.channels.channel{4}.eye = 'both';
streamsMerged{plFixIdx}.info.desc.channels.channel{5}.label = 'azimuthDeg';
streamsMerged{plFixIdx}.info.desc.channels.channel{5}.eye = 'both';
streamsMerged{plFixIdx}.info.desc.channels.channel{6}.label = 'elevationDeg';
streamsMerged{plFixIdx}.info.desc.channels.channel{6}.eye = 'both';

% events:
streamsMerged{plEventIdx}.time_stamps = plEventData.lsl_times';
streamsMerged{plEventIdx}.time_series = horzcat(plEventData.name, plEventData.type)';

streamsMerged{plEventIdx}.info.channel_count = '2';

disp('--- replaced gaze data ---');


%% xsens
% xsensCoMIdx = contains(streamNames, 'CenterOfMass1', 'IgnoreCase', true);
% xsensSegmIdx = contains(streamNames, 'LinearSegmentKinematicsDatagram1', 'IgnoreCase', true);
% xsensJointsIdx = contains(streamNames, 'JointAnglesDatagram1', 'IgnoreCase', true);
xsensCoMIdx = cellfun(@(x) contains(x.info.name, 'CenterOfMass1', 'IgnoreCase', true), streamsMerged);
xsensSegmIdx = cellfun(@(x) contains(x.info.name, 'LinearSegmentKinematicsDatagram1', 'IgnoreCase', true), streamsMerged);
xsensJointsIdx = cellfun(@(x) contains(x.info.name, 'JointAnglesDatagram1', 'IgnoreCase', true), streamsMerged);
if ~isempty(xsensData)
    xsensLSLTime = [];
    xsensCoMData = [];
    xsensSegmentData = [];
    xsensJointData = [];
    for x = 1:length(xsensData.center_of_mass)
        frameRate = xsensData.info{x}{5,2};
        comData = xsensData.center_of_mass{x};
        xsensTime = (0:(size(comData,1)-1)) * (1/frameRate) + streams{xsensCoMIdx}.time_stamps(1);

        % resample existing LSL data to remove time gaps
        lslCoMX = streams{xsensCoMIdx}.time_series(1,:);
        newLSLTime = streams{xsensCoMIdx}.time_stamps(1):1/frameRate:streams{xsensCoMIdx}.time_stamps(end);
        newCoMX = interp1(streams{xsensCoMIdx}.time_stamps, lslCoMX, newLSLTime);

        % use correlation to find offset between offline and LSL data
        [r,lags] = xcorr(newCoMX, comData.CoMPosX);
        [~,maxIdx] = max(r);
        xsensTime = xsensTime + lags(maxIdx)*(1/frameRate);

        % cut off data outside of LSL time
        idx = find(xsensTime >= newLSLTime(1) & xsensTime <= newLSLTime(end));
        xsensLSLTime = [xsensLSLTime xsensTime(idx)];

        tmp = table2array(xsensData.center_of_mass{x});
        tmp(:,1) = [];
        xsensCoMData = [xsensCoMData tmp(idx,:)'];

        tmp = table2array(xsensData.segment_position{x});
        tmp(:,1) = [];
        xsensSegmentData = [xsensSegmentData tmp(idx,:)'];

        tmp = table2array(xsensData.joint_angles{x});
        tmp(:,1) = [];
        xsensJointData = [xsensJointData tmp(idx,:)'];
    end

    streamsMerged{xsensCoMIdx}.time_stamps = xsensLSLTime;
    streamsMerged{xsensCoMIdx}.time_series = xsensCoMData;
    streamsMerged{xsensCoMIdx}.info.desc.channels = {xsensData.center_of_mass{1}.Properties.VariableNames{2:end}};

    streamsMerged{xsensSegmIdx}.time_stamps = xsensLSLTime;
    streamsMerged{xsensSegmIdx}.time_series = xsensSegmentData;
    streamsMerged{xsensSegmIdx}.info.desc.channels = {xsensData.segment_position{1}.Properties.VariableNames{2:end}};

    streamsMerged{xsensJointsIdx}.time_stamps = xsensLSLTime;
    streamsMerged{xsensJointsIdx}.time_series = xsensJointData;
    streamsMerged{xsensJointsIdx}.info.desc.channels = {xsensData.joint_angles{1}.Properties.VariableNames{2:end}};

    disp('--- replaced Xsens data ---');
end


%% add lsl data to eeg struct

% match counter values in eeg and lsl data
[~, idxAmp, idxLSL] = intersect(amp.data(end,:), streams{eegIdx}.time_series);
counterLSLTime = streams{eegIdx}.time_stamps(idxLSL);

newChannelIdx = amp.nbchan+1;

% sample all lsl data at the times of the counter values before adding

% ecg
ampECG = interp1(streamsMerged{farosECGIdx}.time_stamps, streamsMerged{farosECGIdx}.time_series, counterLSLTime, 'cubic');
amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
amp.data(newChannelIdx,idxAmp) = ampECG;
amp.chanlocs(newChannelIdx).labels = 'Faros ECG';
amp.chanlocs(newChannelIdx).type = 'ECG';
newChannelIdx = newChannelIdx + 1;
disp('--- merged ECG into EEG ---');

% photo
% photoIdx = contains(streamNames, 'Photo', 'IgnoreCase', true);
photoIdx = cellfun(@(x) contains(x.info.name, 'Photo', 'IgnoreCase', true), streamsMerged);
ampPhoto = interp1(streamsMerged{photoIdx}.time_stamps, double(streamsMerged{photoIdx}.time_series), counterLSLTime, 'pchip');
amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
amp.data(newChannelIdx, idxAmp) = ampPhoto;
amp.chanlocs(newChannelIdx).labels = 'PhotoSensor';
newChannelIdx = newChannelIdx + 1;
disp('--- merged photo sensor into EEG ---');

% gaze
ampGazeX = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(1,:), counterLSLTime, 'pchip');
ampGazeY = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(2,:), counterLSLTime, 'pchip');
ampFixID = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(3,:), counterLSLTime, 'nearest');
ampBlinkID = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(4,:), counterLSLTime, 'nearest');
amp.data(newChannelIdx:newChannelIdx+3,:) = NaN(4, size(amp.data, 2));
amp.data(newChannelIdx:newChannelIdx+3,idxAmp) = [ampGazeX; ampGazeY; ampFixID; ampBlinkID];
for c = 0:3
    amp.chanlocs(newChannelIdx+c).type = 'Gaze';
end
amp.chanlocs(newChannelIdx).labels = 'GazeX';
amp.chanlocs(newChannelIdx+1).labels = 'GazeY';
amp.chanlocs(newChannelIdx+2).labels = 'FixationID';
amp.chanlocs(newChannelIdx+3).labels = 'BlinkID';
newChannelIdx = newChannelIdx + 4;
disp('--- merged gaze data into EEG ---');

if ~isempty(xsensData)
    % center of mass
    for m = 1:size(streamsMerged{xsensCoMIdx}.time_series,1)
        ampCoM = interp1(streamsMerged{xsensCoMIdx}.time_stamps, streamsMerged{xsensCoMIdx}.time_series(m,:), counterLSLTime, 'pchip');
        amp.data(newChannelIdx+m-1,:) = NaN(1, size(amp.data, 2));
        amp.data(newChannelIdx+m-1, idxAmp) = ampPhoto;
        amp.chanlocs(newChannelIdx+m-1).labels = char(streamsMerged{xsensCoMIdx}.info.desc.channels(m));
        amp.chanlocs(newChannelIdx+m-1).type = 'Xsens';
    end
    newChannelIdx = newChannelIdx + size(streamsMerged{xsensCoMIdx}.time_series,1);
    disp('--- merged center of mass into EEG ---');

    % feet (for step detection)
    leftFootIdx = contains(streamsMerged{xsensSegmIdx}.info.desc.channels, 'LeftFootZ', 'IgnoreCase', true);
    rightFootIdx = contains(streamsMerged{xsensSegmIdx}.info.desc.channels, 'RightFootZ', 'IgnoreCase', true);
    ampLeftFootZ = interp1(streamsMerged{xsensSegmIdx}.time_stamps, streamsMerged{xsensSegmIdx}.time_series(leftFootIdx,:), counterLSLTime, 'pchip');
    ampRightFootZ = interp1(streamsMerged{xsensSegmIdx}.time_stamps, streamsMerged{xsensSegmIdx}.time_series(rightFootIdx,:), counterLSLTime, 'pchip');
    amp.data(newChannelIdx:newChannelIdx+1,:) = NaN(2, size(amp.data, 2));
    amp.data(newChannelIdx:newChannelIdx+1,idxAmp) = [ampLeftFootZ; ampRightFootZ];
    for c = 0:1
        amp.chanlocs(newChannelIdx+c).type = 'Xsens';
    end
    amp.chanlocs(newChannelIdx).labels = char(streamsMerged{xsensSegmIdx}.info.desc.channels(leftFootIdx));
    amp.chanlocs(newChannelIdx+1).labels = char(streamsMerged{xsensSegmIdx}.info.desc.channels(rightFootIdx));
    newChannelIdx = newChannelIdx + 2;
    disp('--- merged foot position into EEG ---');
else
    warning('no Xsens data loaded, merging skipped')
end

% gaze events
% match eeg time and lsl time
lslTimeInEEG = amp.times ./ 1000;
refIdx = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(1));
lslTimeInEEG = lslTimeInEEG - lslTimeInEEG(refIdx) + streams{eegIdx}.time_stamps(1);
for e = 1:size(streamsMerged{plEventIdx}.time_series, 2)
    % find closest matching point for each event
    [m,idx] = min(abs(lslTimeInEEG - streamsMerged{plEventIdx}.time_stamps(e)));
    if m > 0.001
        disp(['Pupil Labs event ' char(streamsMerged{plEventIdx}.time_series(1,e)) ' outside of EEG time']);
        if idx == 1
            m = m*(-1);
        end
        idx = idx + round(m * 1000 / (1000/amp.srate));
    end

    amp.event(e).latency = idx;
    amp.event(e).duration = 1;
    amp.event(e).type = char(streamsMerged{plEventIdx}.time_series(1,e));
end
disp('--- merged gaze events into EEG ---');

% save data
eegMerged = amp;
eegMerged.nbchan = size(eegMerged.data,1);
eegMerged = eeg_checkset(eegMerged);
pop_saveset(eegMerged, 'filename', [subjectFolder '_merged.set'], 'filepath', fullfile(dataPath, subjectFolder));

end




function timealignedgaze = importGaze(filename)
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "sectionId", "recordingId", "timestampns", "gazeXpx", "gazeYpx", "worn", "fixationId", "blinkId", "azimuthdeg", "elevationdeg", "timestamps", "lsl_times"];
opts.VariableTypes = ["double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
% opts = setvaropts(opts, "blinkId", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["sectionId", "recordingId", "blinkId"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["sectionId", "recordingId"], "EmptyFieldRule", "auto");

% Import the data
timealignedgaze = readtable(filename, opts);
end


function fixations = importFixations(filename)
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["sectionId", "recordingId", "fixationId", "startTimestampns", "endTimestampns", "durationms", "fixationXpx", "fixationYpx", "azimuthdeg", "elevationdeg"];
opts.VariableTypes = ["categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["sectionId", "recordingId"], "EmptyFieldRule", "auto");

% Import the data
fixations = readtable(filename, opts);
end

