function convertXsens(subject, dataPath)
% convert xsens data from .xlsx to Matlab .mat format for faster reading

subjectID = ['WaS_' sprintf('%03d', subject)];
files = dir(fullfile(dataPath, subjectID, ['*' num2str(subject) '*.xlsx']));

sheets = {'General Information', 'Center of Mass', 'Segment Position', 'Joint Angles ZXY'};
varNames = {'info', 'center_of_mass', 'segment_position', 'joint_angles'};

xsensData = struct;

for v = 1:length(varNames)
    xsensData.(varNames{v}) = {};
end

for f = 1:length(files)
    for s = 1:length(sheets)
        if s == 1
            xsensData.(varNames{s}){f} = readcell(fullfile(files(f).folder, files(f).name), 'Sheet', sheets{s});
        else
            xsensData.(varNames{s}){f} = readtable(fullfile(files(f).folder, files(f).name), 'Sheet', sheets{s});
        end
        disp(['file ' files(f).name ', sheet ' sheets{s} ' done']);
    end
end

% for f = 1:length(files)
%     ssds = spreadsheetDatastore(fullfile(files(f).folder, files(f).name));
%     sn = sheetnames(ssds,1);
%     for s = 1:length(sheets)
%         if s == 1
%             xsensData.(varNames{s}){f} = readcell(fullfile(files(f).folder, files(f).name), 'Sheet', sheets{s});
%         else
%             idx = find(strcmp(sn, sheets{s}));
%             ssds.Sheets = idx;
%             xsensData.(varNames{s}){f} = readall(ssds);
%         end
%         disp(['file ' files(f).name ', sheet ' sheets{s} ' done']);
%     end
% end

save(fullfile(dataPath, subjectID, [subjectID '_xsens.mat']), 'xsensData');
end