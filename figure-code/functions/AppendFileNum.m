function outputFilename = AppendFileNum(FileName)
% Given a target filepath (FileName), it appends a number to it to create a
% new filename that does not exist in the target dir.
% If file doesn't exist, it just returns the input

if ~isfile(FileName)
    outputFilename = FileName;
    return
end    

[fPath, fName, fExt] = fileparts(FileName);
if isempty(fExt)  % No '.mat' in FileName
  fExt     = '.mat';
  FileName = fullfile(fPath, [fName, fExt]);
end

ii = 1;
while true
    possibleFile = fullfile(fPath, [fName, '_', num2str(ii), fExt]);
    if ~isfile(possibleFile)
        outputFilename = possibleFile;
        return
    end
    ii = ii+1;
end

end