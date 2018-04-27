function [ImageData] = readcih(filename)
%   readcih.m 
%   READCIH Read Photron image data into a structure
%   
%   I=READCIH('c:\Photron\Filename') loads metadata from 
%   'c:\Photron\Filename' into structure I.  
%
%   Remarks
%   -------
%   This function must be handed the common *.cih 
%   Author: Brandon Patterson                   Creation Date: April 27,2018
%
%   Examples
%   --------
%   % Load all metadata
%   I=readmraw('c:\Photron\example.cih'); 
filename = 'E:/brandon/research/acoustic_fountain/acoustic_fountain_videos/Db/D1b1_20180424_135525.cih';
fid=fopen(sprintf('%s',filename),'r');
if fid<1; display([filename ' filenames could not be found']); ImageData=0; return; end

tline = fgetl(fid);

ii = 0;
% mystructargs =cell(100,1);
mystructargs={};

while ischar(tline)
    disp(tline)
    tline = fgetl(fid);
    if strcmp(tline,'#Photron Fastcam Viewer Version Information:')
        break
    end
    
    
    linestr = strsplit(tline,':');
    
    % skip blanklines or useless ones
    if length(linestr) <= 1 || isempty(linestr{end})
        continue
    end
    
    myfieldname = strtrim(linestr{1});
    myfieldname = strrep(myfieldname,' ','_');
    myfieldname = strrep(myfieldname,'(','_');
    myfieldname = strrep(myfieldname,')','_');
    myfieldname = strip(myfieldname,'right','_');
    
    
    if ~isvarname(myfieldname) 
       pause 
    end
    
    myfieldvalue = strtrim(strcat(linestr{2:end}));
    [~,isvalueanumber] = str2num(myfieldvalue);
    if isvalueanumber && ~any(myfieldvalue=='/')
        myfieldvalue = str2double(myfieldvalue);
    end
    
    
    if isempty(myfieldvalue)
        continue
    end
    
    ii = ii + 1;
    
    mystructargs = [mystructargs, myfieldname, myfieldvalue];
end

ImageData = struct(mystructargs{1:(2*ii)});


fclose(fid);


% Read Header Information
% Header=textscan(fid1,'%s','delimiter',':');