if isunix
    [~, username] = system('whoami'); % exists on every unix that I know of
    username =username(1:(end-1));
    if strcmp(username,'brandon'); username='Brandon'; end
elseif ispc
    [~, username] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                         % found it on the net elsewhere, you might want to verify
    cd(sprintf('C:/Users/%s/Documents/MATLAB/',getenv('USERNAME')))
end

% cd(['/hdd/Users/' username '/Dropbox/Research/CFD/MyCode']);
% cd('/mnt/hdd/data/postproc/')

ones(10)*ones(10); %#ok<VUNUS>

% If not running in a terminal, use opengl
if isempty(getenv('TERM'))
    set(0, 'DefaultFigureRenderer', 'OpenGL');
end

% White background for figures
set(0,'defaultfigurecolor',[1 1 1])

% set(0,'defaultfigureposition',[1680 10 560 420]);
% set(0,'DefaultFigurePaperPosition', [1.0 1.0, 6.0, 4.5])
% set(0,'DefaultFigurePaperPositionMode','manual')
set(groot,'defaultLineLineWidth',2)

if exist('mycolororder','file')
    set(0,'DefaultAxesColorOrder',mycolororder())
end



clearvars