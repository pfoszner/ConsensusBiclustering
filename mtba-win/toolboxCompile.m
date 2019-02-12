function toolboxCompile(sys)
% Compile all required c files in ./private/ folder
% Needs gcc installed on linux or Windows
if nargin<1
    if ispc
        sys = 'Win'
    else
        sys = 'Linux'
    end
end

disp('Compiling ...................................................');

% general compile options
%opts = {'-output'};

% compile c functions
% fs = {'flochelp', 'mbimax',...
%   'mfloc', 'misaBF_down', 'misaBF_up', 'misaBF_updown', 'mprintres'};

% for i=1:length(fs), fs{i}=['./private/' fs{i}]; end
% for i=1:length(fs), mex([fs{i} '.c'],opts{:},[fs{i} '.' mexext]); 
if strcmp(sys, 'Linux') 
    % compile BBC
    mex -lgsl -lgslcblas ./private/mBBC.c ./private/BBCalgorithm.c ./private/BBCkmeans.c -outdir ./private/
    % compile floc
    mex ./private/mprintres.c ./private/flochelp.c -outdir ./private/
    mex ./private/mfloc.c ./private/flochelp.c -outdir ./private/
    % compile bimax
    mex ./private/mbimax.c -outdir ./private/
    % compile isa
    mex ./private/misaBF_up.c -outdir ./private/
    mex ./private/misaBF_down.c -outdir ./private/
    mex ./private/misaBF_updown.c -outdir ./private/
    disp('.............................................. Done Compiling');
else if strcmp(sys, 'Win')
        % compile BBC
        gslpath = input('Give the path to your gsl location (eg: ''C:\gsl''):');
        gsllib = strcat(gslpath, '\lib');
        gsllib = ['-L' gsllib];
        gslinc = strcat(gslpath, '\include');
        gslinc = ['-I' gslinc];
        mex(gslinc, gsllib, '-lgsl', '-lgslcblas', '-outdir', 'private\', 'private\mBBC.c', ...
            'private\BBCalgorithm.c', 'private\BBCkmeans.c');
        % compile floc
        mex('private\mprintres.c' , '-outdir', 'private\', 'private\flochelp.c');
        mex('private\mfloc.c', 'private\flochelp.c', '-outdir', 'private\');
        % compile bimax
        mex('-outdir', 'private\', 'private\mbimax.c');
        % compile isa
        mex('-outdir', 'private\', 'private\misaBF_up.c');
        mex('-outdir', 'private\', 'private\misaBF_down.c');
        mex('-outdir', 'private\', 'private\misaBF_updown.c');
        disp('.............................................. Done Compiling');
else
    disp('Sorry you are on your own :(')
    end 
end

