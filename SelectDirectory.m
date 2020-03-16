function dirname = SelectDirectory()
%% create list of measurements
directories = dir('*_*_*');
dirnames    = {directories.name};

%% check if measurements have been pre-analysed to identify soil-spectra
preanalysed = zeros(1,length(dirnames));
for j=1:length(dirnames)
    dirname = dirnames{j};
    files = dir([dirname,'/soil/*.sed']);
    if ~isempty(files)
        preanalysed(j) = 1;
    end
end
dirnames    = dirnames(find(preanalysed));

% ierror  =   isnan(soilspec_)
% interp1(wl(ierror),soilspec(ierror),wl)

%% Choose for which measurements a LUT should be created
[s,v]           =   listdlg('PromptString','Select a file:', 'SelectionMode','single', 'ListString',dirnames);
dirname         =   dirnames{s};