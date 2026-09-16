function p = aglas_paths()
%AGLAS_PATHS  Resolve project directories independently of the working folder.
%
%   p = AGLAS_PATHS() returns a struct with absolute paths to the project root
%   and its data/, results/ and animation/ subfolders, creating them if they do
%   not yet exist, and puts src/common on the MATLAB path.
%
%   Every script calls this instead of using bare filenames such as
%   load('wing_geom.mat'). The original code only worked when the current
%   folder happened to be the one holding the .mat files.

    here = fileparts(mfilename('fullpath'));      % <root>/src/common
    root = fileparts(fileparts(here));            % <root>

    p.root      = root;
    p.src       = fullfile(root, 'src');
    p.common    = fullfile(root, 'src', 'common');
    p.data      = fullfile(root, 'data');
    p.results   = fullfile(root, 'results');
    p.animation = fullfile(root, 'animation');
    p.tests     = fullfile(root, 'tests');

    if ~any(strcmp(p.common, strsplit(path, pathsep)))
        addpath(p.common);
    end

    folders = {p.data, p.results, p.animation};
    for i = 1:numel(folders)
        if ~exist(folders{i}, 'dir')
            mkdir(folders{i});
        end
    end
end
