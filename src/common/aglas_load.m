function aglas_load(file, varargin)
%AGLAS_LOAD  Load stage outputs into the caller, failing loudly if anything is missing.
%
%   AGLAS_LOAD(file, 'var1', 'var2', ...)
%
%   Behaves like load(file, 'var1', ...) called from the calling script, with
%   two differences. A missing file names the stage that writes it, and a
%   requested variable absent from the file is an error. Plain load silently
%   omits a missing variable, so a renamed save in one stage only surfaced as
%   an "undefined variable" error somewhere in the next stage, which is the
%   class of defect the original scripts suffered from.

    if ~exist(file, 'file')
        [~, base] = fileparts(file);
        error('aglas_load:missingFile', ...
              '%s not found. It is written by stage %s; run that stage first (run_aglas runs everything in order).', ...
              file, producer(base));
    end

    % load's own warning for a missing variable is replaced by the error below.
    w = warning('off', 'MATLAB:load:variableNotFound');
    restore = onCleanup(@() warning(w));
    S = load(file, varargin{:});
    missing = setdiff(varargin, fieldnames(S));
    if ~isempty(missing)
        [~, base] = fileparts(file);
        error('aglas_load:missingVariable', ...
              '%s has no variable(s) %s. It may be left over from an older version; rerun stage %s.', ...
              file, strjoin(missing, ', '), producer(base));
    end

    names = fieldnames(S);
    for i = 1:numel(names)
        assignin('caller', names{i}, S.(names{i}));
    end
end

function s = producer(base)
    map = struct('wing_geom', 'c1p', 'structural_model', 'c2p', ...
                 'modal_response', 'c3p', 'stress_field', 'c4m', ...
                 'strain_analysis', 'c5m');
    if isfield(map, base)
        s = map.(base);
    else
        s = '(unknown)';
    end
end
