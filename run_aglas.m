function run_aglas(which_pipeline, run_tests)
%RUN_AGLAS  Run the AGLAS analysis chain end to end.
%
%   run_aglas                     run both pipelines, tests first
%   run_aglas('main')             main pipeline only     (c1p ... c6p)
%   run_aglas('supplementary')    supplementary only     (c1m ... c6m)
%   run_aglas('both')             both, in order
%   run_aglas('tests')            verification suite only
%   run_aglas(pipeline, false)    skip the verification suite
%
%   Each stage runs in its own function workspace, so no stage can leave
%   variables behind for the next one to trip over. That was not possible
%   before: none of the original scripts cleared its workspace, and c4p.m used
%   i as a loop counter over 241 airspeeds, so any driver that called the
%   stages in sequence had its own loop variable overwritten on the first pass.
%
%   Stage order and data flow
%     main            c1p -> wing_geom      -> c2p -> structural_model
%                                           -> c3p -> modal_response
%                                           -> c4p -> aero_data
%                                           -> c5p -> animation GIF
%                                           -> c6p -> fos_data
%     supplementary   c1m -> gust_models
%                     c2m -> mesh_convergence
%                     c3m -> turbulence_response   (needs structural_model)
%                     c4m -> stress_field          (needs modal_response)
%                     c5m -> strain_analysis
%                     c6m -> fos_optimisation
%
%   The supplementary pipeline reads the structural model and gust response
%   produced by the main pipeline, so run 'main' first, or use 'both'.

    if nargin < 1 || isempty(which_pipeline), which_pipeline = 'both'; end
    if nargin < 2 || isempty(run_tests),      run_tests = true;        end

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src', 'common'));
    addpath(fullfile(root, 'src', 'main_pipeline'));
    addpath(fullfile(root, 'src', 'supplementary_pipeline'));
    addpath(fullfile(root, 'tests'));

    which_pipeline = lower(strtrim(which_pipeline));

    banner('AGLAS');
    fprintf(' Gust load response analysis of a flexible wing\n');
    fprintf(' started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    t0 = tic;
    failed = {};

    % ------------------------------------------------------------ tests --
    if run_tests || strcmp(which_pipeline, 'tests')
        banner('Verification suite');
        n_fail = run_all_tests();
        if n_fail > 0
            error('run_aglas:testsFailed', ...
                  ['%d verification checks failed. The physics library does ' ...
                   'not match its closed-form references, so the pipeline ' ...
                   'results would not be trustworthy. Fix these first.'], n_fail);
        end
    end
    if strcmp(which_pipeline, 'tests')
        fprintf('\nTests only. Done in %.1f s.\n', toc(t0));
        return;
    end

    % --------------------------------------------------------- pipelines --
    main_stages = {'c1p','c2p','c3p','c4p','c5p','c6p'};
    supp_stages = {'c1m','c2m','c3m','c4m','c5m','c6m'};

    switch which_pipeline
        case 'main',          stages = main_stages;
        case 'supplementary', stages = supp_stages;
        case 'both',          stages = [main_stages, supp_stages];
        otherwise
            error('run_aglas:unknownPipeline', ...
                  'Unknown pipeline "%s". Use main, supplementary, both or tests.', ...
                  which_pipeline);
    end

    for i = 1:numel(stages)
        name = stages{i};
        banner(sprintf('Stage %d of %d : %s', i, numel(stages), name));
        ts = tic;
        try
            run_stage(name);
            fprintf('\n[ok] %s finished in %.1f s\n', name, toc(ts));
        catch err
            failed{end+1} = sprintf('%s : %s', name, err.message); %#ok<AGROW>
            fprintf('\n[FAILED] %s : %s\n', name, err.message);
        end
    end

    % ----------------------------------------------------------- summary --
    banner('Summary');
    paths = aglas_paths();
    fprintf(' %d stages run in %.1f s\n', numel(stages), toc(t0));
    if isempty(failed)
        fprintf(' All stages completed.\n');
    else
        fprintf(' %d stage(s) failed:\n', numel(failed));
        for i = 1:numel(failed)
            fprintf('   - %s\n', failed{i});
        end
    end
    fprintf('\n Data    : %s\n', paths.data);
    fprintf(' Results : %s\n', paths.results);
    fprintf(' Frames  : %s\n', paths.animation);

    report_headline(paths);
end

% ------------------------------------------------------------------------
function run_stage(stage_name_)
%RUN_STAGE  Execute one pipeline script in an isolated function workspace.
%
%   A script called from inside a function runs in that function's workspace,
%   not the caller's. Every variable the stage creates therefore dies with this
%   call, which is what lets the stages be chained safely. The trailing
%   underscore keeps this function's own variable out of the stage's way.
    eval(stage_name_);
end

% ------------------------------------------------------------------------
function banner(text)
    fprintf('\n');
    fprintf('%s\n', repmat('=', 1, 70));
    fprintf(' %s\n', text);
    fprintf('%s\n', repmat('=', 1, 70));
end

% ------------------------------------------------------------------------
function report_headline(paths)
%REPORT_HEADLINE  Print the numbers a reader of this project actually wants.

    fprintf('\n Headline results\n');
    fprintf(' %s\n', repmat('-', 1, 68));

    show = @(label, fmt, value) fprintf('   %-34s %s\n', label, sprintf(fmt, value));

    f = fullfile(paths.data, 'structural_model.mat');
    if exist(f, 'file')
        d = load(f, 'modes');
        show('first bending frequency', '%.3f Hz', d.modes.freq_hz(1));
        it = find(strcmp(d.modes.type, 'torsion'), 1);
        if ~isempty(it)
            show('first torsion frequency', '%.3f Hz', d.modes.freq_hz(it));
        end
    end

    f = fullfile(paths.data, 'modal_response.mat');
    if exist(f, 'file')
        d = load(f, 'tip_w', 'cfg');
        show('peak tip deflection', '%.4f m', max(abs(d.tip_w)));
        show('  as a fraction of span', '%.2f %%', ...
             100*max(abs(d.tip_w))/d.cfg.geom.L);
    end

    f = fullfile(paths.data, 'aero_data.mat');
    if exist(f, 'file')
        d = load(f, 'flutter_speed', 'flutter_freq', 'divergence_speed');
        show('flutter speed', '%.1f m/s', d.flutter_speed);
        show('flutter frequency', '%.2f Hz', d.flutter_freq);
        show('divergence speed', '%.1f m/s', d.divergence_speed);
    end

    f = fullfile(paths.data, 'fos_data.mat');
    if exist(f, 'file')
        d = load(f, 'peak_root_stress', 'min_FoS');
        show('peak root stress', '%.2f MPa', d.peak_root_stress/1e6);
        show('minimum factor of safety', '%.3f', d.min_FoS);
    end

    f = fullfile(paths.data, 'fos_optimisation.mat');
    if exist(f, 'file')
        d = load(f, 'best_factor', 'base', 'opt');
        show('sized wall thickness factor', '%.4f', d.best_factor);
        show('resulting mass change', '%+.2f %%', ...
             100*(d.opt.mass - d.base.mass)/d.base.mass);
    end
    fprintf(' %s\n', repmat('-', 1, 68));
end
