function failures = run_all_tests()
%RUN_ALL_TESTS  AGLAS verification suite.
%
%   run_all_tests
%   failures = run_all_tests()
%
%   Runs every check in tests/ and prints a summary. Returns the number of
%   failures, so it can gate a build:
%
%       matlab -batch "addpath('tests'); exit(run_all_tests() > 0)"
%
%   The suite checks the code against closed-form solutions wherever one
%   exists, rather than against previously recorded output. Cantilever
%   frequencies come from the Euler-Bernoulli characteristic equation, static
%   deflection and moment from beam statics, tip twist from m*L^2/(2*GJ),
%   Theodorsen's C(k) from published tables, the turbulence spectrum from its
%   own -5/3 inertial slope, and the divergence speed from the classical
%   uniform-wing result GJ*(pi/2L)^2/(c*e*cl_alpha). A regression suite that
%   only compared against stored numbers would have happily locked in every
%   defect this rewrite removed.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'src', 'common'));
addpath(here);

fprintf('=========================================\n');
fprintf(' AGLAS verification suite\n');
fprintf('=========================================\n');

t_start = tic;
state = aglas_test_state();

suites = {@test_section_and_fem, ...
          @test_modal, ...
          @test_loads_and_stress, ...
          @test_gust_models, ...
          @test_response, ...
          @test_aero_and_flutter};

for i = 1:numel(suites)
    name = func2str(suites{i});
    try
        state = suites{i}(state);
    catch err
        state.n = state.n + 1;
        state.fail = state.fail + 1;
        state.failures{end+1} = sprintf('%s raised: %s', name, err.message);
        fprintf('  [FAIL] %-52s %s\n', [name ' (exception)'], err.message);
    end
end

elapsed = toc(t_start);

fprintf('\n=========================================\n');
fprintf(' %d checks, %d passed, %d failed  (%.1f s)\n', ...
        state.n, state.pass, state.fail, elapsed);
if state.fail > 0
    fprintf('\n Failures:\n');
    for i = 1:numel(state.failures)
        fprintf('   - %s\n', state.failures{i});
    end
else
    fprintf(' All checks passed.\n');
end
fprintf('=========================================\n');

failures = state.fail;
end
