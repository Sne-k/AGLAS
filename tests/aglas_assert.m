function state = aglas_assert(state, name, condition, detail)
%AGLAS_ASSERT  Record the outcome of one check.
%
%   state = AGLAS_ASSERT(state, name, condition, detail)
%
%   state is a struct with counters and a log. condition is logical. detail is
%   a short string shown on both pass and failure, normally the measured value
%   against the tolerance, so that a passing run still documents the margin.

    if nargin < 4, detail = ''; end

    state.n = state.n + 1;
    if condition
        state.pass = state.pass + 1;
        mark = 'PASS';
    else
        state.fail = state.fail + 1;
        mark = 'FAIL';
        state.failures{end+1} = sprintf('%s : %s', name, detail);
    end
    fprintf('  [%s] %-52s %s\n', mark, name, detail);
end
