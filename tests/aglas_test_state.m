function state = aglas_test_state()
%AGLAS_TEST_STATE  Fresh counter struct for the test suite.
    state.n = 0;
    state.pass = 0;
    state.fail = 0;
    state.failures = {};
end
