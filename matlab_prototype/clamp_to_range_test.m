%% Test 1: Check to clamp to right boundary
N = 2;
boundary = [(N + 1) (N + 1) 1 1];
assert(sum(clamp_to_range([3 1], boundary) - [3 1]) <= 1e-10, 'Right boundary problem')

assert(sum(clamp_to_range([4 1], boundary) - [3 1]) <= 1e-10, 'Right boundary problem')

assert(sum(clamp_to_range([4 2], boundary) - [3 2]) <= 1e-10, 'Right boundary problem')

%% Test 2: Clamp to left boundary
N = 2;
boundary = [(N + 1) (N + 1) 1 1];

assert(sum(clamp_to_range([1 1], boundary) - [1 1]) <= 1e-10, 'Left boundary problem')

assert(sum(clamp_to_range([0 1], boundary) - [1 1]) <= 1e-10, 'Left boundary problem')

assert(sum(clamp_to_range([0 2], boundary) - [1 2]) <= 1e-10, 'Left boundary problem')

%% Test 3: Clamp to top boundary
N = 2;
boundary = [(N + 1) (N + 1) 1 1];

assert(sum(clamp_to_range([3 3], boundary) - [3 3]) <= 1e-10, 'Top boundary problem')

assert(sum(clamp_to_range([1 4], boundary) - [1 3]) <= 1e-10, 'Top boundary problem')

assert(sum(clamp_to_range([2 4], boundary) - [2 3]) <= 1e-10, 'Top boundary problem')

%% Test 4: Clamp to bottom boundary
N = 2;
boundary = [(N + 1) (N + 1) 1 1];

assert(sum(clamp_to_range([2 1], boundary) - [2 1]) <= 1e-10, 'Right boundary problem')

assert(sum(clamp_to_range([1 0], boundary) - [1 1]) <= 1e-10, 'Right boundary problem')

assert(sum(clamp_to_range([2 -1], boundary) - [2 1]) <= 1e-10, 'Right boundary problem')