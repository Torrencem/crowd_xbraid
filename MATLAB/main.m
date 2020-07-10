% Using only one spatial dimension for simplicity, like we will in TriMGRIT
global space_steps
global time_steps
global m
global rho
global lambda
global q
global h

space_steps = 8;
time_steps = 8;
m = rand(space_steps * time_steps, 1); % Could pick different initial conditions
rho = rand(time_steps, 1);
lambda = rand(space_steps * time_steps + time_steps, 1);
q = sparse(space_steps * (time_steps + 1));

time = 1;
h = time/time_steps;

% Initial and final conditions
q(1:space_steps/2) = ones(space_steps/2, 1);
q(space_steps * time_steps + 1 + space_steps/2 : space_steps * time_steps + space_steps) = ones(space_steps/2, 1);

q = q * (1/h);

iters = 5;
for i=1:iters
    b = -[get_GwL(); get_GlambdaL()];
    disp(norm(b)); % Quick and dirty way to check convergence.
    A = [get_A(), (get_D())'; zeros(space_steps * time_steps + time_steps), get_S()];
    solution = linsolve(A, b);
    alpha = line_search(solution);
    m = m + alpha * solution(1 : space_steps * time_steps);
    rho = rho + alpha * solution(space_steps * time_steps + 1 : space_steps * time_steps + time_steps);
    lambda = lambda + alpha * solution(space_steps * time_steps + time_steps + 1 : (space_steps * time_steps + time_steps) * 2);
end

function D = get_derivative_matrix_time()
    global time_steps
    global h
    D = sparse(time_steps);
    D(1, 1) = 2;
    for i=2:time_steps
        D(i, i) = 1;
        D(i, i-1) = -1;
    end
    D = (1/h) * D;
end

function D = get_derivative_matrix_space()
    global space_steps
    global time_steps
    global h
    D = sparse(space_steps);
    D(1, 1) = 2;
    for i=2:space_steps
        D(i, i) = 1;
        D(i, i-1) = -1;
    end
    D = (1/h) * D;
    Drep = repmat({D}, 1, time_steps);
    D = blkdiag(Drep{:});
end

function A = get_A()
    global space_steps
    global time_steps
end

function D = get_D()
    global space_steps
    global time_steps
end

function S = get_S()
    global space_steps
    global time_steps
end

function GwL = get_GwL()
    global space_steps
    global time_steps
end

function GlambdaL = get_GlambdaL()
    global space_steps
    global time_steps
end

function alpha = line_search(solution)
end