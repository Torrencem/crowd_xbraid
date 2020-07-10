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
m = rand(space_steps * time_steps, 1);
rho = rand(space_steps * time_steps, 1);
lambda = rand(space_steps * time_steps * 2, 1);
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
    A = [get_A(), (get_D())'; zeros(space_steps * time_steps * 2), get_S()];
    solution = linsolve(A, b);
    alpha = line_search(solution);
    m = m + alpha * solution(1 : space_steps * time_steps);
    rho = rho + alpha * solution(space_steps * time_steps + 1 : space_steps * time_steps * 2);
    lambda = lambda + alpha * solution(space_steps * time_steps * 2 + 1 : space_steps * time_steps * 4);
end

function D = get_derivative_matrix_time()
    global space_steps
    global time_steps
    D = sparse(space_steps * time_steps);
    for j=1:space_steps
        D(j, j) = 1;
    end
    for i=2:time_steps
        for j=1:space_steps
            D(j, j) = 1;
            D(j, j - space_steps) = -1;
        end
    end
    D = (1/h) * D;
end

function D = get_derivative_matrix_space()
    global space_steps
    global time_steps
    global h
    D = sparse(space_steps);
    D(1, 1) = 1;
    for i=2:space_steps
        D(i, i) = 1;
        D(i, i-1) = -1;
    end
    D = (1/h) * D;
    Drep = repmat({D}, 1, time_steps);
    D = blkdiag(Drep{:});
end

function As = get_As()
    global space_steps
    global time_steps
    As = sparse(space_steps);
    As(1, 1) = 1;
    for i=2:space_steps
        As(i, i) = 1/2;
        As(i, i-1) = 1/2;
    end
    Asrep = repmat({A}, 1, time_steps);
    As = blkdiag(Asrep{:});
end

function At = get_At()
    global space_steps
    global time_steps
    global space_steps
    global time_steps
    At = sparse(space_steps * time_steps);
    for j=1:space_steps
        At(j, j) = 1/2;
    end
    for i=2:time_steps
        for j=1:space_steps
            At(j, j) = 1/2;
            At(j, j - space_steps) = 1/2;
        end
    end
end

function A = get_A()
    global rho
    global m
    As = get_As();
    At = get_At();
    vec_1 = 2 * As' * At * (1./rho);
    vec_2 = 2 * At' * As * (m.^2);
    vec_3 = 1./(rho.^3);
    A = blkdiag(diag(vec_1), diag(vec_2) * diag(vec_3));
end

function D = get_D()
    D = blkdiag(get_derivative_matrix_space(), get_derivative_matrix_time());
end

function S = get_S()
    D = get_D();
    A = get_A();
    S = -D' * inv(A) * D;
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