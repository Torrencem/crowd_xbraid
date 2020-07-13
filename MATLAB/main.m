% Using only one spatial dimension for simplicity, like we will in TriMGRIT
global space_steps
global time_steps
global m
global rho
global lambda
global q
global h
global trust_radius

space_steps = 20;
time_steps = 12;
m = ones(space_steps * time_steps, 1) + 0.1;
rho = ones(space_steps * time_steps, 1) + 0.1;
lambda = ones(space_steps * time_steps, 1) + 0.1;
q = sparse(space_steps * time_steps, 1);

time = 1;
h = time/time_steps;

% Initial and final conditions
q(1:space_steps/2) = ones(space_steps/2, 1);
q(space_steps/2+1:space_steps) = ones(space_steps/2, 1) * 0.1;
q(space_steps * (time_steps - 1) + 1 : space_steps * (time_steps - 1) + space_steps/2) = ones(space_steps/2, 1) * 0.1;
q(space_steps * (time_steps-1) + 1 + space_steps/2 : space_steps * time_steps) = ones(space_steps/2, 1);

q = q * (1/h);

trust_radius = 1;

iters = 20;
for i=1:iters
    b = -[get_GwL(m, rho, lambda); get_GlambdaL(m, rho)];
    disp(norm(b)/(space_steps*time_steps)); % Quick and dirty way to check convergence.
    Q = get_A();
    A = [get_A(), (get_D())'; zeros(space_steps * time_steps, space_steps * time_steps * 2), get_S()];
    solution = A\b;
    dm = solution(1 : space_steps * time_steps);
    drho = solution(space_steps * time_steps + 1 : space_steps * time_steps * 2);
    dlambda = solution(space_steps * time_steps * 2 + 1 : space_steps * time_steps * 3);
    alpha = line_search(dm, drho, dlambda);
    m = m + alpha * dm;
    rho = rho + alpha * drho;
    lambda = lambda + alpha * dlambda;
end

function D = get_derivative_matrix_time()
    global space_steps
    global time_steps
    global h
    D = sparse(space_steps * time_steps, space_steps * time_steps);
    for j=1:space_steps
        D(j, j) = 1;
    end
    for i=2:time_steps
        for j=1:space_steps
            D(space_steps * (i-1) + j, space_steps * (i-1) + j) = 1;
            D(space_steps * (i-1) + j, space_steps * (i-2) + j) = -1;
        end
    end
    D = (1/h) * D;
end

function D = get_derivative_matrix_space()
    global space_steps
    global time_steps
    global h
    D = sparse(space_steps, space_steps);
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
    As = sparse(space_steps, space_steps);
    As(1, 1) = 1;
    for i=2:space_steps
        As(i, i) = 1/2;
        As(i, i-1) = 1/2; 
    end
    Asrep = repmat({As}, 1, time_steps);
    As = blkdiag(Asrep{:});
end

function At = get_At()
    global space_steps
    global time_steps
    At = sparse(space_steps * time_steps, space_steps * time_steps);
    for j=1:space_steps
        At(j, j) = 1/2;
    end
    for i=2:time_steps
        for j=1:space_steps
            At(space_steps * (i-1) + j, space_steps * (i-1) + j) = 1/2;
            At(space_steps * (i-1) + j, space_steps * (i-2) + j) = 1/2;
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
    D = [get_derivative_matrix_space(), get_derivative_matrix_time()];
end

function S = get_S()
    D = get_D();
    A = get_A();
    S = -D * inv(A) * D';
end

function GwL = get_GwL(m, rho, lambda)
    As = get_As();
    At = get_At();
    D1 = get_derivative_matrix_space();
    D2 = get_derivative_matrix_time();
    GmM = 2 * diag(m) * As' * At * (1./rho) + D1' * lambda;
    GmRho = -diag(1./(rho.^2)) * At' * As * (m.^2) + D2' * lambda;
    GwL = [GmM; GmRho];
end

function GlambdaL = get_GlambdaL(m, rho)
    global q
    D = get_D();
    GlambdaL = get_D() * [m; rho] - q;
end

function reward = reward(x, dm, drho, dlambda)
    global m
    global rho
    global lambda
    b = -[get_GwL(m+x*dm, rho+x*drho, lambda+x*dlambda); get_GlambdaL(m+x*dm, rho+x*drho)];
    reward = norm(b);
end

function alpha = line_search(dm, drho, dlambda)
    global trust_radius

    a = -trust_radius;
    b = trust_radius;
    f = @(x) reward(x, dm, drho, dlambda);
    gr = (1+sqrt(5))/2;
    
    c = b - (b - a)/gr;
    d = a + (b - a)/gr;
    
    while abs(c-d) > 1e-5
        if f(c) < f(d)
            b = d;
        else
            a = c;
        end
        c = b - (b - a)/gr;
        d = a + (b - a)/gr;
    end
    
    alpha = (a+b)/2;
    fprintf("Alpha: %f\n", alpha)
    fprintf("f(Alpha): %f\n", f(alpha));
end