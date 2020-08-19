% A MATLAB version of Haber and Horesh's algorithm for solving crowd model
% problems, useful for rapid prototyping and getting an idea of the results
% that our TriMGRIT implementation should be aiming for.
% 
% The original algorithm can be found in H and H's 2015 paper "A Multilevel
% Method for the Solution of Time Dependent Optimal Transport".
% 
% Using only one spatial dimension for simplicity, like we will in
% TriMGRIT.

%% Global variables, options, and one-line helper functions

global space_steps
global time_steps
global m
global rho
global lambda
global q
global d_time
global d_space
global A_hat
global D
global S

iters = 10;

space_steps = 8;
time_steps = 8;

time = 1;
d_time = time/time_steps;
d_space = 1/space_steps;

showrho = @(x) surf(reshape(x, [space_steps, time_steps+1]));
showm = @(x) surf(reshape(x, [space_steps+1, time_steps]));
showlambda = @(x) surf(reshape(x, [space_steps, time_steps+2]));

expodist = @(x, k) (2^(1/(1+(x-k)^2)-1))^80; % A nice curve to use when building initial and final conditions.

%% Set up m, rho, lambda, and initial conditions

m = zeros((space_steps + 1) * time_steps, 1);
rho = ones(space_steps * (time_steps + 1), 1) * 0.6;
lambda = ones(space_steps * (time_steps + 2), 1) * 0.1;
q = zeros(space_steps * (time_steps + 2), 1);

q(1:space_steps) = arrayfun(@(x) 0.8*expodist(x, 0.3) + 0.3*expodist(x, 0.7) +0.1, 0:1/(space_steps-1):1);
q(space_steps * (time_steps + 1) + 1: space_steps * (time_steps + 2)) = -1*arrayfun(@(x) 0.1*expodist(x, 0.3)+0.8*expodist(x, 0.7), 0:1/(space_steps-1):1);

q = q * (1/d_time);

%% Set up constant matrices

calc_fixed_matrices()

%% Check conservation of mass

if abs(sum(q(1:space_steps)) + sum(q(space_steps * (time_steps + 1) + 1: space_steps * (time_steps + 2)))) > 1e-5
    disp("Conservation of mass is false.  This is not a problem, but might lead to seemingly strange behaviour.");
end


%% Main loop

for i=1:iters
    normcoeff = 0.0001; % Otherwise known as alpha_k

    recalc_matrices(m, rho, normcoeff);
    sizeD = size(D);
    b = D*inv(A_hat)*get_GwL(m, rho, lambda, normcoeff) - get_GlambdaL(m, rho);
    disp(norm(b)/(space_steps*time_steps)); % Quick way to check convergence

    dlambda = S\b;

    dw = -inv(A_hat)*(D'*dlambda+get_GwL(m, rho, lambda, normcoeff));
    dm = dw(1 : (space_steps + 1) * time_steps);
    drho = dw((space_steps + 1) * time_steps + 1 : (space_steps + 1) * time_steps + space_steps * (time_steps + 1));

    alpha = line_search(dm, drho, dlambda, normcoeff);
    m = m + alpha * dm;
    rho = rho + alpha * drho;
    lambda = lambda + alpha * dlambda;
end

%% Output
showrho(rho)

%% Helper functions: generating matrices

function calc_fixed_matrices()
    global D1
    global D2
    global D
    global As
    global At
    global filler_zeros
    global K
    global X
    D1 = get_derivative_matrix_space();
    D2 = get_derivative_matrix_time();
    D = [D1, D2];
    As = get_As();
    At = get_At();
    dim_d = size(D);
    filler_zeros = zeros(dim_d(1));

    K = zeros(space_steps, space_steps+1);
    for i=1:space_steps
        K(i, i) = -1;
        K(i, i+1) = 1;
    end
    K = K/d_space;

    X = zeros(space_steps+1, space_steps);
    for i=1:space_steps
        X(i, i) = 1/4;
        X(i+1, i) = 1/4;
    end
end

function recalc_matrices(m, rho, normcoeff)
    global A_hat
    global S
    A_hat = get_A_hat(m, rho, normcoeff);
    S = get_S(m, rho, normcoeff);
end

function S = get_S(m, rho, normcoeff)
    global K
    global space_steps
    global time_steps
    global d_time
    S = sparse(1,1,0,space_steps*(time_steps+2),space_steps*(time_steps+2));
    S(1:space_steps, 1:space_steps) = inv(Q(0, m, rho, normcoeff))/d_time^2;
    S(1:space_steps, space_steps+1:2*space_steps) = -inv(Q(0, m, rho, normcoeff))/d_time^2;
    for i=1:time_steps
        S(i*space_steps+1:(i+1)*space_steps, (i-1)*space_steps+1:i*space_steps) = -inv(Q(i-1, m, rho, normcoeff))/d_time^2;
        S(i*space_steps+1:(i+1)*space_steps, (i+1)*space_steps+1:(i+2)*space_steps) = -inv(Q(i, m, rho, normcoeff))/d_time^2;
        S(i*space_steps+1:(i+1)*space_steps, i*space_steps+1:(i+1)*space_steps) = inv(Q(i-1, m, rho, normcoeff))/d_time^2+inv(Q(i, m, rho, normcoeff))/d_time^2+K*inv(P(i-1, rho))*K';
    end
    S((time_steps+1)*space_steps+1:end, time_steps*space_steps+1:(time_steps+1)*space_steps) = -inv(Q(time_steps, m, rho, normcoeff))/d_time^2;
    S((time_steps+1)*space_steps+1:end, (time_steps+1)*space_steps+1:end) = inv(Q(time_steps, m, rho, normcoeff))/d_time^2;
    S = -S;
end

function Pi = P(i, rho) %P_i
    global X
    global space_steps
    Pi = 2 * spdiags(X*(1./(rho(i*space_steps+1:(i+1)*space_steps)) + 1./(rho((i+1)*space_steps+1:(i+2)*space_steps))),0,space_steps+1,space_steps+1);
end

function Qi = Q(i, m, rho, normcoeff) %Q_i
    global X
    global space_steps
    global time_steps
    if i == 0
        Qi = 2 * spdiags(X'*(m(i*(space_steps+1)+1:(i+1)*(space_steps+1)).^2),0,space_steps,space_steps) * spdiags(rho(i*space_steps+1:(i+1)*space_steps),0,space_steps,space_steps);
    elseif i == time_steps
        Qi = 2 * spdiags(X'*(m((i-1)*(space_steps+1)+1:i*(space_steps+1)).^2),0,space_steps,space_steps) * spdiags(rho(i*space_steps+1:(i+1)*space_steps),0,space_steps,space_steps);
    else
        Qi = 2 * spdiags(X'*(m(i*(space_steps+1)+1:(i+1)*(space_steps+1)).^2+m((i-1)*(space_steps+1)+1:i*(space_steps+1)).^2),0,space_steps,space_steps) * spdiags(rho(i*space_steps+1:(i+1)*space_steps),0,space_steps,space_steps);
    end
    Qi = Qi + normcoeff*eye(size(Qi));
end


function D = get_derivative_matrix_space()
    global space_steps
    global time_steps
    global d_space
    top_bottom = sparse(space_steps, (space_steps + 1) * time_steps);
    interior_block = sparse(space_steps, space_steps + 1);
    for i=1:space_steps
        interior_block(i, i) = -1;
        interior_block(i, i+1) = 1;
    end
    interior_cell = repmat({interior_block}, 1, time_steps);
    D = [top_bottom; blkdiag(interior_cell{:}); top_bottom];  
    D = (1/d_space) * D;
end


function D = get_derivative_matrix_time()
    global space_steps
    global time_steps
    global d_time
    top = sparse(space_steps, space_steps * (time_steps + 1));
    bottom = sparse(space_steps, space_steps * (time_steps + 1));
    for i=1:space_steps
        top(i, i) = 1;
        bottom(i, space_steps * time_steps + i) = -1;
    end
    
    center = sparse(space_steps * time_steps, space_steps * (time_steps + 1));
    for i=1:time_steps
        for j=1:space_steps
            center((i-1) * space_steps + j, (i-1) * space_steps + j) = -1;
            center((i-1) * space_steps + j, i * space_steps + j) = 1;
        end
    end
    
    D = [top; center; bottom];
    
    D = (1/d_time) * D;
end

function As = get_As()
    global space_steps
    global time_steps
    As = sparse(space_steps, space_steps + 1);
    for i=1:space_steps
        As(i, i) = 1/2;
        As(i, i+1) = 1/2; 
    end
    Asrep = repmat({As}, 1, time_steps);
    As = blkdiag(Asrep{:});
end

function At = get_At()
    global space_steps
    global time_steps
    At = sparse(space_steps * time_steps, space_steps * (time_steps + 1));
    for i=1:time_steps
        for j=1:space_steps
            At(space_steps * (i-1) + j, space_steps * (i-1) + j) = 1/2;
            At(space_steps * (i-1) + j, space_steps * (i) + j) = 1/2;
        end
    end
end

function A = get_A_hat(m, rho, normcoeff)
    global As
    global At
    vec_1 = 2 * As' * At * (1./rho);
    vec_2 = 2 * At' * As * (m.^2);
    vec_3 = 1./(rho.^3);
    A = blkdiag(spdiags(vec_1,0,length(vec_1),length(vec_1)), spdiags(vec_2,0,length(vec_2),length(vec_2)) * spdiags(vec_3,0,length(vec_3),length(vec_3)) + normcoeff * eye(length(vec_2)));
end

%% Helper functions: calculate the gradient of the Lagrangian

function GwL = get_GwL(m, rho, lambda, normcoeff)
    global As
    global At
    global D1
    global D2
    global Gm
    global Grho
    GmL = 2 * spdiags(m,0,length(m),length(m)) * As' * At * (1./rho) + D1' * lambda;
    GrhoL = -spdiags(1./(rho.^2),0,length(rho),length(rho)) * At' * As * (m.^2) + D2' * lambda + 2*normcoeff*rho;
    Gm = GmL;
    Grho = GrhoL;
    GwL = [GmL; GrhoL];
end

function GlambdaL = get_GlambdaL(m, rho)
    global q
    global D
    global Glambda
    GlambdaL = D * [m; rho] - q;
    Glambda = GlambdaL;
end

%% Helper functions: line search
function reward = reward(x, dm, drho, dlambda, normcoeff)
    global m
    global rho
    global lambda
    newm = m + x * dm;
    newrho = rho + x * drho;
    newlambda = lambda + x * dlambda;
    %recalc_matrices(newm, newrho, normcoeff)
    b = -[get_GwL(newm, newrho, newlambda, normcoeff); get_GlambdaL(newm, newrho)];
    reward = norm(b);
end

function alpha = line_search(dm, drho, dlambda, normcoeff)
    trust_radius = 1;
    a = -trust_radius;
    b = trust_radius;

    f = @(x) reward(x, dm, drho, dlambda, normcoeff);
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
