clc
clear

% Using only one spatial dimension for simplicity, like we will in TriMGRIT
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
global At
global preconditioner

space_steps = 30;
time_steps = 20;

m = ones((space_steps + 1) * time_steps, 1) * 0.1;
rho = ones(space_steps * (time_steps + 1), 1) * 0.5;
lambda = ones(space_steps * (time_steps + 2), 1) * 0.1;
q = zeros(space_steps * (time_steps + 2), 1);

time = 1;
d_time = time/time_steps;
d_space = 1/space_steps;
% Initial and final conditions
%q(1+space_steps/2:space_steps) = ones(space_steps/2, 1);
%q(space_steps * (time_steps + 1)+ 1 : space_steps * (time_steps + 3/2)) = ones(space_steps/2, 1);
%q(space_steps * (time_steps + 3/2) + 1 : space_steps * (time_steps + 2)) = zeros(space_steps/2, 1) + 0.1;

q(1:space_steps) = ones(space_steps, 1)*0.5;
q(space_steps * (time_steps + 1) + 1: space_steps * (time_steps + 2)) = -1 * (ones(space_steps, 1)*0.5 + sin((0:1/(space_steps-1):1) * 2*pi)'*0.5 );

q = q * (1/d_time);

calc_fixed_matrices()

iters = 3;
for i=1:iters
    recalc_matrices(m, rho)

    A = [A_hat, D'; zeros(get_zero_matrix_size()), S];
    b = -preconditioner * [get_GwL(m, rho, lambda); get_GlambdaL(m, rho)];
    
    disp(norm(b)/(space_steps*time_steps)); % Quick and dirty way to check convergence.
    
    solution = A\b;

    dm = solution(1 : (space_steps + 1) * time_steps);
    drho = solution((space_steps + 1) * time_steps + 1 : (space_steps + 1) * time_steps + space_steps * (time_steps + 1));
    dlambda = solution((space_steps + 1) * time_steps + space_steps * (time_steps + 1) + 1 : length(solution));

    alpha = line_search(dm, drho, dlambda);

    m = m + alpha * dm;
    rho = rho + alpha * drho;
    lambda = lambda + alpha * dlambda;
    
end

rho_reshaped = reshape(At * rho, [space_steps, time_steps]);
[X, Y] = meshgrid(1/time_steps:1/time_steps:1, 1/space_steps:1/space_steps:1);
surf(X, Y, rho_reshaped)

function dimens = get_zero_matrix_size()
    global S
    global A_hat
    rows = size(S);
    cols = size(A_hat);
    dimens = [rows(1), cols(2)];
end

function calc_fixed_matrices()
    global D1
    global D2
    global D
    global As
    global At
    D1 = get_derivative_matrix_space();
    D2 = get_derivative_matrix_time();
    D = [D1, D2];
    As = get_As();
    At = get_At();
end

function recalc_matrices(m, rho)
    global A_hat
    global S
    global D
    global preconditioner
    A_hat = get_A_hat(m, rho);
    S = -D / A_hat * D';
    dim_D = size(D);
    preconditioner = [eye(dim_D(2)), zeros([dim_D(2), dim_D(1)]); -D/A_hat, eye(dim_D(1))];
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

function A = get_A_hat(m, rho)
    global As
    global At
    vec_1 = 2 * As' * At * (1./rho);
    vec_2 = 2 * At' * As * (m.^2);
    vec_3 = 1./(rho.^3);
    A = blkdiag(diag(vec_1), diag(vec_2) * diag(vec_3));
end

function GwL = get_GwL(m, rho, lambda)
    global As
    global At
    global D1
    global D2
    GmL = 2 * diag(m) * As' * At * (1./rho) + D1' * lambda;
    GrhoL = -diag(1./(rho.^2)) * At' * As * (m.^2) + D2' * lambda;
    GwL = [GmL; GrhoL];
end

function GlambdaL = get_GlambdaL(m, rho)
    global q
    global D
    GlambdaL = D * [m; rho] - q;
end

function reward = reward(x, dm, drho, dlambda)
    global m
    global rho
    global lambda
    global preconditioner
    newm = m + x * dm;
    newrho = rho + x * drho;
    newlambda = lambda + x * dlambda;
    recalc_matrices(newm, newrho)
    b = -preconditioner*[get_GwL(newm, newrho, newlambda); get_GlambdaL(newm, newrho)];
    reward = norm(b);
end

function alpha = line_search(dm, drho, dlambda)
    f = @(x) reward(x, dm, drho, dlambda);
    options = optimoptions(@fminunc, 'Display', 'none');
    alpha = fminunc(f, 1, options); 
    
%     global trust_radius
% 
%     a = -trust_radius;
%     b = trust_radius;

%     f = @(x) reward(x, dm, drho, dlambda);
%     gr = (1+sqrt(5))/2;
%     
%     c = b - (b - a)/gr;
%     d = a + (b - a)/gr;
%     
%     while abs(c-d) > 1e-5
%         if f(c) < f(d)
%             b = d;
%         else
%             a = c;
%         end
%         c = b - (b - a)/gr;
%         d = a + (b - a)/gr;
%     end
%     
%     alpha = (a+b)/2;
%     fprintf("Alpha: %f\n", alpha)
%     fprintf("f(Alpha): %f\n", f(alpha));
end
