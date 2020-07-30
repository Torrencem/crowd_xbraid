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
global At
global As
global filler_zeros

space_steps = 100;
time_steps = 100;
show = @(x) surf(reshape(x, [space_steps, time_steps+1]));
showm = @(x) surf(reshape(x, [space_steps+1, time_steps]));
showl = @(x) surf(reshape(x, [space_steps, time_steps+2]));

m = zeros((space_steps + 1) * time_steps, 1);
rho = ones(space_steps * (time_steps + 1), 1) * 0.5;
lambda = ones(space_steps * (time_steps + 2), 1) * 0.5;
q = zeros(space_steps * (time_steps + 2), 1);

expodist = @(x, k) (2^(80*(1/(1+(x-k)^2)-1)));

time = 1;
d_time = time/time_steps;
d_space = 1/space_steps;
% Initial and final conditions
%q(1+space_steps/2:space_steps) = ones(space_steps/2, 1);
%q(space_steps * (time_steps + 1)+ 1 : space_steps * (time_steps + 3/2)) = ones(space_steps/2, 1);
%q(space_steps * (time_steps + 3/2) + 1 : space_steps * (time_steps + 2)) = zeros(space_steps/2, 1) + 0.1;

q(1:space_steps) = arrayfun(@(x) 0.8*expodist(x, 0.3) + 0.3*expodist(x, 0.7), 0:1/(space_steps-1):1);
q(space_steps * (time_steps + 1) + 1: space_steps * (time_steps + 2)) = -1*arrayfun(@(x) 0.3*expodist(x, 0.3) + 0.8*expodist(x, 0.7), 0:1/(space_steps-1):1);

q = q * (1/d_time);

calc_fixed_matrices()

if abs(sum(q(1:space_steps)) + sum(q(space_steps * (time_steps + 1) + 1: space_steps * (time_steps + 2)))) > 1e-5
    disp("Conservation of mass is false.  This is not a problem, but might lead to seemingly strange behaviour.");
end

iters = 10;
for i=1:iters
    recalc_matrices(m, rho)

    A = [A_hat, D'; D, filler_zeros];
    b = -[get_GwL(m, rho, lambda); get_GlambdaL(m, rho)];
    
    disp(norm(b)/(space_steps*time_steps)); % Quick and dirty way to check convergence.
    
    solution = A\b;
    dm = solution(1 : (space_steps + 1) * time_steps);
    drho = solution((space_steps + 1) * time_steps + 1 : (space_steps + 1) * time_steps + space_steps * (time_steps + 1));
    dlambda = solution((space_steps + 1) * time_steps + space_steps * (time_steps + 1) + 1 : length(solution));

    alpha = line_search(dm, drho, dlambda);
    
    m = m + alpha * dm;
    m = m * 1;
    rho = rho + alpha * drho;
    lambda = lambda + alpha * dlambda;
end
show(rho)

function calc_fixed_matrices()
    global D1
    global D2
    global D
    global As
    global At
    global filler_zeros
    D1 = get_derivative_matrix_space();
    D2 = get_derivative_matrix_time();
    D = [D1, D2];
    As = get_As();
    At = get_At();
    dim_d = size(D);
    filler_zeros = zeros(dim_d(1));
end

function recalc_matrices(m, rho)
    global A_hat
    A_hat = get_A_hat(m, rho);
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
    A = spdiags([vec_1;vec_2.*vec_3],0,length(vec_1)+length(vec_2),length(vec_1)+length(vec_2));
end

function GwL = get_GwL(m, rho, lambda)
    global As
    global At
    global D1
    global D2
    GmL = 2 * spdiags(m,0,length(m),length(m)) * As' * At * (1./rho) + D1' * lambda;
    GrhoL = -spdiags(1./(rho.^2),0,length(rho),length(rho)) * At' * As * (m.^2) + D2' * lambda;
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
    newm = m + x * dm;
    newrho = rho + x * drho;
    newlambda = lambda + x * dlambda;
    %recalc_matrices(newm, newrho)
    b = -[get_GwL(newm, newrho, newlambda); get_GlambdaL(newm, newrho)];
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
%     end5
%     
%     alpha = (a+b)/2;
%     fprintf("Alpha: %f\n", alpha)
%     fprintf("f(Alpha): %f\n", f(alpha));
end
