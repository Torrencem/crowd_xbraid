%% Initial variables
clear
load("Gvalues");

global X
global rho
global m
global space_steps
global time_steps


space_steps = 4;
time_steps = 4;
dx = 1/4;
dt = 1/4;

m = zeros((space_steps + 1) * time_steps, 1);
rho = ones(space_steps * (time_steps + 1), 1) * 0.5;
lambda = ones(space_steps * (time_steps + 2), 1) * 0.5;
q = zeros(space_steps * (time_steps + 2), 1);

dm = sym('dm', [(space_steps + 1) * time_steps 1]);
drho = sym('drho', [space_steps * (time_steps + 1) 1]);
dlambda = sym('dlambda', [space_steps * (time_steps + 2) 1]);

%% Matrices
K = zeros(space_steps, space_steps+1);
for i=1:space_steps
    K(i, i) = -1;
    K(i, i+1) = 1;
end
K = K/dx;

X = zeros(space_steps+1, space_steps);
for i=1:space_steps
    X(i, i) = 1/4;
    X(i+1, i) = 1/4;
end

%% Setting up equations
eqns = [ (1/dt)*drho(1:space_steps) == -Glambda(1:space_steps) ];
for i=0:(time_steps-1)
    eqns = [eqns;
        P(i) * dm(i*(space_steps+1)+1:(i+1)*(space_steps+1)) + K' * dlambda((i+1)*space_steps+1:(i+2)*space_steps) == -Gm(i*(space_steps+1)+1:(i+1)*(space_steps+1));
        
        (-dt*Q(i)*K*inv(P(i))*K'-(1/dt)*eye(space_steps)) * dlambda((i+1)*space_steps+1:(i+2)*space_steps) ==...
        dt * Q(i) * K * inv(P(i)) * Gm(i*(space_steps+1)+1:(i+1)*(space_steps+1)) ...
        - Q(i) * drho((i+1)*space_steps+1:(i+2)*space_steps) ...
        - dt*Q(i) * Glambda((i+1)*space_steps+1:(i+2)*space_steps) ...
        - (1/dt) * dlambda(i*space_steps+1:(i+1)*space_steps) ...
        - Grho(i*space_steps+1:(i+1)*space_steps);
        
        K * dm(i*(space_steps+1)+1:(i+1)*(space_steps+1)) + (1/dt) * drho((i+1)*space_steps+1:(i+2)*space_steps) - (1/dt) * drho(i*space_steps+1:(i+1)*space_steps) == -Glambda((i+1)*space_steps+1:(i+2)*space_steps);
        ];
end
eqns = [eqns;
    -(1/dt) * drho(space_steps*time_steps+1:space_steps*(time_steps+1)) == -Glambda(space_steps*(time_steps+1)+1:space_steps*(time_steps+2));
    Q(time_steps) * drho(space_steps*time_steps+1:space_steps*(time_steps+1)) + (1/dt)*dlambda(time_steps*space_steps+1:(time_steps+1)*space_steps) - (1/dt)*dlambda((time_steps+1)*space_steps+1:(time_steps+2)*space_steps) == -Grho(time_steps*space_steps+1:(time_steps+1)*space_steps);
    ];

%% Solve
solution = solve(eqns, [dm(1:end); drho(1:end); dlambda(1:end)]);

%% Helper functions

function Pi = P(i)
    global X
    global rho
    global space_steps
    Pi = 2 * diag(X*(1./(rho(i*space_steps+1:(i+1)*space_steps)) + 1./(rho((i+1)*space_steps+1:(i+2)*space_steps))));
end

function Qi = Q(i)
    global X
    global rho
    global m
    global space_steps
    global time_steps
    if i == 0
        Qi = 2 * diag(X'*(m(i*(space_steps+1)+1:(i+1)*(space_steps+1)).^2)) * diag(rho(i*space_steps+1:(i+1)*space_steps));
    elseif i == time_steps
        Qi = 2 * diag(X'*(m((i-1)*(space_steps+1)+1:i*(space_steps+1)).^2)) * diag(rho(i*space_steps+1:(i+1)*space_steps));
    else
        Qi = 2 * diag(X'*(m(i*(space_steps+1)+1:(i+1)*(space_steps+1)).^2+m((i-1)*(space_steps+1)+1:i*(space_steps+1)).^2)) * diag(rho(i*space_steps+1:(i+1)*space_steps));
    end
end

