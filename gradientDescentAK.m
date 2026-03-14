function [ropt, niter, gnorm, dr] = gradientDescentAK(r0, TE, a, y, mu, f_d2)

maxiter = 1000;         % Maximum number of allowed iterations
drmin = 1e-6;           % Minimum allowed perturbation
alpha = 1e-8;           % Step size (can be sensitive for different T2*)
tol = 1e-6;             % Termination tolerance

gnorm = inf;            % Onitialize gradient norm
r = r0;                 % Initilize result
niter = 0;              % Iteration counter
dr = inf;               % Perturbation

while and(gnorm >= tol, and(niter <= maxiter, dr >= drmin))
    gradientR = grad(TE, r, a, y, mu, f_d2);
    gnorm = norm(gradientR);
    rnew = r - alpha * gradientR;
    rnew = max(rnew, 1e-6);     % enforce positivity: r = 1/T2* must be > 0
    niter = niter + 1;
    dr = norm(rnew - r);
    r = rnew;
end
ropt = r;
niter = niter - 1;

end

function gradientR = grad(TE, r, a, y, mu, f_d2)

E = exp(-r .* TE);              % exp(-r t)
residual = y - a .* E;          % (s - a exp(-rt))
sensitivity = TE .* E;          % t ⊙ exp(-rt)
data_grad = a * (sensitivity' * residual);   % a (t⊙exp(-rt))^T (s - a exp(-rt))
penalty_grad = mu * (r - f_d2); % μ(r - f - d2)
gradientR = data_grad + penalty_grad;

end
