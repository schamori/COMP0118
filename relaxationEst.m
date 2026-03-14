function [a, r, g, f] = relaxationEst(y, TE, Nrow, Ncol, lambdaA, lambdaR)

AL_iter = 50;
mu = 1e-1;

Npix = size(y, 2);
t = 1;
res = inf;
tol = 1e-3;
tol1 = sqrt(Npix)*tol;
tol2 = sqrt(Npix)*tol;
res_p = inf;
res_d = inf;

a = max(y(:)) * rand(Npix, 1);
r = 0.1*ones(Npix, 1);

g = a;
f = r;

d1 = zeros(size(g));
d2 = zeros(size(f));

puxg = 0.*reshape(g, Nrow, Ncol);
puyg = 0.*reshape(g, Nrow, Ncol);

puxr = 0.*reshape(r, Nrow, Ncol);
puyr = 0.*reshape(r, Nrow, Ncol);

disp(' ------- Running T2* Relaxation Estimation ------- ')

while (t<AL_iter) && ((abs (res_p) > tol1) || (abs (res_d) > tol2))
    g01 = g;
    f02 = f;
    
    %% Solve for g
    tmp = reshape(a - d1, Nrow, Ncol);
    [g_img, puxg, puyg] = chambolle_prox_TV_stop(tmp, lambdaA/mu, puxg, puyg);
    g = g_img(:);
    
    %% Solve for f
    tmp = reshape(r - d2, Nrow, Ncol);
    [f_img, puxr, puyr] = chambolle_prox_TV_stop(tmp, lambdaR/mu, puxr, puyr);
    f = f_img(:);
    f_d2 = f + d2;

    for oo = 1:Npix
        
        %% Solve for a       
        E = exp(-r(oo)*TE(:));
        a(oo) = (E' * y(:,oo) + mu*(g(oo) - d1(oo))) / (E' * E + mu);

        %% Solve for r
        [r(oo), niter, rnorm, dx] = gradientDescentAK(r(oo), TE, a(oo), y(:, oo), mu, f_d2(oo));
    end
   
    %% Update Lagrange multipliers
    res1 = a - g;
    d1 = d1 - res1;
    
    res2 = r - f;
    d2 = d2 - res2;
    
    %% Update mu
    if mod(t, 10) == 1
        res_p = sqrt(norm(res1, 'fro')^2 + norm(res2, 'fro')^2);    % primal residue
        res_d = mu*(norm(g01 - g + f02 - f,'fro'));                 % dual residue
        
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d1 = d1/2;
            d2 = d2/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            d1 = d1*2;
            d2 = d2*2;
        end
    end
    res(1) = norm(res1, 'fro');
    res(2) = norm(res2, 'fro');
    
    if mod(t, 10) == 1
        fprintf(' t = %f, res_p = %f, res_d = %f, mu = %f, niterGD = %f \n', t, res_p, res_d, mu, niter)
    end
    
    t = t + 1;
    
end
disp(' ------- T2* Relaxation Estimation Completed ------- ')

end

%