% EM algorithm
function [Mn, alpha_new] = EMupdates(Phi, A, t, sigma2)
    % E step
    Sn = inv(sigma2.^-1*(Phi'*Phi) + A);
    Mn = sigma2.^-1* Sn * Phi' *t;

    % M step
    % gamma = 1 - alpha .* diag(Sn);
    alpha_new = 1 ./ (Mn.^2 + diag(Sn));
end

