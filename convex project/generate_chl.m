function H = generate_chl(Lp, M, K, d_lambda)
% generate_chl(alpha, phi)
% Generates a mmwave channel with Lp number of path, M bs antennas, and K
% single antenna users. d_lambda is the ratio between the antenna elements separation and the wavelength.
c = 2 * pi * d_lambda;

H = zeros(M, K);
for kk = 1:K
    alpha = 1/sqrt(2) * (randn(Lp, 1) + 1j * randn(Lp, 1));
    phi = 2 * pi * rand(1, Lp);
    [X, Y] = meshgrid(phi, 0:(M-1));
    A = exp(1j  * c * sin(X) .* Y);
    H(:, kk) = 1/sqrt(Lp) * A * alpha;
end

end