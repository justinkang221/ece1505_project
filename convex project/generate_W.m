function W = generate_W(params)
    Qb = params.Qb;
    Nrf = params.Nrf;
    M = params.M;
    L = params.L;
    c = randi(2^Qb - 1,L * Nrf, 1);
    [X, Y] = meshgrid(0:(M-1), c);
    Wrf = exp(1j * 2 * pi/2^Qb * X .* Y);
    % Switching Matrix (a random permutation matix)
    idx = randperm(M);
    A = zeros(L * Nrf, M);
    for ii = 1:(L * Nrf)
        A(ii, idx(ii)) = 1;
    end
    W = zeros(L * Nrf, L * M);
    for ii = 1:L
        temp = A((ii-1)* Nrf + (1:Nrf), :).*Wrf((ii-1)* Nrf + (1:Nrf), :);
        W((ii-1)* Nrf + (1:Nrf), (ii-1)* M + (1:M)) = temp;
    end
end