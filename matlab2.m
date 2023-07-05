clear all
close all
clc

A = 1;
z0 = [0, 0, 1];
mu1 = 4 * pi * 1e-7; %магнитная проницаемость среды 
mu2 = 4 * pi * 1e-7; %магнитная проницаемость цилиндра
n = 1.5; %показатель преломления цилиндра
lambda = 633e-9; %длина волны
k1 = 2 * pi / lambda; %волновое число для среды
k2 = k1 * n; %волновое число для цилиндра
m = -15 : 1 : 15; % порядок уравнения в виде скаляра
Rs = [lambda / 2, lambda, lambda * 2]; %радиусы цилиндра

for rs = 1 : 3
    R = Rs(1, rs);
    %Коэффициенты:
    b = ((k1 * mu2 * besselh(m, 2, k1 * R) .* (besselj(m - 1,k1 * R) - (m) / (k1 * R) .* besselj(m, k1 * R))) ...
    - (k1 * mu2 * (besselh(m - 1, 2, k1 * R) - (m) / (k1 * R) .* besselh(m, 2, k1 * R)) .* besselj(m, k1 * R))) ./ ...
    ((k2 * mu1 * (besselj(m - 1, k2 * R) - (m) / (k2 * R) .* besselj(m, k2 * R)) .* besselh(m, 2, k1 * R)) - ...
    (k1 * mu2 * (besselh(m - 1, 2, k1 * R) - (m) / (k1 * R) .* besselh(m, 2, k1 * R)) .* besselj(m, k2 * R)));

    c = ((k2 * mu1 * (besselj(m - 1, k2 * R) - (m) / (k2 * R) .* besselj(m, k2 * R)) .* besselj(m, k1 * R)) ...
    - (k1 * mu2 * besselj(m, k2 * R) .* (besselj(m - 1, k1 * R) - (m) / (k1 * R) .* besselj(m, k1 * R)))) ./ ...
    ((k1 * mu2 * besselj(m, k2 * R) .* (besselh(m - 1, 2, k1 * R) - (m) / (k1 * R) .* besselh(m, 2, k1 * R))) - ...
    (k2 * mu1 * (besselj(m - 1, k2 * R) - (m) / (k2 * R) .* besselj(m, k2 * R)) .* besselh(m, 2, k1 * R)));

    size = 200; %размерность осей
    x = linspace(-5 * R, 5 * R, size);
    y = linspace(-5 * R, 5 * R, size);
    %Матрицы интенсивностей:
    Plus_matrix= zeros(size, size);
    Minus_matrix = zeros(size, size);

    for i = 1 : size
        for j = 1 : size
            % Преобраз точек декартовой системы в точки полярной системы:
            [fi, r] = cart2pol(x(i),y(j));
            %Распределение интенсивности внутреннего и внешнего поля дифракции:
            E_plus = z0 * A * sum(b .* (-1i) .^ m .* besselj(m, k2 * r) .* exp(1i * m * fi)); %внутреннее
            E_minus = z0 * A * sum(c .* (-1i) .^ m .* besselh(m, 2, k1 * r) .* exp(1i * m * fi)); %внешнее
            if r >= R
                Minus_matrix(i, j) = norm(E_minus); 
            else
                Plus_matrix(i, j) = norm(E_plus);
            end
        end
    end

    subplot(2, 3, rs);
    imagesc(x, y, rot90(Minus_matrix));
    title('E-', R);
    xlabel('X, м');
    ylabel('Y, м');
    colorbar;
    axis xy;

    subplot(2, 3, rs + 3);
    imagesc(x, y, rot90(Plus_matrix));
    title('E+', R);
    xlabel('X, м');
    ylabel('Y, м');
    colorbar;
    axis xy;

end
