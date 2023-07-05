%5
%clc %очистка экрана консоли
clear all %очистка значений переменных
close all %закрытие всех открытых файлов



for NA = [0.65, 0.8, 0.95] 
    %расчет максимального угла, определяемого числовой аппертурой
    alpha = asin(NA); % NA/n; n = 1
   
    %анонимные функции для расчета значения l(theta), формулы 5.10 5.11
    %l(theta)-начальное распределение напряженности электрического поля в координатах выходного зрачка
    ltan = @(theta) sin(theta).*exp(-(sin(theta)/sin(alpha)).^2); %5.10
    lsin = @(theta) tan(theta).*exp(-(tan(theta)/tan(alpha)).^2); %5.11
   
    %анонимные функции для рачета значения Т1,2(theta)
    %T1(theta) - функция аподизации для апланатического объектива
    T1 = @(theta) cos(theta).^0.5;
    %T2(theta) - функция аподизации для плоской линзы
    T2 = @(theta) cos(theta).^-1.5;
    %поперечное сечение
    r = linspace(-1, 1, 100);
    %фокальная плоскость располагается в 0
    z = 0;

    %/////////Для апланатического объектива/////////////
   
    %анонимная функция для расчета по формуле (5.8)
    calculate_Er = @(theta, r, z) lsin(theta) .* T1(theta) .* sin(2 * theta) .*...
    exp(1i * (2 * pi / 0.532) .* z .* cos(theta)) .* besselj(1, (2 * pi / 0.532) .* r .* sin(theta));
    %анонимная функция для расчета по формуле (5.9)
    calculate_Ez = @(theta, r, z) lsin(theta).*T1(theta).*(sin(theta).^2) .*...
    exp(1i * (2 * pi / 0.532) .* z .* cos(theta)) .* besselj(0, (2 * pi / 0.532) .* r .* sin(theta));
   
    %расчет интегралов с помощью анонимных функций
    %уравнения Ричардса-Вольфа для радиально поляризованного света
    Er = integral(@(theta) calculate_Er(theta, r, z), 0, alpha, 'ArrayValued', true);
    Ez = 2i * integral(@(theta) calculate_Ez(theta, r, z), 0, alpha, 'ArrayValued', true);
    %по полученным значениям состовляющих считается интенсивность
    I = abs(Er).^2 + abs(Ez).^2;
   
    %вывод графиков
    figure
    subplot(2,1,1)
    plot(r, I);
    title(strcat('Апланатический объектив NA = ',num2str(NA)));
    xlabel('y, мкм');
    ylabel('Intensity');
    grid on;

   
    %/////////Для дифракционной линзы/////////////
       
    %анонимная функция для расчета по формуле (5.8)
    calculate_Er = @(theta, r, z) ltan(theta) .* T2(theta) .* sin(2 * theta) .*...
    exp(1i * (2 * pi / 0.532) .* z .* cos(theta)) .* besselj(1, (2 * pi / 0.532) .* r .* sin(theta));
    %анонимная функция для расчета по формуле (5.9)
    calculate_Ez = @(theta, r, z) ltan(theta) .* T2(theta) .* (sin(theta).^2) .*...
    exp(1i * (2 * pi / 0.532) .* z .* cos(theta)) .* besselj(0, (2 * pi / 0.532) .* r .* sin(theta));
   
    %расчет интегралов с помощью анонимных функций
    %уравнения Ричардса-Вольфа для радиально поляризованного света
    Er = integral(@(theta) calculate_Er(theta, r, z), 0, alpha,'ArrayValued',true);
    Ez = 2i * integral(@(theta) calculate_Ez(theta, r, z), 0, alpha, 'ArrayValued',true);
    %по полученным значениям состовляющих считается интенсивность
    I = abs(Er).^2 + abs(Ez).^2;
   
    %вывод графиков
    subplot(2,1,2)
    plot(r, I);
    title(strcat('Дифракционная линза NA = ',num2str(NA)));
    xlabel('y, мкм');
    ylabel('Intensity');
    grid on;


end
