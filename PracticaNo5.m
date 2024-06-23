clc;
clear all;
close all;

% Solicitar inputs al usuario
%a = input("Coordena x del centro de la celda principal: ");
% b = input("Coordena y del centro de la celda principal: ");
a = 0;
b = 0;
% c = input("Ingresa el valor del lado del hexagono (km): ");
c = 3;
d = input("Ingresa la cantidad de usuarios por hexágono: ");

% Parámetros del sistema
alpha = 10;
apotema = sqrt(3) * c / 2;
desvia = 7;

%ad = makedist('Normal', 'mu', 0, 'sigma', desvia);

% Dibujar el hexágono central y los hexágonos alrededor
[x, y, rx, ry, centers] = DibujarHexagonos_y_usuarios(a, b, c, d);

% Graficar los hexágonos y los usuarios
figure();
hold on;
grid on;
for i = 1:length(centers)
    plot(x{i}, y{i}, 'r', 'LineWidth', 2); % Graficar cada hexágono
    plot(rx{i}, ry{i}, 'x'); % Graficar los usuarios en cada hexágono
    plot(centers(i, 1), centers(i, 2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'b'); % Marcar el centro de cada hexágono
end
title('Usuarios dentro de los hexágonos');
xlabel('Eje X');
ylabel('Eje Y');
legend('Hexágonos', 'Usuarios', 'Centros');

% Función para dibujar hexágonos y generar usuarios
function [vectores_x, vectores_y, randomx, randomy, centers] = DibujarHexagonos_y_usuarios(a, b, c, d)
    apotema = sqrt(3) * c / 2;
    L = linspace(0, 2 * pi, 7);

    % Centros de los hexágonos
    centers = [
        a, b;
        a + 2 * apotema * cosd(30), b + 2 * apotema * sind(30);
        a + 2 * apotema * cosd(90), b + 2 * apotema * sind(90);
        a + 2 * apotema * cosd(150), b + 2 * apotema * sind(150);
        a + 2 * apotema * cosd(210), b + 2 * apotema * sind(210);
        a + 2 * apotema * cosd(270), b + 2 * apotema * sind(270);
        a + 2 * apotema * cosd(330), b + 2 * apotema * sind(330)
    ];

    vectores_x = cell(1, 7);
    vectores_y = cell(1, 7);
    randomx = cell(1, 7);
    randomy = cell(1, 7);

    for i = 1:7
        a_center = centers(i, 1);
        b_center = centers(i, 2);
        vectores_x{i} = a_center + c * cos(L);
        vectores_y{i} = b_center + c * sin(L);

        % Generar d usuarios aleatorios dentro del hexágono
        rx_aux = (a_center - c) + (2 * c) * rand(d, 1);
        ry_aux = (b_center - apotema) + (2 * apotema) * rand(d, 1);
        p = inpolygon(rx_aux, ry_aux, vectores_x{i}, vectores_y{i});
        rx_aux = rx_aux(p);
        ry_aux = ry_aux(p);

        randomx{i} = rx_aux;
        randomy{i} = ry_aux;

        while length(randomx{i}) < d
            p = false;
            while ~p
                rx_aux = (a_center - c) + (2 * c) * rand(1, 1);
                ry_aux = (b_center - apotema) + (2 * apotema) * rand(1, 1);
                p = inpolygon(rx_aux, ry_aux, vectores_x{i}, vectores_y{i});
            end
            randomx{i}(end + 1) = rx_aux;
            randomy{i}(end + 1) = ry_aux;
        end
    end
end