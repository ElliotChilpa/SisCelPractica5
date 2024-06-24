clc;
clear all;
close all;

% Definir coordenadas del centro de la celda principal
a = 0;
b = 0;

% Definir el lado del hexágono en km
c = 3;

% Solicitar al usuario la cantidad de usuarios por hexágono
d = input("Ingresa la cantidad de usuarios por hexágono: ");

% Exponente de pérdidas por distancia
alpha = 10;

% Desviación estándar de pérdidas por ensombrecimiento
desvia = 7;
ad = makedist('Normal', 'mu', 0, 'sigma', desvia);

% Calculo del apotema del hexágono
apotema = (sqrt(3) * c) / 2;

% Función para dibujar hexágonos y generar usuarios
[x, y, rx, ry, centers] = DibujarHexagonos_y_usuarios(a, b, c, d);

% Colores para los hexágonos
colors = lines(7); % Utiliza la paleta de colores 'lines' para tener colores distintos

% Graficar los hexágonos y los usuarios
figure();
hold on;
grid on;
for i = 1:length(centers)
    plot(x{i}, y{i}, 'Color', colors(i, :), 'LineWidth', 2); % Graficar cada hexágono con un color diferente
    plot(rx{i}, ry{i}, 'x', 'Color', colors(i, :)); % Graficar los usuarios en cada hexágono con el mismo color
    plot(centers(i, 1), centers(i, 2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', colors(i, :)); % Marcar el centro de cada hexágono
end
title('Usuarios dentro de los hexágonos');
xlabel('Eje X');
ylabel('Eje Y');
legend('Hexágonos', 'Usuarios', 'Centros');

% Calcular y mostrar pérdidas de señal
figure();
hold on;
grid on;
for i = 1:length(centers)
    for j = 1:length(rx{i})
        d_ij = sqrt((rx{i}(j) - centers(i, 1))^2 + (ry{i}(j) - centers(i, 2))^2);
        Omega_ij = random(ad);
        L_ij = 10 * alpha * log10(d_ij) + Omega_ij;
        plot(d_ij, L_ij, 'o', 'Color', colors(i, :));
    end
end
title('Pérdidas de señal según el modelo lognormal');
xlabel('Distancia (km)');
ylabel('Pérdida de señal (dB)');
legend('Pérdidas de señal por hexágono');

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

    % Inicializar variables de salida
    vectores_x = cell(1, 7);
    vectores_y = cell(1, 7);
    randomx = cell(1, 7);
    randomy = cell(1, 7);

    % Generar hexágonos y usuarios aleatorios
    for i = 1:7
        a_center = centers(i, 1);
        b_center = centers(i, 2);
        vectores_x{i} = a_center + c * cos(L);
        vectores_y{i} = b_center + c * sin(L);

        % Generar d usuarios aleatorios dentro del hexágono
        [rx, ry] = GenerarUsuariosHexagono(a_center, b_center, c, apotema, d, vectores_x{i}, vectores_y{i});
        randomx{i} = rx;
        randomy{i} = ry;
    end
end

% Función para generar usuarios aleatorios dentro de un hexágono
function [rx, ry] = GenerarUsuariosHexagono(a_center, b_center, c, apotema, d, vectores_x, vectores_y)
    rx_aux = (a_center - c) + (2 * c) * rand(d, 1);
    ry_aux = (b_center - apotema) + (2 * apotema) * rand(d, 1);
    p = inpolygon(rx_aux, ry_aux, vectores_x, vectores_y);
    rx = rx_aux(p);
    ry = ry_aux(p);

    while length(rx) < d
        p = false;
        while ~p
            rx_aux = (a_center - c) + (2 * c) * rand(1, 1);
            ry_aux = (b_center - apotema) + (2 * apotema) * rand(1, 1);
            p = inpolygon(rx_aux, ry_aux, vectores_x, vectores_y);
        end
        rx(end + 1) = rx_aux;
        ry(end + 1) = ry_aux;
    end
end
