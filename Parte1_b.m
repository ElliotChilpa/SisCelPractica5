clc
clear all
close all

a = input("Coordena X de la celda central: ");
b = input("Coordena Y de la celda central: ");
c = input("Valor del lado del hexágono: ");
num_usuarios = input("Número de usuarios por hexágono: ");

figure()
hold on
grid on

% Dibujar hexágonos
colors = ['r', 'b', 'g', 'm', 'y', 'c', 'k'];
apotema = sqrt(3) * c / 2;
L = linspace(0, 2 * pi, 7);

% Matrices para almacenar las coordenadas de los usuarios
randomx_total = [];
randomy_total = [];
hex_centers = [];

for i = 1:7
    if i == 1
        aaux = 0;
        baux = 0;
    else
        aaux = 2 * apotema * cosd(60 * (i - 2) + 30);
        baux = 2 * apotema * sind(60 * (i - 2) + 30);
    end
    
    % Guardar las coordenadas del centro del hexágono
    hex_centers = [hex_centers; a + aaux, b + baux];
    
    % Llamar a la función para dibujar el hexágono y los usuarios
    [x, y, randomx, randomy] = DibujarHexagono_y_usuarios(a + aaux, b + baux, c, num_usuarios);
    
    plot(x, y, 'LineWidth', 2, 'Color', colors(i))
    plot(randomx, randomy, '.', 'Color', colors(i))
    plot(a + aaux, b + baux, 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', colors(i))
    
    % Guardar las coordenadas de los usuarios
    randomx_total = [randomx_total; randomx];
    randomy_total = [randomy_total; randomy];
end

title('Distribución uniforme de usuarios dentro de la celda hexagonal central y adyacentes')
xlabel('X')
ylabel('Y')
axis equal

% Parámetros de la simulación
Potx = 10; % Potencia de Transmisión en dBm
Gtx = 7; % Ganancia de la antena transmisora en dB
Grx = 0; % Ganancia de la antena receptora en dB
f = 2.4e9; % Frecuencia en Hz
c0 = 3e8; % Velocidad de la luz en m/s
lambda = c0 / f; % Longitud de onda en metros
L0 = 20 * log10(4 * pi * 1 / lambda); % Pérdida en la referencia de 1 metro

% Matrices y vectores para almacenar resultados
num_usuarios_total = length(randomx_total); % Número total de usuarios
perdidas = zeros(num_usuarios_total, 7); % Matriz para almacenar las pérdidas
distancias = zeros(num_usuarios_total, 7); % Matriz para almacenar las distancias
omega_total = zeros(num_usuarios_total, 7); % Matriz para almacenar las pérdidas lognormales
celda_asociada = zeros(num_usuarios_total, 1); % Vector para almacenar la celda asociada a cada usuario

% Calcular distancias y pérdidas para todos los usuarios
for i = 1:num_usuarios_total
    min_dist = inf;
    for j = 1:7
        % Calcular la distancia entre el usuario i y el centro del hexágono j
        dist = sqrt((randomx_total(i) - hex_centers(j, 1))^2 + (randomy_total(i) - hex_centers(j, 2))^2);
        distancias(i, j) = dist;
        
        % Generar las pérdidas adicionales con distribución Gaussiana
        sigma = 7.0; % Desviación estándar para las pérdidas adicionales
        
        % Se repite el proceso para sigma = 0 dB, 7 dB y 14 dB
        omega = sigma * randn();
        omega_total(i, j) = omega;
        
        % Calcular las pérdidas usando la fórmula
        perdidas(i, j) = L0 + 10 * 3 * log10(dist) + omega;
        
        % Determinar la celda asociada (la de menor pérdida)
        if perdidas(i, j) < min_dist
            min_dist = perdidas(i, j);
            celda_asociada(i) = j;
        end
    end
end

% Seleccionar arbitrariamente 5 usuarios
selected_users_indices = randperm(num_usuarios_total, 5);

% Mostrar la información para los 5 usuarios seleccionados
for k = 1:5
    i = selected_users_indices(k);
    fprintf('Usuario %d\n', i);
    fprintf('Coordenadas: (%.2f, %.2f)\n', randomx_total(i), randomy_total(i));
    fprintf('Celda asociada: %d\n', celda_asociada(i));
    for j = 1:7
        fprintf('d%d,%d: %.2f, Ω%d,%d: %.2f, L%d,%d: %.2f\n', i, j, distancias(i, j), i, j, omega_total(i, j), i, j, perdidas(i, j));
    end
    fprintf('\n');
end

% Marcar los 5 usuarios seleccionados aleatoriamente en la gráfica
figure(1) % Seleccionar la figura creada anteriormente
hold on
for k = 1:5
    idx = selected_users_indices(k); % Índice del usuario seleccionado
    plot(randomx_total(idx), randomy_total(idx), 'kx', 'MarkerSize', 10, 'LineWidth', 2)
end

% Definición de la función DibujarHexagono_y_usuarios
function [vectores_x, vectores_y, randomx, randomy] = DibujarHexagono_y_usuarios(a, b, c, num_usuarios)
    apotema = sqrt(3) * c / 2;
    L = linspace(0, 2 * pi, 7);

    % Coordenadas del hexágono
    vectores_x = a + c * cos(L);
    vectores_y = b + c * sin(L);

    % Distribuir usuarios aleatoriamente dentro del hexágono
    randomx = [];
    randomy = [];
    while length(randomx) < num_usuarios
        rx_aux = (a - c) + 2 * c * rand(num_usuarios, 1);
        ry_aux = (b - apotema) + 2 * apotema * rand(num_usuarios, 1);
        p = inpolygon(rx_aux, ry_aux, vectores_x, vectores_y);
        randomx = [randomx; rx_aux(p)];
        randomy = [randomy; ry_aux(p)];
        % Limitar el número de usuarios al requerido
        if length(randomx) > num_usuarios
            randomx = randomx(1:num_usuarios);
            randomy = randomy(1:num_usuarios);
        end
    end
end
