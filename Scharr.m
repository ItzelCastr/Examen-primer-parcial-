%Equipo 1
%Integrantes:
%Arroyo Hernández Irving Antonio
%Castro Ramírez Itzel Jenifer
%Cruz San Nicolas Brayan
%Jimenéz Granillo Tania Dyane
%Meraz Ramírez Zoe Jovana
 

%Cargamos las imágenes
Ima_1 = imread('157032.jpg');
Ima_2 = imread('Ruido.jpg');
% Convertimos a escala de grises usando convolución
ImaGris_1 = conv_Gris(Ima_1);
ImaGris_2 = conv_Gris(Ima_2);
function ImaGris = conv_Gris(imagenRGB)
    [filas, columnas, canales] = size(imagenRGB);
    if canales == 1
        ImaGris = imagenRGB; % Si ya es gris, no hace nada
        return;
    end
    % Inicializamos la imagen en escala de grises
    ImaGris = zeros(filas, columnas);
    % Pesos para la conversión
    peso_R = 0.2989; 
    peso_G = 0.5870;
    peso_B = 0.1140;
    % Recorremos cada píxel manualmente
    for i = 1:filas
        for j = 1:columnas
            R = double(imagenRGB(i, j, 1));
            G = double(imagenRGB(i, j, 2));
            B = double(imagenRGB(i, j, 3));
            ImaGris(i, j) = R * peso_R + G * peso_G + B * peso_B;
        end
    end
    % Convertimos a uint8
    ImaGris = uint8(ImaGris);
end



% 1.- Aplicación del Filtro del Promedio a la Imagen 1
% Fórmula del Kernel del Promedio:
% K[s, t] = 1 / (tam_Kernel_Prom * tam_Kernel_Prom)
% Definimos el tamaño del kernel del promedio (5x)
tam_Kernel_Prom = 5;
radioKernelPromedio = floor(tam_Kernel_Prom / 2);
% Conversión de la imagen a double
imagenGrisDouble1 = double(ImaGris_1);
[fil_Ima, col_Ima] = size(imagenGrisDouble1);
% Inicializar la imagen promedio con ceros
Ima_Promedio = zeros(fil_Ima, col_Ima);
% Crear la matriz del kernel del promedio
kernel_promedio_base = [1 1 1 1 1;
                        1 1 1 1 1;
                        1 1 1 1 1;
                        1 1 1 1 1;
                        1 1 1 1 1];
kernel_promedio = kernel_promedio_base / sum(kernel_promedio_base(:));
disp('Matriz del Filtro del Promedio con Kernel:');
disp(kernel_promedio);
% Iterar sobre los píxeles de la imagen
for y = 1:fil_Ima
    for x = 1:col_Ima
        suma_contador = 0;
        contador = 0; % Para contar los píxeles válidos en el vecindario
        % Iterar sobre el vecindario del píxel actual
        for s = -radioKernelPromedio:radioKernelPromedio
            for t = -radioKernelPromedio:radioKernelPromedio
                indiceY = y + s;
                indiceX = x + t;
                % Verificar si los índices están dentro de los límites de la imagen
                if (indiceY >= 1 && indiceY <= fil_Ima && indiceX >= 1 && indiceX <= col_Ima)
                    suma_contador = suma_contador + imagenGrisDouble1(indiceY, indiceX);
                    contador = contador + 1;
                end
            end
        end
        % Calcular el promedio y asignarlo al píxel actual
        Ima_Promedio(y, x) = suma_contador / contador;
    end
end

% 2.- Aplicar el kernel Gaussiano por convolución
% G(x, y) =  1 / (2 * pi * σ²) * exp(-(x² + y²) / (2 * σ²)
% Donde:
% G(x, y) es el valor del kernel Gaussiano en la coordenada (x, y) relativa al centro del kernel.
% σ (sigma) es la desviación estándar de la distribución Gaussiana.
% exp es la función exponencial natural.
% Definir los parámetross
%imagenSuavizada = imgaussfilt(imagenGrisDouble2, sigma);
sigma = 0.8;  % Desviación estándar
kernel_size = 5; % Tamaño del kernel (debe ser impar para un centro definido)
radio_kernel = floor(kernel_size / 2);
% Crear el kernel gaussiano manualmente usando la fórmula
kernel_gaussiano = zeros(kernel_size);
sum_kernel = 0; % Para normalizar
for i = -radio_kernel:radio_kernel
    for j = -radio_kernel:radio_kernel
        kernel_gaussiano(i + radio_kernel + 1, j + radio_kernel + 1) = exp(-(i^2 + j^2) / (2 * sigma^2));
        sum_kernel = sum_kernel + kernel_gaussiano(i + radio_kernel + 1, j + radio_kernel + 1);
    end
end
% Normalizar el kernel
kernel_gaussiano = kernel_gaussiano / sum_kernel;
disp('Matriz del Kernel Gaussiano:');
disp(kernel_gaussiano);
% Obtener las dimensiones de la imagen
[alto_imagen, ancho_imagen] = size(imagenGrisDouble1); % Usar imagenGrisDouble1
% Crear una imagen vacía para almacenar la imagen convolucionada
imagen_Gaussiano = zeros(alto_imagen, ancho_imagen);
% Realizar la convolución manualmente
for i = 1:alto_imagen
    for j = 1:ancho_imagen
        suma_contador = 0;
        for m = -radio_kernel:radio_kernel
            for n = -radio_kernel:radio_kernel
                indice_y = i + m;
                indice_x = j + n;
                % Verificar límites de la imagen
                if (indice_y >= 1 && indice_y <= alto_imagen && indice_x >= 1 && indice_x <= ancho_imagen)
                    peso = kernel_gaussiano(m + radio_kernel + 1, n + radio_kernel + 1);
                    valor_pixel_vecino = imagenGrisDouble1(indice_y, indice_x);
                    suma_contador = suma_contador + peso * valor_pixel_vecino;
                end
            end
        end
        imagen_Gaussiano(i, j) = suma_contador;
    end
end


% Mostramos las imágenes de la figura 1
figure;imshow(ImaGris_1); title('Entrada de Imagen 1');
figure; imshow(uint8(Ima_Promedio)); title('Imagen 1 con Filtro de Promedio');
figure; imshow(uint8(imagen_Gaussiano)); title('Imagen 1 con Kernel Gaussiano por convolución');


% 3.- Aplicación del Filtro de la Mediana a la Imagen 2
% Definir el tamaño del vecindario para la mediana (debe ser impar)
tamanoVecindarioMediana = 3;
radioVecindarioMediana = floor(tamanoVecindarioMediana / 2);
% Convertir la imagen a double
imagenGrisDouble2 = double(ImaGris_2);
[fil_Ima2, col_Ima2] = size(imagenGrisDouble2);
% Inicializar la imagen resultante con ceros
imagenMediana2 = zeros(fil_Ima2, col_Ima2);
% Iterar sobre los píxeles de la imagen
for y = 1:fil_Ima2
    for x = 1:col_Ima2
        vecindario = []; % Inicializar un array para almacenar los valores del vecindario
        % Iterar sobre el vecindario del píxel actual
        for s = -radioVecindarioMediana:radioVecindarioMediana
            for t = -radioVecindarioMediana:radioVecindarioMediana
                indiceY = y + s;
                indiceX = x + t;
                % Verificar si los índices están dentro de los límites de la imagen
                if (indiceY >= 1 && indiceY <= fil_Ima2 && indiceX >= 1 && indiceX <= col_Ima2)
                    vecindario = [vecindario, imagenGrisDouble2(indiceY, indiceX)];
                end
            end
        end
        % Calcular la mediana del vecindario
        vecindarioOrdenado = sort(vecindario);
        numeroElementos = length(vecindarioOrdenado);
        if rem(numeroElementos, 2) == 1
            % Número de elementos impar, la mediana es el elemento central
            medianaVecindario = vecindarioOrdenado((numeroElementos + 1) / 2);
        else
            % Número de elementos par, la mediana es el promedio de los dos elementos centrales
            medio1 = vecindarioOrdenado(numeroElementos / 2);
            medio2 = vecindarioOrdenado((numeroElementos / 2) + 1);
            medianaVecindario = (medio1 + medio2) / 2;
        end
        % Asignar la mediana al píxel actual
        imagenMediana2(y, x) = medianaVecindario;
    end
end


% 4. Derivada respecto a X (Con Kernel Scharr):
[filas, columnas] = size(imagenMediana2); % Obtiene dimensiones de la imagen suavizada.
derivadaX_Scharr = zeros(filas, columnas); % Inicializa matriz para derivada en X.
kernel_dx_scharr = [3 0 -3; 10 0 -10; 3 0 -3]; % Kernel de Scharr para derivada horizontal.
disp('Matriz del Kernel Scharr para la Derivada en X:'); % Muestra el kernel.
disp(kernel_dx_scharr);
radio_kernel = floor(size(kernel_dx_scharr, 1) / 2); % Radio del kernel.
for y = 1:filas % Itera sobre filas.
    for x = 1:columnas % Itera sobre columnas.
        suma_contador_x = 0; % Inicializa suma para la derivada en X.
        for ky = -radio_kernel:radio_kernel % Itera sobre filas del kernel.
            for kx = -radio_kernel:radio_kernel % Itera sobre columnas del kernel.
                indiceY = y + ky; % Índice de fila del vecino.
                indiceX = x + kx; % Índice de columna del vecino.
                % Verificar límites de la imagen
                if (indiceY >= 1 && indiceY <= filas && indiceX >= 1 && indiceX <= columnas) % Si el vecino está dentro de la imagen.
                    peso_x = kernel_dx_scharr(ky + radio_kernel + 1, kx + radio_kernel + 1); % Peso del kernel.
                    valor_pixel = imagenMediana2(indiceY, indiceX); % Valor del píxel vecino.
                    suma_contador_x = suma_contador_x + peso_x * valor_pixel; % Acumula la suma ponderada.
                end
            end
        end
        derivadaX_Scharr(y, x) = suma_contador_x; % Asigna el valor de la derivada en X.
    end
end

% 5. Derivada respecto a Y (Con Kernel Scharr):
derivadaY_Scharr = zeros(filas, columnas); % Inicializa matriz para derivada en Y.
kernel_dy_scharr = [3 10 3; 0 0 0; -3 -10 -3]; % Kernel de Scharr para derivada vertical.
disp('Matriz del Kernel Scharr para la Derivada en Y:'); % Muestra el kernel.
disp(kernel_dy_scharr);
radio_kernel = floor(size(kernel_dy_scharr, 1) / 2); % Radio del kernel.
for y = 1:filas % Itera sobre filas.
    for x = 1:columnas % Itera sobre columnas.
        suma_contador_y = 0; % Inicializa suma para derivada en Y.
        for ky = -radio_kernel:radio_kernel % Itera sobre filas del kernel.
            for kx = -radio_kernel:radio_kernel % Itera sobre columnas del kernel.
                indiceY = y + ky; % Índice de fila del vecino.
                indiceX = x + kx; % Índice de columna del vecino.
                % Verificar límites de la imagen
                if (indiceY >= 1 && indiceY <= filas && indiceX >= 1 && indiceX <= columnas) % Si el vecino está dentro de la imagen.
                    peso_y = kernel_dy_scharr(ky + radio_kernel + 1, kx + radio_kernel + 1); % Peso del kernel.
                    valor_pixel = imagenMediana2(indiceY, indiceX); % Valor del píxel vecino.
                    suma_contador_y = suma_contador_y + peso_y * valor_pixel; % Acumula la suma ponderada.
                end
            end
        end
        derivadaY_Scharr(y, x) = suma_contador_y; % Asigna el valor de la derivada en Y.
    end
end

% 6. Magnitud del Gradiente (usando las derivadas de Scharr):
magnitudGradiente_Scharr = zeros(filas, columnas); % Inicializa matriz para magnitud del gradiente.
for y = 1:filas % Itera sobre filas.
    for x = 1:columnas % Itera sobre columnas.
        magnitudGradiente_Scharr(y, x) = sqrt(derivadaX_Scharr(y, x)^2 + derivadaY_Scharr(y, x)^2); % Calcula la magnitud del gradiente.
    end
end

% Aplicamos valor absoluto a las derivadas de X e Y para visualización
derivadaX_Scharr_abs = abs(derivadaX_Scharr);
derivadaY_Scharr_abs = abs(derivadaY_Scharr);
% Normalizar las derivadas para mejor visualización
derivadaX_Scharr_norm = (derivadaX_Scharr_abs - min(derivadaX_Scharr_abs(:)))/ (max(derivadaX_Scharr_abs(:)) - min(derivadaX_Scharr_abs(:))); % Normaliza a [0,1]

derivadaY_Scharr_norm = (derivadaY_Scharr_abs - min(derivadaY_Scharr_abs(:))) / (max(derivadaY_Scharr_abs(:)) - min(derivadaY_Scharr_abs(:))); % Normaliza a [0,1]

magnitudGradiente_Scharr_norm = (magnitudGradiente_Scharr - min(magnitudGradiente_Scharr(:))) / (max(magnitudGradiente_Scharr(:)) - min(magnitudGradiente_Scharr(:))); % Normaliza a [0,1]



% Mostrar las imágenes de la figura 2
figure;imshow(ImaGris_2); title('Entrada de Imagen 2');
figure;imshow(uint8(imagenMediana2));title('Imagen 2 con Filtro de la Mediana');

% Mostrar las imágenes
figure;imshow(derivadaX_Scharr_norm, []); title('Derivada respecto a X (Scharr)');
figure;imshow(derivadaY_Scharr_norm, []); title('Derivada respecto a Y (Scharr)');
figure; imshow(magnitudGradiente_Scharr_norm, []); title('Magnitud del Gradiente (Scharr)');