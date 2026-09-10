%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% PhC_1D.m                                                         x      
% Matriz de Transferência para Cristal Fotônico 1D                 x
%     Nº par de camadas A = Ti2O5 & B = SiO2                       x
%     Índice de Ref. na região visível: n1 ~ 2.1 e n2 ~ 1.46       x
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

format long 
n1 = 2.14;                  % IR, camada A (Ti2O5)
n2 = 1.46;                  % IR, camada B (SiO2)
n0 = 1.0;                   % O típico é ar, mas pode mudar
n_out = 1.0;                % Pode ser o ar ou outro substrato
N = 6;                      % Nº de pares de camadas
lamb0 = 532.0;              % Comprimento de onda de Bragg em nm
d1 = lamb0/(4*n1);
d2 = lamb0/(4*n2);        

% Matriz de Transferência de Interface (MTI) para a interface n0 -|- n1
n01_p = (n0 + n1)/(2*n0); n01_m = (n0 - n1)/(2*n0);
B01 = [n01_p, n01_m; n01_m, n01_p]; 

% MTI para a interface n1 -|- n2
n12_p = (n1 + n2)/(2*n1); n12_m = (n1 - n2)/(2*n1);
B12 = [n12_p, n12_m; n12_m, n12_p]; 

% MTI para a interface n2 -|- n1 
n21_p = (n2 + n1)/(2*n2); n21_m = (n2 - n1)/(2*n2);
B21 = [n21_p, n21_m; n21_m, n21_p]; 

% MTI de saída para a interface n1 -|- n_out
n1t_p = (n1 + n_out)/(2*n1); n1t_m = (n1 - n_out)/(2*n1);
B1t = [n1t_p, n1t_m; n1t_m, n1t_p]; 

% MTI de saída para a interface n2 -|- n_out
n2t_p = (n2 + n_out)/(2*n2); n2t_m = (n2 - n_out)/(2*n2);
B2t = [n2t_p, n2t_m; n2t_m, n2t_p]; 

rt = fopen('rt0.dat', 'w+');          % para salvar os dados no arquivo rt0.dat
wavelength = 400:1:700;               % intervalo do comprimento de onda (nm)
R0 = zeros(1,length(wavelength));     % para armazenar os valores de reflexividade
T0 = zeros(1,length(wavelength));     % para armazenar os valores de transmissividade

for k = 1:length(wavelength)
    lambda = wavelength(k); 

% Matriz de propagação na camada de Ti2O5
    alfa1 = exp(1j*2*pi*n1*d1/lambda);
    A1 = [alfa1, 0; 0, 1/(alfa1)]; 

% Matriz de propagação na camada de SiO2 
    alfa2 = exp(1j*2*pi*n2*d2/lambda);
    A2 = [alfa2, 0; 0, 1/(alfa2)]; 


    %============================================================
    % Espelho M0                                                x
    %   n0-|-n1-|-|-n2-|...|-n1-|-|-n2-|-|-n1-|-|-n2-|-n_out    x
    %   B01-(A1-B12-A2-B21)...(A1-B12-A2-B21)-A1-B12-A2-B2t     x
    %   B01-(A1-B12-A2-B21)^(N-1)-A1-B12-A2-B2t                 x
    %============================================================

    % Método da Matriz de Transferência (MMT) para o Cristal Fotônico 1D
    M0 = B01*(A1*B12*A2*B21)^(N-1)*A1*B12*A2*B2t; 

    % Espectros de reflexão & transmissão (R0 & T0) do Cristal Fotônico
    R0(k) = abs(M0(2,1)/M0(1,1))^2;
    T0(k) = n_out/n0*abs(1/M0(1,1))^2; 

    % Salvar dados em três colunas para lambda, R0, T0 
    fprintf(rt, '%f %12.8f %12.8f\n', lambda, R0(k), T0(k));
end

tick1 = 400:50:700;
figure(1)
  plot(wavelength, R0, 'r', wavelength, T0, 'b', 'linewidth', 2)
  legend('R_0', 'T_0'); grid minor
  set(gca, 'FontSize', 18);
  axis([400 700 0 1.05]);
  set(gca, 'XTick', tick1);
  set(gca, 'YTick', 0:0.2:1);
  xlabel('Wavelength (nm)');
  ylabel('R & T');
  title('A = Ti2O5, B = SiO2, N = 6');

fclose(rt);
%==========================================================================