%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% LinearCavity.m                                                   x      
% Matriz de Transferência para Cristal Fotônico defeituoso         x
%                                                                  x
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

format long 
n1 = 2.14;                  % IR, camada A (Ti2O5)
n2 = 1.46;                  % IR, camada B (SiO2)
nc = 1.60;                  % IR, camada C (defeito) 
n_in = 1.0;                   
n_out = 1.0;               

lamb0 = 532.0;             % Comprimento de onda de Bragg (nm)
d1 = lamb0/(4*n1);
d2 = lamb0/(4*n2);
dc = 20000;

N1 = 8;
N2 = 6;

% Matriz de Transferência de Interface (MTI) para a interface n0 -|- n1
n01_p = (n_in + n1)/(2*n_in); n01_m = (n_in - n1)/(2*n_in);
B01 = [n01_p, n01_m; n01_m, n01_p]; 

% MTI para a interface n1 -|- n2
n12_p = (n1 + n2)/(2*n1); n12_m = (n1 - n2)/(2*n1);
B12 = [n12_p, n12_m; n12_m, n12_p]; 

% MTI para a interface n2 -|- n1 
n21_p = (n2 + n1)/(2*n2); n21_m = (n2 - n1)/(2*n2);
B21 = [n21_p, n21_m; n21_m, n21_p]; 

% MTI para a interface n2 -|- nc
n2c_p = (n2 + nc)/(2*n2); n2c_m = (n2 - nc)/(2*n2);
B2c   = [n2c_p, n2c_m; n2c_m, n2c_p];

% MTI para a interface nc -|- n2
nc2_p = (nc + n2)/(2*nc); nc2_m = (nc - n2)/(2*nc);
Bc2 = [nc2_p, nc2_m; nc2_m, nc2_p];

% MTI de saída para a interface n1 -|- nL
n1t_p = (n1 + n_out)/(2*n1); n1t_m = (n1 - n_out)/(2*n1);
B1t = [n1t_p, n1t_m; n1t_m, n1t_p]; 

% MTI de saída para a interface n2 -|- nL
n2t_p = (n2 + n_out)/(2*n2); n2L_m = (n2 - n_out)/(2*n2);
B2t = [n2t_p, n2L_m; n2L_m, n2t_p]; 

%==============================================================
% Modelo                                                      x
%                Grade M1               Grade M2              x
%        --->-|A|B|A|B|...|A|B|====|B|A|B|A|...|B|A|-->-      x
%                 pares N1      dc       pares N2             x 
%==============================================================

wavelength = 450:0.002:650;               % intervalo do comprimento de onda (nm)
R01 = zeros(1,length(wavelength));         
R02 = zeros(1,length(wavelength));  
R = zeros(1, length(wavelength));
T = zeros(1, length(wavelength)); 

for k = 1:length(wavelength)
    lambda = wavelength(k); 

% Matriz de propagação na camada de Ti2O5
    alfa1 = exp(1j*2*pi*n1*d1/lambda);
    A1 = [alfa1, 0; 0, 1/(alfa1)]; 

% Matriz de propagação na camada de SiO2 
    alfa2 = exp(1j*2*pi*n2*d2/lambda);
    A2 = [alfa2, 0; 0, 1/(alfa2)];

    alfaC = exp(1j*2*pi*nc*dc/lambda);
       Ac = [alfaC, 0; 0, 1/(alfaC)]; 

    %============================================================
    % Espelho M1                                                x
    %   n0-|-n1-|-|-n2-|...|-n1-|-|-n2-|-|-n1-|-|-n2-|=nc       x
    %   B01-(A1-B12-A2-B21)...(A1-B12-A2-B21)-A1-B12-A2-B2c     x
    %   B01-(A1-B12-A2-B21)^(N-1)-A1-B12-A2-B2c                 x
    %============================================================
    M1 = B01*(A1*B12*A2*B21)^(N1-1)*A1*B12*A2*B2c;
    % Método da Matriz de Transferência (MMT) para o espelho M1
    M01 = B01*(A1*B12*A2*B21)^(N1-1)*A1*B12*A2*B2t; 

    % MMT para o espelho M2
    M02 = B01*(A1*B12*A2*B21)^(N2-1)*A1*B12*A2*B2t; 

    %============================================================
    % Espelho M2                                                x
    %   nc=|-n2-|-|-n1-|...|-n2-|-|-n1-|-|-n2-|-|-n1-|-n_out    x
    %    Bc2-A2-B21-A1-B12...A2-B21-A1-B12-A2-B21-Al-B2t        x
    %   Bc2-(A2-B21-A1-B12)...(A2-B21-A1-B12)-A2-B21-A1-B2t     x
    %   Bc2-(A2-B21-A1-B12)^(N2-1)-A2-B21-A1-B2t                x
    %============================================================
    M2 = Bc2*(A2*B21*A1*B12)^(N2-1)*A2*B21*A1*B1t;

    % MTT para cavidade com meio ativo entre M1, M2
    Md = M1*Ac*M2;

    % Reflexão e transmissão de M1
    R01(k) = abs(M01(2,1)/M01(1,1))^2;
    R02(k) = abs(M02(2,1)/M02(1,1))^2; 

    % Reflexão e transmissão de M2
    R(k) = abs(Md(2,1)/Md(1,1))^2;
    T(k) = n_out/n_in*abs(1/Md(1,1))^2; 

end

atick1 = 450:50:650;
figure(1)
  plot(wavelength, R01, 'r', wavelength, R02, '-.b', 'linewidth', 2)
  legend('R_1', 'R_2');
  set(gca, 'FontSize', 18);
  axis([450 650 0 1.05]);
  set(gca, 'XTick', atick1);
  set(gca, 'YTick', 0:0.2:1);
  grid minor
  xlabel('Wavelength (nm)');
  ylabel('R & T');
  title('Reflexão R_1 e R_2');

  figure(2)
  plot(wavelength, T, 'g','linewidth', 2)
  legend('T');
  set(gca, 'FontSize', 18);
  axis([450 650 0 1.05]);
  set(gca, 'XTick', atick1);
  set(gca, 'YTick', 0:0.2:1);
  grid minor
  xlabel('Wavelength (nm)');
  ylabel('T');
  title('Modos na Cavidade: L = 2 cm');
%==========================================================================