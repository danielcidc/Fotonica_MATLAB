%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%  DFB_1C.m   Modelo de cavidade DFB (Distributed Feedback laser)                  x
%  Matriz de Transferência para DFB com 2 Grades (Gratings) e 1 Cavidade           x
%  n = 1.5634, Delta_n = 1e-4                                                      x
%  L1 = 20 mm (N1 ~ 41000); L2 ~ 15mm (N2 ~ 30000)                                 x
%  Lc1 = 50 micron                                                                 x
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
format long

N1 = 41000;                              % Nº de períodos
N2 = 30000;
neff = 1.5634;                           % Índice de refração efetivo
d1 = 1.5352/(4*neff);
d2 = 1.5352/(4*neff);
dB = d1 + d2;                            % Comprimento de um período em micrômetros
lambdaB = 2*neff*dB;                      % Comprimento de onda de Bragg em micrômetros
Lc = 50.0;                               % Comprimento da cavidade em micrômetros 

% Indíce
n1 = neff;                               % Índice efetivo de modo de propagação 
dn = 1*1e-4;                             % Delta n
n2 = n1 - dn;
nc = n1;                                 % neff na cavidade = n1

n0 = 1.0;
n_out = 1.0;

% Interface de entrada (primeira da esquerda ar-vidro): n0-|n1 
n01_p = (n0 + n1)/(2*n0); n01_m = (n0 - n1)/(2*n0);
B01 = [n01_p, n01_m; n01_m, n01_p]; 

% Na grade: n1-|n2
n12_p = (n1 + n2)/(2*n1); n12_m = (n1 - n2)/(2*n1);
B12 = [n12_p, n12_m; n12_m, n12_p]; 

% Na grade: n2-|n1
n21_p = (n2 + n1)/(2*n2); n21_m = (n2 - n1)/(2*n2);
B21 = [n21_p, n21_m; n21_m, n21_p];

% Interface de saída (substrato vidro-ar mais à direita): n1-|n_out
n1t_p = (n1 + n_out)/(2*n1); n1t_m = (n1 - n_out)/(2*n1);
B1t = [n1t_p, n1t_m; n1t_m, n1t_p]; 

% Matriz unitária
B11 = [1, 0; 0, 1];

%------------------------------------------------------------------------
%
%    --->-||||||||||||||||||===||||||||||||||||-->-
%  
%                M1         c1        M2   
%
%-------------------------------------------------------------------------
dfb2 = fopen('rt1c.dat','w+');         % salvando dados

wavelength = 1.53505:0.0000001:1.53520;
R01 = zeros(1,length(wavelength));
R02 = zeros(1,length(wavelength));
R = zeros(1,length(wavelength));
T = zeros(1,length(wavelength));

for k = 1:length(wavelength)
    lambda = wavelength(k);

% Grades 
    alfa1 = exp(1i*2*pi*n1*d1/lambda);  
    alfa2 = exp(1i*2*pi*n2*d2/lambda);  
    A1 = [alfa1, 0; 0, 1/(alfa1)];
    A2 = [alfa2, 0; 0, 1/(alfa2)];
    
    %=============================================================    
    % Primeira Grade M1                                          x
    % Ar --|-n1-|-(|-n2-|-|-n1-|)...-(|-n2-|--|-n1-|)-n1        x
    %     B01-A1-(B12-A2-B21-A1)-...-(B12-A2-B21-A1)-B11         x
    %     B01-A1-(B12-A2-B21-A1)^Na-B11                          x
    %=============================================================
    M1 = B01*A1*(B12*A2*B21*A1)^N1*B11;  
    %=============================================================    
    % Segunda Grade M2                                           x
    % n1--|n1|---(|n2|---|n1|)-...-(|n2|---|n1|--n_out           x
                                                                 
    %    B11-A1-(B12-A2-B21-A1)-...-(B12-A2-B21-A1)-B1t          x
    %    B11-A1-(B12-A2-B21-A1)^N2-B1t                           x
    %=============================================================
    M2 = B11*A1*(B12*A2*B21*A1)^N2*B11;
    % Matriz de propagação na cavidade Lc
    alfaC = exp(1i*2*pi*n1*Lc/lambda);  
    AC = [alfaC, 0; 0, 1/(alfaC)];
    % Matriz de transferência para FG1 & FG2 
    M01 = B01*A1*(B12*A2*B21*A1)^N1*B1t; 
    R01(k) = abs(M01(2,1)/M01(1,1))^2;                 
    % FG1's R

    M02 = B01*A1*(B12*A2*B21*A1)^N2*B1t; 
    R02(k) = abs(M02(2,1)/M02(1,1))^2;

    % Matriz de Transferência de toda a cavidade
    M = M1*AC*M2*B1t;
    R(k) = abs(M(2,1)/M(1,1))^2;
    T(k) = n_out/n0*abs(1/M(1,1))^2;
fprintf(dfb2,'%f %12.8f %12.8f %12.8f12.8f\n', lambda,R01(k),R02(k),R(k),T(k));
end

atick1 = 1.5351:0.00005:1.53520;
figure(1)
    plot(wavelength,R01,'-.r', wavelength,R02,'-.b','linewidth', 3)
    legend('R_1', 'R_2');
    set(gca,'FontSize',24);
    axis([1.5351 1.53520 0 1.05]);
    set(gca,'XTick',atick1);
    set(gca,'YTick',0:0.2:1);
    grid minor
    xlabel('Wavelength (\mum)');
    ylabel('Reflectivity R');
    title( 'Gradeamento R_1, R_2');

figure(2)
    plot(wavelength,R,'r', wavelength,T,'b','linewidth', 3)
    legend('R', 'T');
    set(gca,'FontSize',24);
    axis([1.5351 1.53520 0 1.05]);
    set(gca,'XTick',atick1);
    set(gca,'YTick',0:0.2:1);
    grid minor
    xlabel('Wavelength (\mum)');
    ylabel('R & T');
    title( 'Cavidade DFB monomodo');

fclose(dfb2);
%%===================== fim =========================================