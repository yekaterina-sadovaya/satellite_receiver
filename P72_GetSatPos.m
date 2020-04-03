function [SatPos, GPSTime, TProp] = P72_GetSatPos(Data,  inGPSTime, inTProp,Params) %#ok<INUSD>
% inGPSTime - t, tc
% Фунция производит вычисление координат спутника в момент времени
% inGPSTime и учитывает время распространения сигнала inTProp при переводе
% координат в систему ECEF
%
% Входные переменные
%   Data - структура, содержащая, как минимум, параметры подкадров 1, 2 и
%     3;
%   inGPSTime - время испускания сигнала;
%   inTProp - время распространения сигнала.
%
% Выходные переменные
%   SatPos - массив (8х1) координат и параметров спутника для пересчёта
%       координат:
%         [x; y; z; ... % координаты в прямоугольной системе координат
%         xs_k; ys_k; i_k; ... % исходные координаты спутника
%         Omega_k; % исходное значение Omega_k
%         ZaZa]; % параметр для пересчёта Omega_k в Omega_k_TProp
%   GPSTime - скорректированное время испускания сигнала;
%   TProp - скорректированное время распространения сигнала.

% страница 104

SatPos =  struct( 'x', [], 'y', [], 'z',[], 'xs_k', [], 'ys_k',[], 'i_k', [], 'Omega_k', []);

mu = 3.986005*10^14;
OMEGA_e = 7.2921151467*10^-5;
A = Data.sqrtA^2;
n_0 = sqrt(mu/A^3);
t_c = inGPSTime;
% % % t_k = t_c-Data.t_oe;

if t_c-Data.t_oe > 302400
% % %     t_k = t_c-Data.t_oe - 604800;
         t_c = t_c - 604800;
end
if t_c-Data.t_oe < -302400
% % %     t_k = t_k + 604800;
            t_c = t_c + 604800;
end
t_k = t_c-Data.t_oe;
n = n_0 + Data.Delta_n*pi;
M = Data.M_0*pi + n*t_k;

% M = E - Data.e*sin(E);
%Решаем это уравнение

Eold = M;
error = 1;
while error > 1*10^-12
      E = M + Data.e*sin(Eold);
      error  = abs(E-Eold);
      Eold = E;
end
%Это мы нашли E с использованием времени t_k, но это грубая оценка, его
%нужно скорректировать из-за влияния различных факторов на форму орбиты
%спутника

F = -4.442807633*10^-10;
delta_t_r = F * Data.e*Data.sqrtA*sin(E);
 delta_t = Data.a_f0 + Data.a_f1*(t_c - Data.t_oc)+Data.a_f2*(t_c - Data.t_oc)^2+delta_t_r;
t_corrected = t_c - delta_t; %оно просто t в книге
    
% nu =atan((sqrt(1-Data.e^2)*sin(E)/(1-Data.e*cos(E)))/((cos(E)-Data.e)/(1-Data.e*cos(E))));
v1 = acos((cos(E)-Data.e)/(1-Data.e*cos(E)));
v2 = asin(((1-Data.e^2)^0.5)*sin(E)/(1-Data.e*cos(E)));
nu = v1*sign(v2);
%E = acos((Data.e+cos(nu))/(1+Data.e*cos(nu)));Зачем нужно это
%уравнение????????????

Phi = nu + Data.omega*pi;
delta_u = Data.C_us*sin(2*Phi) + Data.C_uc*cos(2*Phi);
delta_r = Data.C_rs*sin(2*Phi) + Data.C_rc*cos(2*Phi);
delta_i = Data.C_is*sin(2*Phi) + Data.C_ic*cos(2*Phi);

u = Phi + delta_u;
r = A*(1-Data.e*cos(E)) + delta_r;
i = Data.i_0*pi + delta_i + Data.IDOT*(t_corrected-Data.t_oe);
%  i = Data.i_0 + delta_i + Data.IDOT*(t_k);
x1 = r*cos(u);
y1 = r*sin(u);

OMEGA = Data.Omega_0*pi + (Data.DOmega*pi-OMEGA_e)*(t_corrected-Data.t_oe) - OMEGA_e*Data.t_oe;
%  OMEGA = Data.Omega_0 + (Data.DOmega-OMEGA_e)*(t_k) - OMEGA_e*Data.t_oe;

x = x1*cos(OMEGA) - y1*cos(i)*sin(OMEGA);
y = x1*sin(OMEGA) + y1*cos(i)*cos(OMEGA);
z = y1*sin(i);

SatPos.x = x;
SatPos.y = y;
SatPos.z = z;
SatPos.xs_k = x1;
SatPos.ys_k = y1;
SatPos.i_k = i;
SatPos.Omega_k = OMEGA;
GPSTime = t_corrected;
TProp = inTProp + delta_t;%сигнал испускается раньше на дельта t, значит время распростронения будет больше на дельта t, это из документа
end







