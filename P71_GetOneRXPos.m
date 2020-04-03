function UPos = P71_GetOneRXPos(Es, inGPSTime, inTimeShifts, ...
    SampleNums, Params)
%
% Функция расчёта одного набора координат приёмника
%
% Входные переменные
%   Es - cell-массив с эфемеридами спутников;
%   inGPSTimes - моменты времени, в которые были испущены сигналы со
%       спутников;
%   inTimeShifts - отличия значений задержки распространения сигналов
%       спутников от общей постоянной составляющей задержки
%       распространения;
%   SampleNums - номера отсчётов, в которые пришли сигналы спутников с
%       метками времени inGPSTimes.
%
% Выходные переменные
UPos = struct('x',[]);
%   UPos - результат-структура с полями:
%       x, y, z - координаты в прямоугольной ДСК
%       T0 - общая составляющая сдвига по времени
%       tGPSs, SampleNums - значения времени GPS для заданных отсчётов 
%       Lat, Lon, Alt - широта, долгота, высота
%       SatsPoses - координаты спутников, таблица со столбцами
%           x, y, z - координаты в прямоугольной ДСК;
%           xs_k, ys_k, i_k - координаты перед преобразованиями СК;
%           Lat, Lon, Alt - широта, долгота, высота;
%           El, Az - угол склонения и азимут;
%       NumIters, MaxNumIters - выполненное и максимальное число итераций;
%       Delta, MaxDelta - достигнутое и максимальное значение оценки
%           изменения координат пользователя между соседними итерациями
%           (м);
%       inGPSTimes, GPSTimes, inTimeShifts, TimeShifts - сохранение
%           параметров и скорректированных параметров.

%% УСТАНОВКА ПАРАМЕТРОВ
    % Максимальное число итераций
        MaxNumIters = Params.P71_GetOneRXPos.MaxNumIters;
    % Максимальное изменение координат пользователя между соседними
    % итерациями (м). Если фактическое изменение меньше, то цикл
    % останавливается
        MaxDelta = Params.P71_GetOneRXPos.MaxDelta;

%% УСТАНОВКА КОНСТАНТ
    % Скорость света, м/с
        c = 299792458;
    % Радиус Земли, м
        R = 6356863;
        %Вычисляются координаты спутника первый раз
        T_0 = 68*10^-3;
%  inGPSTime = [];       
% inGPSTime = [405444, 405443.995572825, 405443.994012708, 405443.996125122, 405443.988963343, 405443.992138319];

for i = 1:size(Es,2)
    inTProp(i) = T_0+inTimeShifts(i);
    [SatPos(i), GPSTime(i), TProp(i)] = P72_GetSatPos(Es{i},  inGPSTime(i), inTProp(i),Params);
    x_s(i) = SatPos(i).x;
    y_s(i) = SatPos(i).y;
    z_s(i) = SatPos(i).z;
    TimeShifts(i) = inTimeShifts(i) + TProp(i)-inTProp(i);
end



r = 1*10^-2;
% чтоб в цикл зайти
dx=100;
dy = 100;
dz = 100;
dT_0 = 100;

   


%Начальные координаты приемника на шаге итерации 0
x_m = mean(x_s);
z_m = mean(z_s);
y_m = mean(y_s);

x_0(1) = R*(x_m/sqrt(x_m^2+y_m^2+z_m^2));
y_0(1) = R*(y_m/sqrt(x_m^2+y_m^2+z_m^2));
z_0(1) = R*(z_m/sqrt(x_m^2+y_m^2+z_m^2));

x_0(2) = -R*(x_m/sqrt(x_m^2+y_m^2+z_m^2));
y_0(2) = -R*(y_m/sqrt(x_m^2+y_m^2+z_m^2));
z_0(2) = -R*(z_m/sqrt(x_m^2+y_m^2+z_m^2));

if sqrt((x_0(1)-x_m)^2+(y_0(1)-y_m)^2+(z_0(1)-z_m)^2)<sqrt((x_0(2)-x_m)^2+(y_0(2)-y_m)^2+(z_0(2)-z_m)^2)
   x=x_0(1);
   y=y_0(1);
   z=z_0(1); 
else 
   x=x_0(2);
   y=y_0(2);
   z=z_0(2); 
end

% % %     for i = 1:size(Es,2)
% % %  t_t(i) = sqrt((x-x_s(i))^2+(y-y_s(i))^2+(z-z_s(i))^2) ;
% % %      OutSatPos(i) = P73_RenewSatPos(SatPos(i), TProp(i), Params);
% % % % %     OutSatPos(i) = P73_RenewSatPos(SatPos(i), t_t(i), Params);
% % %     x_s(i) = OutSatPos(i).x;
% % %     y_s(i) = OutSatPos(i).y;
% % %     z_s(i) = OutSatPos(i).z;
% % %     SatPos(i).x = x_s(i);
% % %     SatPos(i).y = y_s(i);
% % %     SatPos(i).z = z_s(i);
% % %     end
count = 0;
while sqrt(dx^2+dy^2+dz^2+(c*dT_0)^2) >= r
    
    
    
    B=[];
    A = [];
    for ii = 1:size(Es,2)
        T(ii) = sqrt((x_s(ii)-x)^2+(y_s(ii)-y)^2+(z_s(ii)-z)^2)/c-T_0;
        % % %         B = [B;(inTimeShifts(ii)-T(ii))];
        B = [B;(TimeShifts(ii)-T(ii))];
        A = [A;[-(x_s(ii)-x)/(c*(T(ii)+T_0)),-(y_s(ii)-y)/(c*(T(ii)+T_0)),-(z_s(ii)-z)/(c*(T(ii)+T_0)), -c]];
    end
    B = c.*B;
    
    X =A\B;
    
    dx = X(1);
    dy = X(2);
    dz = X(3);
    dT_0 = X(4);
    if sqrt(dx^2+dy^2+dz^2+(c*dT_0)^2) >= r
    % здесь докручиваем Землю
    % но если уже сошлось, то ничего не прибавляем и не докручиваем
        x = x + dx;
        y = y + dy;
        z = z + dz;
        T_0 = T_0 + dT_0;
        % % % TProp = TProp +  dT_0;
        for i = 1:size(Es,2)
            
            OutSatPos(i) = P73_RenewSatPos(SatPos(i), T(i), Params);
            
            x_s(i) = OutSatPos(i).x;
            y_s(i) = OutSatPos(i).y;
            z_s(i) = OutSatPos(i).z;
            SatPos(i).x = x_s(i);
            SatPos(i).y = y_s(i);
            SatPos(i).z = z_s(i);
        end
        
    end
    for i = 1:size(Es,2)
        [SatPos(i).El, SatPos(i).Az] = P75_CalculateSatElAz([SatPos(i).x,SatPos(i).y,SatPos(i).z], [x,y,z], Params);
        [SatPos(i).Lat, SatPos(i).Lon, SatPos(i).Alt] = P74_Cartesian2Spherical([SatPos(i).x,SatPos(i).y,SatPos(i).z], Params);
    end
    count = count+1;
    % disp(B);
end
UPos.x = x;
UPos.y = y;
UPos.z = z;
UPos.T_0 = T_0;
UPos.inGPSTime = inGPSTime;
UPos.GPSTime = GPSTime;
UPos.TProp = TProp;
UPos.inTimeShifts = inTimeShifts;
UPos.SampleNums = SampleNums;
UPos.NumIters = count;
UPos.SatPos = SatPos;

[UPos.Lat, UPos.Lon, UPos.Alt] = P74_Cartesian2Spherical([x,y,z], Params);
