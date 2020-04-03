function Res = P60_GatherSatsEphemeris(inRes, Params) %#ok<INUSD>

% Функция сбора навигационной информации для спутников, у которых было
% найдено хотя бы одно значение TOW_Count_Message
%
% Входные переменные
%   inRes - структура с результатами модели, объявленная в Main;
%
% Выходные переменные
%   Res - структура, которая отличается от inRes добавлением нового поля,
%       описание которого дано ниже в коде.

% Пересохранение результатов
Res = inRes;

%% ИНИЦИАЛИЗАЦИЯ РЕЗУЛЬТАТА
% Ephemeris = cell(N, Res.Search.NumSats);

% Количество строк cell-массива Ephemeris совпадает с количеством
% подкадров спутников с одинаковым значением TOW (естественно,
% учитываются только те спутники, у которых SatsData.isSat2Use = 1).
% Элементами Ephemeris являются структуры, содержащие значения всех
% параметров первого, второго и третьего подкадров, а также порядковый
% номер подкадра, номер первого CA подкадра, передаваемое значение TOW,
% для которого верна эта информация. 
%?????????????Если полную информацию собрать не
% удалось, то элемент cell-массива должен быть пустым.

%% УСТАНОВКА ПАРАМЕТРОВ

%% РАСЧЁТ ПАРАМЕТРОВ
% Имена всех полей структур, являющихся элементами Ephemeris
ENames = { ...
    ... % Эти поля не относятся к навигационной информации
    'SFNum', ... % Порядковый номер подкадра спутника,
    ... % соответствующего текущей строке (подкадру) Ephemeris
    'CANum', ... % Номер CA-кода спутника, с которого начинается
    ... % подкадр с порядковым номером SFNum
    'TOW', ... % Значение TOW, передаваемое в подкадре с порядковым
    ... % номером SFNum. Это значение одинаковое для всех
    ... % элементов одной строки Ephemeris
    ...
    ... % Поля с навигационной информацией
    'WeekNumber', ...
    'CodesOnL2', ...
    'URA', ...
    'URA_in_meters', ...
    'SV_Health_Summary', ...
    'SV_Health', ...
    'IODC', ...
    'Flag_L2_P_Code', ...
    'T_GD', ...
    't_oc', ...
    'a_f2', ...
    'a_f1', ...
    'a_f0', ...
    'IODE', ...
    'C_rs', ...
    'Delta_n', ...
    'M_0', ...
    'C_uc', ...
    'e', ...
    'C_us', ...
    'sqrtA', ...
    't_oe', ...
    'Fit_Interval_Flag', ...
    'AODO', ...
    'C_ic', ...
    'Omega_0', ...
    'C_is', ...
    'i_0', ...
    'C_rc', ...
    'omega', ...
    'DOmega', ...
    'IODE', ...
    'IDOT', ...
    };

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
% Строка состояния

% Определим порядковые номера спутников, для котрых мы будем пытаться
% собирать эфемериды
i = 1;
iiii = 1;
iiiii = 1;
for k = 1:Res.Search.NumSats
    if Res.SatsData.isSat2Use
        for kk = 1:size(Res.SatsData.HOW{k},2)
            Ephemeris{kk,i}.(ENames{1}) = kk; % порядковый номер сабфрейма (но не его действительный номер)
            Ephemeris{kk,i}.(ENames{2}) = Res.BitSync.CAShifts(k) + Res.SubFrames.BitShift(k)*20 + (kk-1)*20*300+1; % порядковый номер СА кода - мы хотим в СА выразить шифт до каждого подкадра
            Ephemeris{kk,i}.(ENames{3}) = Res.SatsData.HOW{k}(kk).TOW;
            
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  1
                for ii = kk:kk+4
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 4:16
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF1{k}(kk).(ENames{iii});
                        end
                    end
                end
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  2
                for ii = kk-1:kk+3
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 17:27
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF2{k}(kk).(ENames{iii});
                        end
                    end
                end
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  3
                for ii = kk-2:kk+2
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 28:36
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF3{k}(kk).(ENames{iii});
                        end
                    end
                end
            end 
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  4 & kk-3<0
               fill_gaps_4SF(iiii,1) =   k;
               fill_gaps_4SF(iiii,2) =   kk;
               iiii = iiii +1;
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  5 & kk-4<0
               fill_gaps_5SF(iiiii,1) =   k;
               fill_gaps_5SF(iiiii,2) =   kk;
               iiiii = iiiii+1;
            end
        end     
    end
    

                
    % Для каждого спутника определим порядковый номер подкадра, в котором
    % встречается первое значение TOW общее с остальными спутниками
    
    % Заготовим результат
    
    % Теперь попробуем для каждого подкадра каждого спутника определить
    % эфемериды
    
    % Добавим новое поле с результатами в Res
    
    % Строка состояния
    i = i+1;
end
    %Тут конечно нужно сравнить IODC и IODE и для SF4 сделать
                for i = 1:size(fill_gaps_5SF,1)
                        for iii = 4:16
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF1{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+1).(ENames{iii});
                        end
                        for iii = 17:27
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF2{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+2).(ENames{iii});
                        end

                        for iii = 28:36
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF3{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+3).(ENames{iii});
                        end
                end
Res.Ephemeris = Ephemeris;

end



function E = MakeEmptyE(ENames)
    % Создадим все поля

    % Установим в поля, не относящиеся к навигационным данным, произвольные
    % параметры, чтобы тест isGathered проходил успешно (см. CheckAndAddE)
        E.SFNum = -1;
        E.CANum = -1;
        E.TOW   = -1;
end

function [outE, isNew] = CheckAndAddE(inE, SFNum, SFData, ENames)

% В зависимости от номера подкадра мы сравниваем значение либо IODC, либо
% IODE, имеющееся в InE с тем же значением в SFData, потом сравниваем
% значения IODC и IODE в inE

% Если нужно создавать новую E, то сделаем это, в противном случае
% скопируем E со входа
        
% Обновим пустые поля outE значениями из SFData
end

function isGathered = CheckE(E, ENames)
% Проверим, остались ли пустые поля
end