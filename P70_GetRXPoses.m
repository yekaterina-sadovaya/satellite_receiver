function Res = P70_GetRXPoses(inRes, Params)
%который превышает подкадр и numsubframe будет не работать
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
% Positioning = struct(  ...
%   'RXPoses', cell(N, M), ...
%   'CAStep', CAStep, ...
%   'isCommonRxTime', isCommonRxTime ...
% );

% Количество строк N cell-массива RXPoses совпадает с количеством
%   подкадров спутников с одинаковыми значениями TOW. Количество
%   столбцов определяется количеством вычислений координат на
%   длительности одного подкадра и зависит от параметра CAStep,
%   определяемого ниже.
% CAStep - шаг в периодах CA-кода между соседними вычислениями
%   координат.
% isCommonRxTime - параметр, определяющий как вычислять параметры
%   спутников: в одно время приёмника или в разное.

%% УСТАНОВКА ПАРАМЕТРОВ
% Шаг в периодах CA-кода между соседними вычислениями координат. Всего
% в подкадре 6000 периодов CA-кода, поэтому, например, CAStep = 1000
% приведёт к вычислению 6 координат за один подкадр.
CAStep = Params.P70_GetRXPoses.CAStep;
CALen = 1023 * Res.File.R;
% Вариант вычисления координат.
% isCommonRxTime = 1 - координаты спутников вычисляются в одинаковый
%   момент  времени приёмника, соответствующий разным меткам
%   времени GPS
% isCommonRxTime = 0 - координаты спутников вычисляются в разные
%   моменты времени приёмника, соответствующие одинаковой метке
%   времени GPS
isCommonRxTime = Params.P70_GetRXPoses.isCommonRxTime;

% Порядковые номера спутников, учитываемых при вычислении координат:
% 'all' - все спутники;
% 'firstX' - первые Х спутников, например 'first5';
% [1, 2, 5, 7] - конкретные номера.
SatNums2Pos = Params.P70_GetRXPoses.SatNums2Pos;


Positioning = cell(size(Res.Ephemeris,1), 6000/CAStep);

% struct(  ...
%   'RXPoses', cell(size(Res.Ephemeris,1), 6000/CAStep), ...
%   'CAStep', CAStep, ...
%   'isCommonRxTime', isCommonRxTime ...
% );
%% РАСЧЁТ ПАРАМЕТРОВ
% Интервал дискретизации сигнала
dt = 1/Res.File.Fs;

% Определим конкретные номера спутников
if ischar(SatNums2Pos)
    if strcmp(SatNums2Pos, 'all')
        CurSatNums2Pos = 1:Res.Search.NumSats;
    else
        Buf = str2double(SatNums2Pos(6:end));
        CurSatNums2Pos = 1:Buf;
    end
else
    CurSatNums2Pos = SatNums2Pos;
end

%% РАСЧЁТ КООРДИНАТ
%Определяем стартовый шивт, смотрим на шифты СА кода первых подкадров,
%идем по первой строчке, потому что в накопе эфемерид мы сделали так, что
%HOW одинакова во всех этих подкадрах первой строки

kk = 1;
NumSubframe = 1;
for k = CurSatNums2Pos
    CAShifs(kk) = Res.Ephemeris{NumSubframe,k}.CANum;
    kk = kk+1;
end

[~,min_index]=min(CAShifs);
NumSubframe = 1;

n=1;
for N = 1:size(Res.Ephemeris,1) % N - цикл по количеству подкадров
    for M = 1:6000/CAStep % 6000/CAStep вычисляет сколько раз вычислим позицию за подкадр 
        
        SampleNum = Res.Track.SamplesShifts{min_index}(min(CAShifs))+1+(n-1)*CALen*CAStep; % Res.Track.SamplesShifts{min_index}(min(CAShifs)) берем порядковый номер
        
        kk=1;
        Es = [];
        
        for k = CurSatNums2Pos
            
            RefTOW(kk) = Res.Ephemeris{N,k}.TOW;
            RefCANum(kk) = Res.Ephemeris{N,k}.CANum;
            
            TCA = 10^-3;
            [~,NumberCA(kk)] = min(abs(Res.Track.SamplesShifts{k}+1-SampleNum)); % порядковый номер СА
            SampleNums(kk) = Res.Track.SamplesShifts{k}(NumberCA(kk)) ; % шифт в сэмплах до этого СА
            d(kk) = Res.Track.SamplesShifts{k}(NumberCA(kk)) +1 - SampleNum;
            inGPSTime(kk) = GettGPS (RefCANum(kk), RefTOW(kk), dt,TCA,NumberCA(kk),d(kk)); %оно же t_c
            
            Es{kk} = Res.Ephemeris{N,k};
            kk = kk+1;
        end
        
        % находим первый СА код в подкадре (референс) т.к. мы не можем
        % просто по номеру СА смотреть. Потому что у нас запись и там по
        % ним рассинхронизация -- номер одного СА не соответствует СА
        % другого спутника. А подкадр точно испустился синхронно
        [min_value_ref,min_index_ref] = min(RefCANum);
        % min_value_ref - нaчало подкадра спутнника с минимальным СА. 
        % значит нужно искать разницу текущего номера СА и референсного
        diff_CA =  NumberCA(min_index_ref) - min_value_ref ;
        kk=1;
        % inTimeShifts - это T(i)
        for k = CurSatNums2Pos
            inTimeShifts(kk) = (Res.Track.SamplesShifts{k}(RefCANum(kk)+diff_CA) +1- SampleNum).*dt;
            kk = kk+1;
        end
        
        UPos = P71_GetOneRXPos(Es, inGPSTime, inTimeShifts, ...
            SampleNums, Params);
        Positioning{N,M} = struct(  ...
            'RXPoses', struct(), ...
            'CAStep', CAStep, ...
            'isCommonRxTime', isCommonRxTime ...
            );
        
        Positioning{N,M}.RXPoses = UPos;
        
        RXPoses{N,M} = UPos;
    
        n=n+1;   
    end
end

Res.Positioning = Positioning;

P76_ExportResults(RXPoses, Params);

end

function inGPSTime = GettGPS(RefCANum, RefTOW, dt,TCA,NumberCA,d)

% Функция определяет tGPS для отсчёта сигнала SampleNum

% Входные переменные
%   SampleNum - номер отсчёта записи, для которого надо расчитать время
%       GPS;
%   SamplesNums - номера первых отсчётов CA-кодов текущего спутника;
%   RefCANum - номер CA-кода, который является первым в подкадре, в котором
%       передаётся значение RefTOW;
%   RefTOW - значение RefTOW;
%   dt - интервал дискретизации записи.
%
% Выходные переменные
%   tGPS - время GPS в отсчёт SampleNum.

%  t = (TOW-1) + две дельты т
% TOW дается на конец подкадра в котором он передается
% находим время отправки как сумму TOW; времени от начала подкадра (по 1-му
% СА) до начала второго времени; второе время определили ранее когда
% считали разность массива шифт и текущего сэмпла

% Константы

inGPSTime = (RefTOW-1)*6 + (NumberCA - RefCANum)*TCA - d*dt;

end