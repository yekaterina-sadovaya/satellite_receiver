function Res = P10_NonCohSearchSats(inRes, Params)
% 1) Посмотреть что будет если брать не модуль а просто квадраты
% Функция некогерентного поиска спутников в файле-записи
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
Search = struct( ...
    'NumSats',       [0], ... % Скаляр, количество найденных спутников
    'SatNums',       [], ... % массив 1хNumSats с номерами найденных
    ... % спутников
    'SamplesShifts', [], ... % массив 1хNumSats, каждый элемент -
    ... % количество отсчётов, которые нужно пропустить в файле-
    ... % записи до начала первого периода CA-кода соответствующего
    ... % спутника
    'FreqShifts',    [], ... % массив 1хNumSats со значениями частотных
    ... % сдвигов найденных спутников в Гц
    'CorVals',       [], ... % массив 1хNumSats вещественных значений
    ... % пиков корреляционных функций нормированных на среднее
    ... % значение, по которым были найдены спутники
    'AllCorVals',    zeros(1, 32) ... % массив максимальных значений
    ... % всех корреляционных функций
    );

%% УСТАНОВКА ПАРАМЕТРОВ
% Количество периодов, учитываемых при обнаружении.
NumCA2Search = Params.P10_NonCohSearchSats.NumCA2Search;
%         NumCA2Search = 1;
% Массив центральных частот анализируемых диапазонов, Гц
CentralFreqs = Params.P10_NonCohSearchSats.CentralFreqs;

% Порог обнаружения
SearchThreshold = Params.P10_NonCohSearchSats.SearchThreshold;

%% СОХРАНЕНИЕ ПАРАМЕТРОВ
Search.NumCA2Search    = NumCA2Search;
Search.CentralFreqs    = CentralFreqs;
Search.SearchThreshold = SearchThreshold;

%% РАСЧЁТ ПАРАМЕТРОВ
% Количество рассматрвиаемых частотных диапазонов
NumCFreqs = length(CentralFreqs);
%если частата дискретизации больше, нужно тут задавать, потому что в файле записано с ошибкой
% Res.File.R = 4;
% Длина CA-кода с учётом частоты дискретизации
CALen = 1023 * Res.File.R;
%Флаг использования корня при расчете модуля, проверяем что будет
FlagSqrt = 1;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ
Search.NumSats = 0;

[Signal,~] = ReadSignalFromFile(Res.File, 0, NumCA2Search*CALen+CALen-1);

dt=10^-3/CALen;
n = 0:CALen-1;

%цикл по количеству спутников
for i = 1:32
    correlation_result = [];
    summa = [];
    norm_summa = [];
    Signal_shift = [];
    CACode_one=[];
    CACode = [];
    CACode_one = GenCACode(i, 1);
    CACode_one(CACode_one==0)=-1;
    for r = 1:Res.File.R
        CACode(r,:) = CACode_one;
    end

    
    for k = 1:NumCFreqs
         CACode_shift = [];
        CACode_shift = reshape(CACode,1,Res.File.R*size(CACode,2)).*exp(-n*1j*2*pi*CentralFreqs(k)*dt);
        %можно было сдвигать сигнал, но я сдвигала СА код
        Signal_shift(k,:) = Signal;
        correlation_result(k,:) = conv(Signal_shift(k,:),fliplr(CACode_shift),'valid');
            if NumCA2Search~=1
                if FlagSqrt == 1
                    summa(k,1:CALen) = sum(abs(reshape(correlation_result(k,1:CALen*NumCA2Search),CALen,NumCA2Search)'));
                else
                    ImValues = reshape(correlation_result(k,1:CALen*NumCA2Search),CALen,NumCA2Search);
                    SqrImValue = ImValues.*conj(ImValues); % это просто сумма квадратов
                    summa(k,1:CALen) = sum(SqrImValue');
                end
            else
                if FlagSqrt == 1
                    summa(k,:) = abs(correlation_result(k,:));
                else
                    ImValues = correlation_result(k,:);
                    SqrImValue = ImValues.*conj(ImValues);
                    summa(k,:) = sum(SqrImValue');
                end
            end
    end

    norm_summa = summa/mean(mean(summa));

    max_norm_summa = max(max(norm_summa));
 
    %все записываю в матрицу а потом по максимальной мощности сортирую
    %столбцы
    if  max_norm_summa > SearchThreshold
        [frequency_point,sample_point] = find(norm_summa == max_norm_summa);
        matrix(Search.NumSats+1,1) = max_norm_summa;
        matrix(Search.NumSats+1,2) = mod(sample_point,CALen)-1;
        matrix(Search.NumSats+1,3) = CentralFreqs(frequency_point);
        matrix(Search.NumSats+1,4) = i;
        Search.NumSats = Search.NumSats + 1;
    end
    Search.AllCorVals(i) =  max_norm_summa;
    
        figure();
        surf(1:size(norm_summa,2),CentralFreqs,norm_summa);
end
%тут сортирую
if size(matrix,1)>0
    matrix = sortrows(matrix, -1); % same as sortrows(matrix, 1, 'descend')
    Search.FreqShifts = matrix(:,3)';
    Search.SamplesShifts = matrix(:,2)';
    Search.SatNums = matrix(:,4)';
    Search.CorVals =  matrix(:,1)';
end

Res.Search  = Search;

figure ();
plot(Search.AllCorVals,'*-');
hold on;
plot([1,length(Search.AllCorVals)],[SearchThreshold,SearchThreshold])
grid on;
