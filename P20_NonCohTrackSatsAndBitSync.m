function Res = P20_NonCohTrackSatsAndBitSync(inRes, Params)
% Функция некогерентного трекинга спутников и битовой синхронизации
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
Track = struct( ...
    'SamplesShifts', {cell(Res.Search.NumSats, 1)}, ...
    'CorVals',       {cell(Res.Search.NumSats, 1)} ...
    );
% Каждая ячейка cell-массивов SamplesShifts и CorVals является массивом
%   1xN, где N - количество периодов CA-кода соответствующего спутника,
%   найденных в файле-записи (N может быть разным для разных
%   спутников).
% Каждый элемент массива SamplesShifts{k} - количество отсчётов,
%   которые надо пропустить в файле-записи до начала соответствующего
%   периода CA-кода.
% Каждый элемент массива CorVals{k} - комплексное значение корреляции
%   части сигнала, содержащей соответствующий период CA-кода, с опорным
%   сигналом.

BitSync = struct( ...
    'CAShifts', zeros(Res.Search.NumSats, 1), ...
    'Cors', zeros(Res.Search.NumSats, 19) ...
    );
% Каждый элемент массива CAShifts - количество периодов CA-кода,
%   которые надо пропустить до начала бита.
% Каждая строка массива Cors - корреляции, по позиции минимума которых
%   определяется битовая синхронизация.

%% УСТАНОВКА ПАРАМЕТРОВ
% Количество периодов CA-кода между соседними синхронизациями по
% времени (NumCA2NextSync >= 1, NumCA2NextSync = 1 - синхронизация для
% каждого CA-кода)
NumCA2NextSync = Params.P20_NonCohTrackSatsAndBitSync.NumCA2NextSync;

% Половина количества дополнительных периодов CA-кода, используемых для
% синхронизации по времени
HalfNumCA4Sync = Params.P20_NonCohTrackSatsAndBitSync.HalfNumCA4Sync;

% Количество учитываемых значений задержки/набега синхронизации по
% времени
HalfCorLen = Params.P20_NonCohTrackSatsAndBitSync.HalfCorLen;

% Период, с которым производится отображение числа обработанных
% CA-кодов
NumCA2Disp = Params.P20_NonCohTrackSatsAndBitSync.NumCA2Disp;

% Максимальное число обрабатываемых CA-кодов (inf - до конца файла!)
MaxNumCA2Process = Params.P20_NonCohTrackSatsAndBitSync.MaxNumCA2Process;

% Количество бит, используемых для битовой синхронизации
NBits4Sync = Params.P20_NonCohTrackSatsAndBitSync.NBits4Sync;

%% СОХРАНЕНИЕ ПАРАМЕТРОВ
Track.NumCA2NextSync   = NumCA2NextSync;
Track.HalfNumCA4Sync   = HalfNumCA4Sync;
Track.HalfCorLen       = HalfCorLen;
Track.MaxNumCA2Process = MaxNumCA2Process;

BitSync.NBits4Sync     = NBits4Sync;

%% РАСЧЁТ ПАРАМЕТРОВ
% Длина CA-кода с учётом частоты дискретизации
CALen = 1023 * Res.File.R;

% Количество периодов CA-кода, приходящихся на один бит
CAPerBit = 20;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ТРЕКИНГ
[~,File] = ReadSignalFromFile(Res.File, 0, 0);
all_period_number = floor(File.SamplesLen/CALen);
dt=10^-3/CALen;
Shifts_EPL = -HalfCorLen:HalfCorLen; % EPL -- early; promt; late 

for SatNum  = 1:Res.Search.NumSats
    dNCO =0; % NCO - numerical control oscillator - ГУН цифровой
    Time_shifts_nacop = 0;
    CACode = GenCACode(Res.Search.SatNums(SatNum),1);
    CACode(CACode==0)=-1;
    CACode = repmat(CACode,Res.File.R,1);
    CACode = reshape(CACode,1,[]);
    dt=10^-3/CALen;
    
    Shift = Res.Search.SamplesShifts(SatNum);
    dNCO = exp(-1j*2*pi*Res.Search.FreqShifts(SatNum)*dt);
    Init_NCO = 1;
    for k = 1:all_period_number
        Signal = [];
        [Signal, ~] = ReadSignalFromFile(Res.File, Shift, CALen);
        Signal = Signal.*Init_NCO.*dNCO.^(0:CALen-1);
        if ~isnan(Signal)
            
            Track.CorVals{SatNum}(k) =  Signal*CACode';
            Track.SamplesShifts{SatNum}(k) =  Shift;
            
            Time_shifts(k) = 0;

            if mod(k,NumCA2NextSync +1)==0 % шаг точечной оценки - 100; заходим сюда в каждую точку
                
                EPL_sum = [];
                ii = 1;
                
                for sh = Shifts_EPL
                    EPL = [];
                    [EPL, ~] = ReadSignalFromFile(Res.File, Shift-HalfNumCA4Sync*CALen+sh, 2*HalfNumCA4Sync*CALen);
                    EPL = EPL.*dNCO.^(0:2*HalfNumCA4Sync*CALen-1).*Init_NCO*dNCO^sh;
                    EPL_sum(ii) = sum(abs(CACode*reshape(EPL,CALen,2*HalfNumCA4Sync)));
                    ii = ii+1;
                end
 
                [~,max_index]=max(EPL_sum) ;
                Shift = Shift + Shifts_EPL(max_index);
                % Time_shifts(k) - это маленький шифт на каждом k-ом шаге,
                % Time_shifts_nacop(k) - это накопл енный шифт (это и есть кумулятивная прямая кривая)
                %  Еще есть другой способ нахождения кумулятивной прямой (в тетради)
                Time_shifts(k) = Shifts_EPL(max_index);
            end
            
            if Time_shifts(k)~=0
                %при перемножении степени складываются!
                % Init нужен, чтобы сдвигать на time shift; когда его нет,
                % то в цикл не заходит
                Init_NCO = Init_NCO*dNCO^Time_shifts(k);
            end
            
            if k~=1
                Time_shifts_nacop(k) = Time_shifts_nacop(k-1) + Time_shifts(k);
            end
            Shift = Shift + CALen;
            
            if mod(k,NumCA2Disp) == 0
                disp(strcat('Обработано: ',num2str(k),' периодов CA кода'));
            end
        end
    end
    
    %% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - БИТОВАЯ СИНХРОНИЗАЦИЯ
    %dCA_all{SatNum} = Track.CorVals{SatNum}(1:end-1).*conj(Track.CorVals{SatNum}(2:end));
    for nn = 0:NBits4Sync
        % ищем разницу 0/pi
        % и сразу идут строчки по nn, которые потом складываем 
        dCA{SatNum}(nn+1,:) = Track.CorVals{SatNum}(nn*20+(2:20)).*conj(Track.CorVals{SatNum}(nn*20+(1:20-1)));
    end
    sum_dCA{SatNum} = abs(sum(dCA{SatNum})); % складываем корреляции, но не фазы а уже комплексные числа и где стык - будет минимум
                                             % 1*exp(0) = 1; 1*exp(pi) = -1
    BitSync.Cors(SatNum,:)=sum_dCA{SatNum};
    [~,BitSync.CAShifts(SatNum)] = min( sum_dCA{SatNum} );
    
    
    figure ();
    plot(abs(Track.CorVals{SatNum}),'.');
    title(num2str(Res.Search.SatNums(SatNum)));
    
    
    figure ();
    subplot(3,1,1);
    plot(sum_dCA{SatNum});
    title(num2str(Res.Search.SatNums(SatNum)));
    subplot(3,1,2);
    plot(angle(Track.CorVals{SatNum})/pi,'.');
    subplot(3,1,3);
    plot(Time_shifts_nacop); 
    
end
% Добавим новое поле с результатами в Res
Res.Track = Track;
Res.BitSync = BitSync;
end



