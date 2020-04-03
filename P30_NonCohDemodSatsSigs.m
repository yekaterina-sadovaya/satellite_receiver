function Res = P30_NonCohDemodSatsSigs(inRes, Params)
%
% Функция некогерентной демодуляции сигналов спутников
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
    Demod = struct( ...
        'Bits', {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждый элемент cell-массива Bits - массив 4хN значений 0/1, где
    %   N - количество демодулированных бит.

%% УСТАНОВКА ПАРАМЕТРОВ

%% РАСЧЁТ ПАРАМЕТРОВ
    % Количество периодов CA-кода, приходящихся на один бит
        CAPerBit = 20;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
i=1;

      for k = 1:Res.Search.NumSats
%       for k = 11
          sequence_1 = [];
          sequence_2 = [];
          sequence_3 = [];
          sequence_4 = [];
          sequence_5 = [];
          sequence_6 = [];
          
        dCA_all{i} = conj(Res.Track.CorVals{i}(1:end-1)).*Res.Track.CorVals{i}(2:end); % это между соседними отсчетами для графиков
        ddCA_all{i} = conj(dCA_all{i}(1:end-1)).*dCA_all{i}(2:end);
        
        dCA_all_20{i} = conj(Res.Track.CorVals{i}((Res.BitSync.CAShifts(k,1)+1):end-CAPerBit)).*Res.Track.CorVals{i}((Res.BitSync.CAShifts(k,1)+CAPerBit+1):end);
        
        % получили облака - суммировали:
        dCA_all_1{i} = sum(reshape(dCA_all_20{i}(1:floor(length(dCA_all_20{i})/CAPerBit)*CAPerBit),CAPerBit,floor(length(dCA_all_20{i})/CAPerBit)));
        ddCA_all_1{i}=conj(dCA_all_1{i}(1:end-1)).*dCA_all_1{i}(2:end); % устраняем зависимость по частоте
        
        sequence_1 = zeros(1,length(ddCA_all_1{i}));
        sequence_1(real(ddCA_all_1{i})<0) = 1;
        sequence_1(real(ddCA_all_1{i})>0) = 0;
        sequence_2 = sequence_1;
        sequence_1 = [1,sequence_1];
        sequence_2 = [0,sequence_2];
        sequence_3 = [1,mod(cumsum(sequence_1),2)];
        sequence_4 = [0,mod(cumsum(sequence_1),2)];
        sequence_5 = [1,mod(cumsum(sequence_2),2)];
        sequence_6 = [0,mod(cumsum(sequence_2),2)];
        Demod.Bits{i}(1,:) = mod(cumsum(sequence_3),2);
        Demod.Bits{i}(2,:) = mod(cumsum(sequence_4),2);
        Demod.Bits{i}(3,:) = mod(cumsum(sequence_5),2);
        Demod.Bits{i}(4,:) = mod(cumsum(sequence_6),2);
%         Demod.Bits{i}=zeros(1,length(demodulated_1{i}));
%         Demod.Bits{i}(find(real(d_demodulated_1{i})<0)) = 0;
%         Demod.Bits{i}(find(real(d_demodulated_1{i})>0)) = 1;
        
            figure ();
            subplot(5,1,1);
            plot(angle(dCA_all{i})/pi,'.');
            title(strcat(num2str(Res.Search.SatNums(k)),', Фазы после первого дифф.'));
            subplot(5,1,2);
            plot(angle(ddCA_all{i})/pi,'.');
            title(strcat(num2str(Res.Search.SatNums(k)),', Фазы после второго дифф.'));
            subplot(5,1,3);
            plot(angle(dCA_all_20{i})/pi,'.');
            title(strcat(num2str(Res.Search.SatNums(k)),', Фазы через 20, первое дифф.'));
            subplot(5,1,4);
            plot(angle(dCA_all_1{i})/pi,'.');
            title(strcat(num2str(Res.Search.SatNums(k)),', Средняя фаза бита, первое дифф.'));
            subplot(5,1,5);
            plot(angle(ddCA_all_1{i})/pi,'.');
            title(strcat(num2str(Res.Search.SatNums(k)),', Средняя фаза бита, второе дифф.'));
            
        i  = i+1;
    end
     Res.Demod = Demod;
end