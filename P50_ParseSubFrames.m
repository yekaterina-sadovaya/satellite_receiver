function Res = P50_ParseSubFrames(inRes, Params)
% 1) nan
%2) сосчитать параметры

% Функция демодуляции сигналов спутников
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
SatsData = struct( ...
    'isSat2Use', zeros(1, Res.Search.NumSats), ...
    'TLM', {cell(Res.Search.NumSats, 1)}, ...
    'HOW', {cell(Res.Search.NumSats, 1)}, ...
    'SF1', {cell(Res.Search.NumSats, 1)}, ...
    'SF2', {cell(Res.Search.NumSats, 1)}, ...
    'SF3', {cell(Res.Search.NumSats, 1)}, ...
    'SF4', {cell(Res.Search.NumSats, 1)}, ...
    'SF5', {cell(Res.Search.NumSats, 1)} ...
    );
% Элементами всех cell-массивов (TLM, HOW, SF1, SF2, SF3, SF4, SF5)
% являются структуры-массивы (1хN) с результатами парсинга, где N -
% количество обработанных для спутника подкадров. Если какое то поле не
% расшифровано из-за того, что не сошлось CRC, то его значение должно
% быть установлено в nan. isSat2Use - массив флагов, указывающих,
% было ли расшифровано хотя бы одно поле HOW.TOW_Count_Message, т.е.
% имеет ли смысл в дальнейшем изучать содержимое подкадров (конечно,
% isSat2Use = 0, если у этого спутника isSubFrameSync = 0).

%% УСТАНОВКА ПАРАМЕТРОВ
SatNums = Res.Search.SatNums;
NumSats = Res.Search.NumSats;
isSubFrameSync = Res.SubFrames.isSubFrameSync;
Words = Res.SubFrames.Words;
%% РАСЧЁТ ПАРАМЕТРОВ

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ С УСПЕШНОЙ
% ПОДКАДРОВОЙ СИНХРОНИЗАЦИЕЙ
for k = 1:Res.Search.NumSats
    count = 0;
    if isSubFrameSync(k)
        for kk =1:size(Words{k},1)
            SatsData.TLM{k}(kk) = ParseTLM(Words{k}{kk,1});
            SatsData.HOW{k}(kk) = ParseHOW(Words{k}{kk,2});

            Fun = str2func(strcat('ParseSF',num2str(SatsData.HOW{k}(kk).Subframe_ID)));
            
            Var  = Fun(Words{k}(kk,:));
            
            % создали структуру SatsData.SF' и подгрузили в неё Var в
            % зависимости от ID
            eval(['SatsData.SF',num2str(SatsData.HOW{k}(kk).Subframe_ID),'{k}(kk)',' = Var;']);
            if length(SatsData.HOW{k}(kk).TOW)~=0
               count = count+1;
            end
        end
    else
        SatsData.isSat2Use(k)=0;
    end
    if count>0
       SatsData.isSat2Use(k)=1;
    end
end
Res.SatsData = SatsData;
end

function [Bits, isCRC] = Words2BitFrame(SubFrames)
% Из (1х8) cell-массива Words составим кадр, т.е. добавим нулевые биты CRC
% и нулевые первые два слова. Это удобно для анализа кода по спецификации.
% Также составим массив флагов, указывающих на то, сошлось CRC в конкретном
% слове или нет
end

function Data = ParseSF1(Words)
%p.103
% Парсинг подкадра №1
    for i = 3:10
        if i ==3
            Data = struct( ...
                'WeekNumber',       [bi2de(fliplr(Words{i}(1:10)))], ...
                ...
                'CodesOnL2',        [bi2de(fliplr(Words{i}(11:12)))], ...
                ...
                'URA',       [bi2de(fliplr(Words{i}(13:16)))], ...
                ...
                'URA_in_meters',       [], ...
                ...
                'SV_Health_Summary', [Words{i}(17)], ...
                ...
                'SV_Health',       [bi2de(fliplr(Words{i}(18:22)))], ...
                ...
                'IODC_MSB',       [Words{i}(23:24)] );
            switch Data.CodesOnL2
                case 0
                    Data.CodesOnL2 = 'Reserved';
                case 1
                    Data.CodesOnL2 = 'P code ON';
                case 2
                    Data.CodesOnL2 = 'C/A code ON';
            end

            switch Data.URA
                case 0
                    Data.URA_in_meters = '0.00 < URA <= 2.40';
                case 1
                    Data.URA_in_meters = '2.40 < URA <= 3.40';
                case 2
                    Data.URA_in_meters = '3.40 < URA <= 4.85';
                case 3
                    Data.URA_in_meters = '4.85 < URA <= 6.85';
                case 4
                    Data.URA_in_meters = '6.85 < URA <= 9.65';
                case 5
                    Data.URA_in_meters = '9.65 < URA <= 13.65';
                case 6
                    Data.URA_in_meters = '13.65 < URA <= 24.00';
                case 7
                    Data.URA_in_meters = '24.00 < URA <= 48.00';
                case 8
                    Data.URA_in_meters = '48.00 < URA <= 96.00';
                case 9
                    Data.URA_in_meters = '96.00 < URA <= 192.00';
                case 10
                    Data.URA_in_meters = '192.00 < URA <= 384.00';
                case 11
                    Data.URA_in_meters = '384.00 < URA <= 768.00';
                case 12
                    Data.URA_in_meters = '768.00 < URA <= 1536.00';
                case 13
                    Data.URA_in_meters = '1536.00 < URA <= 3072.00';
                case 14
                    Data.URA_in_meters = '3072.00 < URA <= 6144.00';
                case 15
                    Data.URA_in_meters = '6144.00 < URA <= or no accuracy prediction is available - standard positioning service users are advised to use the SV at their own risk.';
            end
            
            switch Data.SV_Health_Summary
                case 0
                    Data.SV_Health_Summary = 'all NAV data are OK';
                case 1
                    Data.SV_Health_Summary = 'some or all NAV data are bad';
            end
            
            switch Data.SV_Health
                case 0
                    Data.SV_Health = 'All Signals OK';
                case 1
                    Data.SV_Health = 'All Signals Weak';
                case 2
                    Data.SV_Health = 'All Signals Dead';
                case 3
                    Data.SV_Health = 'All Signals Have No Data Modulation';
                case 4
                    Data.SV_Health = 'L1 P Signal Weak';
                case 5
                    Data.SV_Health = 'L1 P Signal Dead';
                case 6
                    Data.SV_Health = 'L1 P Signal Has No Data Modulation';
                case 7
                    Data.SV_Health = 'L2 P Signal Weak';
                case 8
                    Data.SV_Health = 'L2 P Signal Dead';
                case 9
                    Data.SV_Health = 'L2 P Signal Has No Data Modulation';
                case 10
                    Data.SV_Health = 'L1 C Signal Weak';
                case 11
                    Data.SV_Health = 'L1 C Signal Dead';
                case 12
                    Data.SV_Health = 'L1 C Signal Has No Data Modulation';
                case 13
                    Data.SV_Health = 'L2 C Signal Weak';
                case 14
                    Data.SV_Health = 'L2 C Signal Dead';
                case 15
                    Data.SV_Health = 'L2 C Signal Has No Data Modulation';
                case 16
                    Data.SV_Health = 'L1 & L2 P Signal Weak';
                case 17
                    Data.SV_Health = 'L1 & L2 P Signal Dead';
                case 18
                    Data.SV_Health = 'L1 & L2 P Signal Has No Data Modulation';
                case 19
                    Data.SV_Health = 'L1 & L2 C Signal Weak';
                case 20
                    Data.SV_Health = 'L1 & L2 C Signal Dead';
                case 21
                    Data.SV_Health = 'L1 & L2 C Signal Has No Data Modulation';
                case 22
                    Data.SV_Health = 'L1 Signal Weak';
                case 23
                    Data.SV_Health = 'L1 Signal Dead';
                case 24
                    Data.SV_Health = 'L1 Signal Has No Data Modulation';
                case 25
                    Data.SV_Health = 'L2 Signal Weak';
                case 26
                    Data.SV_Health = 'L2 Signal Dead';
                case 27
                    Data.SV_Health = 'L2 Signal Has No Data Modulation';
                case 28
                    Data.SV_Health = 'SV Is Temporarily Out (Do not use this SV during current pass)';
                case 29
                    Data.SV_Health = 'SV Will Be Temporarily Out (Use with caution)';
                case 30
                    Data.SV_Health = 'One Or More Signals Are Deformed, However The Relevant URA Parameters Are Valid';
                case 31
                    Data.SV_Health = 'More Than One Combination Would Be Required To Describe Anomalies';
            end
        end
        
        
        if i == 4
            Data.Flag_L2_P_Code = [Words{i}(1)] ;

            switch Data.Flag_L2_P_Code
                case 0
                    Data.Flag_L2_P_Code= 'The NAV data stream was commanded ON on the P-code of the L2 channel.';
                case 1
                    Data.Flag_L2_P_Code = 'The NAV data stream was commanded OFF on the P-code of the L2 channel.';
            end
            %Бит со 2-24 зарезервированы
        end
        
        
        if i == 5        
            %Бит со 1-24 зарезервированы
        end
        if i == 6
            %Бит со 1-24 зарезервированы
        end
        if i == 7

            Data.T_GD = 2^(-31)* comp2de(fliplr(Words{i}(17:24)));
            %Бит со 1-16 зарезервированы
        end
        if i == 8
            Data.IODC_LSB = Words{i}(1:8);
            Data.t_oc = 2^(4)*bi2de(fliplr(Words{i}(9:24)));
            Data.IODC = bi2de(fliplr([Data.IODC_MSB,Data.IODC_LSB]));
        end
        if i == 9
            Data.a_f2 = 2^(-55)* comp2de(fliplr(Words{i}(1:8)));
            Data.a_f1 = 2^(-43)* comp2de(fliplr(Words{i}(9:24)));
        end
        if i == 10
            Data.a_f0 =  2^-31*comp2de(fliplr(Words{i}(1:22)));
            %Бит со 23-24 неинформационные биты для вычисления пэрити
        end
        
    end
end

    function Data = ParseSF2(Words)
        %p.113
        % Парсинг подкадра №2
        for i = 3:10
            if i == 3
 
                Data = struct( ...
                    'IODE',       [bi2de(fliplr(Words{i}(1:8)))], ...
                    ...
                    'C_rs',       [2^-5*comp2de(fliplr(Words{i}(9:24)))] ...
                    ...
                    );
            end
            if i == 4
                Data.Delta_n = 2^-43*comp2de(fliplr(Words{i}(1:16)));
                Data.M_0_MSB = Words{i}(17:24);
            end
            if i == 5

                Data.M_0_LSB = Words{i}(1:24);
                Data.M_0 = 2^-31*comp2de(fliplr([Data.M_0_MSB,Data.M_0_LSB]));
            end
            if i == 6

                Data.C_uc = 2^-29*comp2de(fliplr(Words{i}(1:16)));
                Data.e_MSB = Words{i}(17:24);
 
            end
            if i == 7
                Data.e_LSB = Words{i}(1:24);
                Data.e = 2^-33*bi2de(fliplr([Data.e_MSB,Data.e_LSB]));
            end
            if i == 8
                Data.C_us = 2^-29*comp2de(fliplr(Words{i}(1:16)));
                Data.sqrt_A_MSB = Words{i}(17:24);

            end
            if i == 9
                Data.sqrt_A_LSB = Words{i}(1:24);
                Data.sqrtA = 2^-19*bi2de(fliplr([Data.sqrt_A_MSB,Data.sqrt_A_LSB]));
            end
            if i == 10
  
                Data.t_oe = 2^4*bi2de(fliplr(Words{i}(1:16)));
                Data.Fit_Interval_Flag = Words{i}(17);
                Data.AODO = 900* bi2de(fliplr(Words{i}(18:22)));

                switch Data.Fit_Interval_Flag
                    case 0
                        Data.Fit_Interval_Flag = '4 hours';
                    case 1
                        Data.Fit_Interval_Flag = 'greater than 4 hours';
                end
 
            end
            
        end
    end

    function Data = ParseSF3(Words)
        %
        % Парсинг подкадра №3
        for i = 3:10

            if i == 3

                Data = struct( ...
                    'C_ic',       [2^-29*comp2de(fliplr(Words{i}(1:16)))], ...
                    ...
                    'OMEGA_0_MSB',       [Words{i}(17:24)] ...
                    ...
                    );

            end
            if i == 4

                Data.OMEGA_0_LSB = Words{i}(1:24);
                Data.Omega_0 = 2^-31*comp2de(fliplr([Data.OMEGA_0_MSB,Data.OMEGA_0_LSB]));
            end
            if i == 5

                Data.C_is = 2^-29*comp2de(fliplr(Words{i}(1:16)));
                Data.i_0_MSB = Words{i}(17:24);
   
            end
            if i == 6
                Data.i_0_LSB = Words{i}(1:24);
                Data.i_0 = 2^-31*comp2de(fliplr([Data.i_0_MSB,Data.i_0_LSB]));
            end
            if i == 7
 
                Data.C_rc = 2^-5*comp2de(fliplr(Words{i}(1:16)));
                Data.omega_MSB = Words{i}(17:24);
 
            end
            if i == 8

                Data.omega_LSB = Words{i}(1:24);
                Data.omega = 2^-31*comp2de(fliplr([Data.omega_MSB,Data.omega_LSB]));
            end
            if i == 9

                Data.DOmega = 2^-43*comp2de(fliplr(Words{i}(1:24)));

            end
            if i == 10
                Data.IODE =bi2de(fliplr( Words{i}(1:8)));
                Data.IDOT = 2^-43*comp2de(fliplr( Words{i}(9:22)));
            end
        end
        
    end

function Data = ParseSF4(Words)
%
% Парсинг подкадра №4 - реализован только для (SV_Page_ID = 56)

                Data = struct( ...
                    'Almanac',  []);
end

function Data = ParseSF5(Words)
%
% Парсинг подкадра №5

                Data = struct( ...
                    'Almanac',  []);
% Парсинг не реализован

end

function Data = ParseTLM(Word)
%
% Парсинг слова TLM
Data = struct( ...
    'TLM_Message',       [Word(9:22)], ... 
    ...
    'TLM_Integrity_Status_Flag',       [Word(23)] ...
    );
%Бит 24 зарезервирован 
switch Data.TLM_Integrity_Status_Flag
    case 0
        Data.TLM_Integrity_Status_Flag = 'The legacy level of integrity assurance';
    case 1
        Data.TLM_Integrity_Status_Flag = 'An enhanced level of integrity assurance';
end
end

function Data = ParseHOW(Word)
%
% Парсинг слова HOW
Data = struct( ...
    'TOW',       [bi2de(fliplr(Word(1:17)))], ...%это TOW_MSB оно же HOW
    ...
    'Alert_Flag',       [Word(18)], ...
    ...
    'Anti_Spoof_Flag',       [Word(19)], ...
    ...
    'Subframe_ID',       [bi2de(fliplr(Word(20:22)))] ...
    ...
    );
switch Data.Alert_Flag
    case 0
        Data.Alert_Flag = 'No alert';
    case 1
        Data.Alert_Flag = 'The signal URA may be worse than indicated in subframe 1';
end
switch Data.Anti_Spoof_Flag
    case 0
        Data.Anti_Spoof_Flag = 'The A-S mode is OFF';
    case 1
        Data.Anti_Spoof_Flag = 'The A-S mode is ON';
end

% Бит с 23 по 24 Solved for bits to preserve parity check with zeros in
% bits 29 and 30???????
end

function Out = comp2de(In)
%Неотрицательное число такое же как в прямом коде, отризательное - делается
%сначала обратный код, т.е конвертируются все значения и вперед добавляется
%1, далее
% Функция перевода двоичного дополнительного кода в десятичное число

% доп. код -- наиболее распространённый способ представления отрицательных целых чисел в компьютерах. 
% Он позволяет заменить операцию вычитания на операцию сложения и сделать операции сложения и 
% вычитания одинаковыми для знаковых и беззнаковых чисел, 
% чем упрощает архитектуру ЭВМ.

%последнее число смотрим, потому что мы уже перевернули ранее
%правильную запись fliplr
if In(length(In))==0 % если последний 0 - положит.
   Out = bi2de(In(1:length(In)-1));
else
   % Крутить не надо потому что две операции и они взаимоуничтожились
   Out_obr = de2bi(bi2de(In(1:length(In)-1))-1);
   if length(In)-1~=length(Out_obr)
       Out_obr=[Out_obr,zeros(1,length(In)-1-length(Out_obr))];
   end
   Out_bi = zeros(1,length(Out_obr));
   Out_bi(Out_obr==0)=1;
   Out = -bi2de(Out_bi);
end

end
% function DE = Scale(BI,bits_number,scale_factor)
% %Неотрицательное число такое же как в прямом коде, отризательное - делается
% %сначала обратный код, т.е конвертируются все значения и вперед добавляется
% %1, далее
% % Функция перевода двоичного дополнительного кода в десятичное число
% BI(1:length(BI)-bits_number())
% 
% end