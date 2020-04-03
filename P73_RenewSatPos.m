function OutSatPos = P73_RenewSatPos(SatPos, TProp, Params) %#ok<INUSD>
%
% ѕересчЄт координат спутника по новому значению времени распространени€
% сигнала
%
% ¬ходные переменные
%   SatPos - массив координат спутника (8х1) из P72_GetSatPos;
%   TProp - врем€ распространени€ сигнала. оно же t_t
%
% ¬ыходные переменные
%   OutSatPos - массив (3x1) скорректированных координат спутника.

% WGS 84 value of the earth's rotation rate (rad/sec)
     OutSatPos = struct('x',[]);
     OMEGA_e = 7.2921151467*10^-5;
     OMEGA = SatPos.Omega_k - OMEGA_e*TProp ;
    
     OutSatPos.x = SatPos.xs_k *cos(OMEGA) - SatPos.ys_k *cos(SatPos.i_k)*sin(OMEGA);
     OutSatPos.y = SatPos.xs_k*sin(OMEGA) + SatPos.ys_k*cos(SatPos.i_k)*cos(OMEGA);
     OutSatPos.z = SatPos.ys_k*sin(SatPos.i_k);
% % 
% %      OMEGA_e = 7.2921151467*10^-5;
% %      OMEGA = SatPos(7)- OMEGA_e*TProp ;
% %     
% %      OutSatPos(1) = SatPos(4)*cos(OMEGA) - SatPos(5) *cos(SatPos(6))*sin(OMEGA);
% %      OutSatPos(2) = SatPos(4)*sin(OMEGA) + SatPos(5)*cos(SatPos(6))*cos(OMEGA);
% %      OutSatPos(3) = SatPos(5)*sin(SatPos(6));
     
end