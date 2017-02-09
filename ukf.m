function [xhat_new, P_new, G] = ukf(f,h,B,Q,R,y,xhat,P)
% 列ベクトルに整形
sprintf('a');
xhat = xhat(:);
y = y(:);
% 事前推定値
[xhatm,Pm] = ut(f,xhat,P);
Pm         = Pm + B*Q*B';
[yhatm,Pyy,Pxy] = ut(h,xhatm,Pm);
% カルマンゲイン行列
G = Pxy/(Pyy+R);
% 事後推定値
xhat_new = xhatm + G*(y-yhatm); % 状態
P_new    = Pm - G*Pxy';         % 誤差共分散
end
