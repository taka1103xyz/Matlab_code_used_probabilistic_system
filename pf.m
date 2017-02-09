function [x_new,p_new] = pf(f,Q,y,p,M,sigma)
% f=A h=c B=b Q=noize R=noize y=Observation M=particle 
% sigma=particle dispersion
% 記憶領域の確保
p_new = zeros(1,M);
% 重み
weight = (1.0/M)*ones(1,M);
% サンプリング
p = f*p + sqrtm(Q)*randn(1,M);
% 尤度の計算
for i=1:M
    weight(1,i) = 1.0/(sqrt(2.0*pi)* sigma) * exp(-((p(1,i)-y)^2/(2.0*sigma^2)));
end
% 重みの正規化
weight_sum = sum(weight);
weight = weight/weight_sum;
% リサンプリング
for i = 1 : M
   p_new(1,i) = p(find(rand(1) <= cumsum(weight),1));
end
% 重み付き平均を算出
x_new = sum(p_new) / M;
end
