function [x_new,p_new] = pf(f,Q,y,p,M,sigma)
% f=A h=c B=b Q=noize R=noize y=Observation M=particle 
% sigma=particle dispersion
% �L���̈�̊m��
p_new = zeros(1,M);
% �d��
weight = (1.0/M)*ones(1,M);
% �T���v�����O
p = f*p + sqrtm(Q)*randn(1,M);
% �ޓx�̌v�Z
for i=1:M
    weight(1,i) = 1.0/(sqrt(2.0*pi)* sigma) * exp(-((p(1,i)-y)^2/(2.0*sigma^2)));
end
% �d�݂̐��K��
weight_sum = sum(weight);
weight = weight/weight_sum;
% ���T���v�����O
for i = 1 : M
   p_new(1,i) = p(find(rand(1) <= cumsum(weight),1));
end
% �d�ݕt�����ς��Z�o
x_new = sum(p_new) / M;
end
