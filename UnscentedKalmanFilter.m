%% ���ݒ�
% �����p�����[�^�l�̐ݒ�
rho = 1.23; g = 9.81; eta = 6e3;
M = 3e4; a = 3e4;
% �ϑ�����
T = 0.5;
EndTime = 30;
time = 0:T:EndTime;
N = EndTime/T+1;
n = 3;
R = 4e3;
Q = 0;
B = [0;0;0];
% �V�X�e��
f = @(x) [x(1)+T*x(2);
    x(2)+T*(0.5*rho*exp(-x(1)/eta)*x(2)^2*x(3)-g);
    x(3)];
h = @(x) sqrt(M^2+(x(1)-a).^2);
% f �̃��R�r�A�� (EKF�ŕK�v)
A = @(x) [1 T 0;
    -0.5*T*rho/eta*exp(-x(1)/eta)*x(2)^2*x(3) ...
    1+T*rho*exp(-x(1)/eta)*x(2)*x(3) ...
    0.5*T*rho*exp(-x(1)/eta)*x(2)^2;
    0 0 1];
% h �̃��R�r�A��(EKF�ŕK�v)
C = @(x) [(x(1)-a)/sqrt(M^2+(x(1)-a)^2); 0; 0];
%%�ϑ��f�[�^�̐���
w = randn(N,1)*sqrtm(R);
% �L���̈�̊m��
x = zeros(N,n); y = zeros(n,1);
% �������
x(1,:) = [90000;-6000;0.003]';
y(1) = h(x(1,:))+w(1);
% ���ԍX�V
for k = 2:N
    x(k,:) = f(x(k-1,:));
    y(k) = h(x(k,:))+w(k);
end
%% EKF�A���S���Y��
% ����l�L���̈�̊m��
xhat_ekf=zeros(N,3);
% �����l�̐ݒ�
xhat_ekf(1,:) = x(1,:);             % ��Ԑ���l
P_ekf = [9e3 0 0;0 4e5 0;0 0 0.4];  % �덷�����U
%����l�̔����X�V
for k = 2:N
    [xhat_ekf(k,:),P_ekf] = ...
        ekf(f,h,A,B,C,Q,R,y(k),xhat_ekf(k-1,:),P_ekf);
end
%% UKF�A���S���Y��
xhat_ukf=zeros(N,3);
% �����l
xhat_ukf(1,:) = x(1,:);
P_ukf = [9e3 0 0;0 4e5 0; 0 0 0.4];
% ����l�̔����X�V
for k = 2:N
    [xhat_ukf(k,:), P_ukf] = ...
        ukf(f,h,B,Q,R,y(k),xhat_ukf(k-1,:),P_ukf);
end
%% ���ʂ̕\��
% ��Ԃ̐��l�̐}��
figure, clf
for p=1:3
    subplot(3,1,p)
    plot(time,x(:,p));
    xlabel('Time [s]'), ylabel(sprintf('x%d',p))
end
% ����l�̐}��
figure(2),clf
for p=1:3
    subplot(3,1,p)
    plot(time,x(:,p),'k', ...
        time,xhat_ekf(:,p),'b-', ...
        time,xhat_ukf(:,p),'r:');
    xlabel('Time [s]'),ylabel(sprintf('x%d',p))
    legend('true','ekf','ukf')
end
% RMSE �̕\��
RMSE = @(x) sqrt(mean(x.^2));
fprintf('%10s %10s %10s\n','variable','RMSE(ekf)','RMSE(ukf)');
for p=1:3
    vname = sprintf('x%d', p);
    fprintf('%10s % 10.5f % 10.5f\n', ...
        vname, RMSE(xhat_ekf(:,p)-x(:,p)), ...
        RMSE(xhat_ukf(:,p)-x(:,p)));
end
