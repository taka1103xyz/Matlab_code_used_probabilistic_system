% �p�����[�^�ݒ�
A=1;
b=1;
c=1;
P=0;
Q=0.01;
R=0.02;
N=300;

% �s��̒�`
x=zeros(N,3);
y=zeros(N,3);
xhat=zeros(N,3);

% �m�C�Y�̐���
v = randn(1,3)*sqrtm(Q);
w = randn(1,3)*sqrtm(R);

% �O���t�̍쐬
figure(1),clf;

% �����l�ݒ�
y(1,:) = c*x(1,:)+w;

% �J���}���t�B���^
for k=2:N
 % �V���ȃm�C�Y�̐���
 t=k*pi*0.1;
 u=[sin(t),cos(t),k];
 v = randn(1,3)*sqrtm(Q);
 w = randn(1,3)*sqrtm(R);
 
 % �^�l�E�ϑ��n�E����l�̑��
 x(k,:)=b*u+v; %A*x(k-1,:)++v
 y(k,:)=c*x(k,:)+w; %+w
 [xhat(k,:),P,G]=kf(A,b,0,c,Q,R,0,y(k,:),xhat(k-1,:),P);

 % �v���b�g
 %plot(1:N,y,'k:',1:N,x,'r--',1:N,xhat,'b--');
 
 %pause;
end
% �v���b�g
scatter3(x(:,1),x(:,2),x(:,3),'k');
hold on;
scatter3(xhat(:,1),xhat(:,2),xhat(:,3),'h');
% ���x����\��
%xlabel('No. of samples');
%legend('measured','true','estimates');