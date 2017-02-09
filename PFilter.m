% �p�����[�^�ݒ�
A = 1;
b = 1;
c = 1;
P = 0;
Q = 1;
R = 2;
dt = 1;                 % ���ԊԊu
M = 100;                % �p�[�e�B�N����
N = 100;                % �v����
sigma = 2;            

% �s��̒�`
x=zeros(N,1);
y=zeros(N,1);
xkhat=zeros(N,1);
xphat=zeros(N,1);

% �m�C�Y�̐���
v = randn*sqrtm(Q);
w = randn*sqrtm(R);

% �p�[�e�B�N���s��
p = zeros(1,M);
p_bar = zeros(2,M);

% �d��
weight = (1.0/M)*ones(1,M);

% �O���t�̍쐬
figure(1),clf;

% �����l�ݒ�
y(1,1)=c'*x(1,1)'+w(1);

% �J���}���t�B���^
for k=1:N
    if(k>1)
        % �V���ȃm�C�Y�̐���   
        v = randn*sqrtm(Q);
        w = randn*sqrtm(R);
 
        % �^�l�E�ϑ��n�E����l�̑��
        x(k,1)=A*x(k-1,1)'+b*v;
        y(k,1)=c'*x(k,1)'+w;
        [xkhat(k,1),P,G]=kf(A,b,0,c,Q,R,0,y(k,1),xkhat(k-1,1),P);
    end
 
    % ��������p�[�e�B�N��
 
    % �T���v�����O
    p = A*p + b*sqrtm(Q)*randn(1,M);
 
    % �ޓx�̌v�Z
    for i=1:M
        weight(1,i) = 1.0/(sqrt(2.0*pi)* sigma) * exp(-((p(1,i)-y(k,1))^2/(2.0*sigma^2)));
    end
 
    % �d�݂𐳋K��
    weight_sum = sum(weight);
    weight = weight/weight_sum;
 
    % p_bar�ɓ���Ă���
    p_bar = [p;weight];
     
    % ���T���v�����O
    for i = 1 : M
        p(i) = p(find(rand(1) <= cumsum(weight),1));
    end
    
    % �d�ݕt�����ς��Z�o
    xphat(k,1) = sum(p) / M;
 
    if(k>1)
        % �v���b�g
        plot(1:N,y,'k:',1:N,x,'r--',1:N,xkhat,'b--',1:N,xphat,'g--');
        % pause;
    end
end

% ���x����\��
xlabel('No. of samples');
legend('measured','true','estimates kalman','estimate particle');