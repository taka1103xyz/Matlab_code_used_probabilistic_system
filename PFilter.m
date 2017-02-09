% パラメータ設定
A = 1;
b = 1;
c = 1;
P = 0;
Q = 1;
R = 2;
dt = 1;                 % 時間間隔
M = 100;                % パーティクル数
N = 100;                % 計測数
sigma = 2;            

% 行列の定義
x=zeros(N,1);
y=zeros(N,1);
xkhat=zeros(N,1);
xphat=zeros(N,1);

% ノイズの生成
v = randn*sqrtm(Q);
w = randn*sqrtm(R);

% パーティクル行列
p = zeros(1,M);
p_bar = zeros(2,M);

% 重み
weight = (1.0/M)*ones(1,M);

% グラフの作成
figure(1),clf;

% 初期値設定
y(1,1)=c'*x(1,1)'+w(1);

% カルマンフィルタ
for k=1:N
    if(k>1)
        % 新たなノイズの生成   
        v = randn*sqrtm(Q);
        w = randn*sqrtm(R);
 
        % 真値・観測地・推定値の代入
        x(k,1)=A*x(k-1,1)'+b*v;
        y(k,1)=c'*x(k,1)'+w;
        [xkhat(k,1),P,G]=kf(A,b,0,c,Q,R,0,y(k,1),xkhat(k-1,1),P);
    end
 
    % こっからパーティクル
 
    % サンプリング
    p = A*p + b*sqrtm(Q)*randn(1,M);
 
    % 尤度の計算
    for i=1:M
        weight(1,i) = 1.0/(sqrt(2.0*pi)* sigma) * exp(-((p(1,i)-y(k,1))^2/(2.0*sigma^2)));
    end
 
    % 重みを正規化
    weight_sum = sum(weight);
    weight = weight/weight_sum;
 
    % p_barに入れておく
    p_bar = [p;weight];
     
    % リサンプリング
    for i = 1 : M
        p(i) = p(find(rand(1) <= cumsum(weight),1));
    end
    
    % 重み付き平均を算出
    xphat(k,1) = sum(p) / M;
 
    if(k>1)
        % プロット
        plot(1:N,y,'k:',1:N,x,'r--',1:N,xkhat,'b--',1:N,xphat,'g--');
        % pause;
    end
end

% ラベルを貼る
xlabel('No. of samples');
legend('measured','true','estimates kalman','estimate particle');