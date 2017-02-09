% パラメータ設定
A=1;
b=1;
c=1;
P=0;
Q=0.01;
R=0.02;
N=300;

% 行列の定義
x=zeros(N,3);
y=zeros(N,3);
xhat=zeros(N,3);

% ノイズの生成
v = randn(1,3)*sqrtm(Q);
w = randn(1,3)*sqrtm(R);

% グラフの作成
figure(1),clf;

% 初期値設定
y(1,:) = c*x(1,:)+w;

% カルマンフィルタ
for k=2:N
 % 新たなノイズの生成
 t=k*pi*0.1;
 u=[sin(t),cos(t),k];
 v = randn(1,3)*sqrtm(Q);
 w = randn(1,3)*sqrtm(R);
 
 % 真値・観測地・推定値の代入
 x(k,:)=b*u+v; %A*x(k-1,:)++v
 y(k,:)=c*x(k,:)+w; %+w
 [xhat(k,:),P,G]=kf(A,b,0,c,Q,R,0,y(k,:),xhat(k-1,:),P);

 % プロット
 %plot(1:N,y,'k:',1:N,x,'r--',1:N,xhat,'b--');
 
 %pause;
end
% プロット
scatter3(x(:,1),x(:,2),x(:,3),'k');
hold on;
scatter3(xhat(:,1),xhat(:,2),xhat(:,3),'h');
% ラベルを貼る
%xlabel('No. of samples');
%legend('measured','true','estimates');