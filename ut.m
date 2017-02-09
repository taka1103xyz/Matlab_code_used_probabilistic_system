function [ ym, Pyy, Pxy ] = ut( f,xm,Pxx )
% 列ベクトルに整形
xm = xm(:);
% mapcols(f,x): xの各列をfで写像する関数
mapcols = @(f,x) cell2mat(cellfun(f,mat2cell(x,size(x,1),ones(1,size(x,2))),'UniformOutput',false));
% 定数                    
n = length(xm);
kappa = 3-n;
w0 = kappa/(n+kappa);
wi = 1/(2*(n+kappa));
W = diag([w0;wi*ones(2*n,1)]);
%% U変換
% シグマポイントの生成
L = chol(Pxx);
X = [xm';
    ones(n,1)*xm'+sqrt(n+kappa)*L;
    ones(n,1)*xm'-sqrt(n+kappa)*L];
% シグマポイントに対応する y を計算
Y = mapcols(f,X')';
% y の期待値
ym = sum(W*Y)';
% 共分散行列
Yd = bsxfun(@minus,Y,ym');
Xd = bsxfun(@minus,X,xm');
Pyy = Yd'*W*Yd;
Pxy = Xd'*W*Yd;
end

