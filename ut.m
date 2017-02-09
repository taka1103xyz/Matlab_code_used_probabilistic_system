function [ ym, Pyy, Pxy ] = ut( f,xm,Pxx )
% ��x�N�g���ɐ��`
xm = xm(:);
% mapcols(f,x): x�̊e���f�Ŏʑ�����֐�
mapcols = @(f,x) cell2mat(cellfun(f,mat2cell(x,size(x,1),ones(1,size(x,2))),'UniformOutput',false));
% �萔                    
n = length(xm);
kappa = 3-n;
w0 = kappa/(n+kappa);
wi = 1/(2*(n+kappa));
W = diag([w0;wi*ones(2*n,1)]);
%% U�ϊ�
% �V�O�}�|�C���g�̐���
L = chol(Pxx);
X = [xm';
    ones(n,1)*xm'+sqrt(n+kappa)*L;
    ones(n,1)*xm'-sqrt(n+kappa)*L];
% �V�O�}�|�C���g�ɑΉ����� y ���v�Z
Y = mapcols(f,X')';
% y �̊��Ғl
ym = sum(W*Y)';
% �����U�s��
Yd = bsxfun(@minus,Y,ym');
Xd = bsxfun(@minus,X,xm');
Pyy = Yd'*W*Yd;
Pxy = Xd'*W*Yd;
end

