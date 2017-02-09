
f = @(x) x + 3*cos(x/10);
h = @(x) x^3;
a = @(x) 1 - 3/10*sin(x/10);
b = 1;
c = @(x) 3*x^2;

N = 300; P = 1; Q = 1; R = 100;

v = zeros(N,1); w = zeros(N,1);
x = zeros(N,1); y = zeros(N,1); xhat = zeros(N,1);

x(1) = 10;
y(1) = h(x(1));
xhat(1) = x(1) + 1;

figure(1),clf

for k = 2:N
    v(k,1) = randn(1)*sqrtm(Q);
    w(k,1) = randn(1)*sqrtm(R);
    x(k) = f(x(k-1)) + b*v(k-1);
    y(k) = h(x(k)) + w(k);
    [xhat(k,:),P] = ekf(f,h,a,b,c,Q,R,y(k),xhat(k-1,:),P);
    
end

subplot(2,1,1)
    plot(1:N,y,'k')
    xlabel('k'),ylabel('y')
    subplot(2,1,2) 
    plot(1:N,x,'r:',1:N, xhat,'b-')
    xlabel('k'),ylabel('x')
    legend('true','estimate','Location','SouthEast')