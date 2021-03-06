function [xhat_new, P_new, G] = ekf(f,h,A,B,C,Q,R,y,xhat,P)
   xhat = xhat(:);
   y = y(:);
   xhatm = f(xhat);
   Pm = A(xhatm)*P*A(xhatm)' + B*Q*B';
   G = Pm*C(xhatm)/(C(xhatm)'*Pm*C(xhatm)+R);
   xhat_new = xhatm + G*(y-h(xhatm));
   P_new = (eye(size(A(xhatm))) - G*C(xhatm)')*Pm;
end

