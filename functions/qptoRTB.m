function x = qptoRTB(x,V,T1,C,a,mu)
x = V*((T1*C*x)+[a-mu;0;0;0;a-mu;0]);
end