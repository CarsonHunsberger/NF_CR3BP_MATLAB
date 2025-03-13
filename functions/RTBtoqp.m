function x = RTBtoqp(x,Cinv,T1inv,Vinv,a,mu)

x = Cinv*T1inv*Vinv*x - Cinv*T1inv*[a-mu;0;0;0;a-mu;0];

end