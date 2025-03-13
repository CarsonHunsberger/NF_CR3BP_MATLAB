function statedot = CR3BP(t,state,mu)

statedot(1,1) = state(4,1);
statedot(2,1) = state(5,1);
statedot(3,1) = state(6,1);

x = state(1); y = state(2); z = state(3);
xdot = state(4); ydot = state(5); %zdot = state(6);

r1 = sqrt((x+mu)^2+y^2+z^2); r2 = sqrt((x+mu-1)^2+y^2+z^2);

statedot(4,1) = 2*ydot+x-(1-mu)*(x+mu)/(r1^3)-mu*(x+mu-1)/(r2^3);
statedot(5,1) = -2*xdot+y-(1-mu)*y/(r1^3)-mu*y/(r2^3);
statedot(6,1) = -(1-mu)*z/(r1^3)-mu*z/(r2^3);
end