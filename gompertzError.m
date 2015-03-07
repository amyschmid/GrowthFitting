function e=gompertzError(x,t,m)
b=x(1);
c=x(2);
h=exp(-exp(b-c.*t));
e=m'*m-(h'*m)^2/(h'*h);