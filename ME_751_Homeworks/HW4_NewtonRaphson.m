clc;
clear all;

k = 0;
n = 1000; 

xk = zeros(n,1);

x0 = -1;

xk(1,1) = x0;
x = 1;
i = 2;

syms y;

while x > 1e-15
    g = (y^6)-y-1;
    dg = 6*(y^5)-1;
    f(y)=g;
    df(y) = dg; 
    
    xk(i) = xk(i-1)-(f(xk(i-1))/df(xk(i-1)));
    x = abs(xk(i)-xk(i-1));
    i = i + 1;
    
    if i==1000
        x = 0;
    end 
    
end

xk = xk(1:(i-1),:);

