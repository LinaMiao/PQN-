% This is for the question 9.10 from bv_cvxhook
% pure Newon method: Newton?s method with fixed step size t = 1 can diverge
% if the initial point is not close to x?. consider two examples.

% (a) f(x) = log(ex + e?x) has a unique minimizer x? = 0. Run Newton?s 
% method with fixed step size t = 1, starting at x(0) = 1 and at x(0) = 1.1.

x = -2:0.1:2;
y = @(x) log(exp(x)+exp(-x));
g = @(x) (exp(x)-exp(-x))./(exp(x)+exp(-x));
h = @(x) ((exp(x)+exp(-x)).^2-(exp(x)-exp(-x)).^2)./(exp(x)+exp(-x)).^2;
tol = 10^-12;
for x0 = [1 1.1]
    itn = 0;
    while abs(g(x0)) >= tol && itn<=10;
        x0 = x0 - h(x0)\g(x0);
        itn = itn+1;
        fprintf('%2.0f  %3.2e  %3.2e\n',itn,g(x0),x0);
    end
end
figure;plot(x,y(x));hold on; plot(x,g(x),'r');legend('f','divf');

% (b) f (x) = ? log x + x has a unique minimizer x? = 1. Run Newton?s method
% with fixed step size t = 1, starting at x(0) = 3.
clear;
x = -2:0.1:2;
y = @(x) -log(x) + x;
g = @(x) -1./x + 1;
h = @(x) 1/x^2;
tol = 10^-12;
for x0 = [3]
    itn = 0;
    while abs(g(x0)) >= tol && itn<=10;
        x0 = x0 - h(x0)\g(x0);
        itn = itn+1;
        fprintf('%2.0f  %3.2e  %3.2e\n',itn,g(x0),x0);
    end
end
figure;plot(x,y(x));hold on; plot(x,g(x),'r');legend('f','divf');
