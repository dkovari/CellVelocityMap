figure();
x = linspace(-2,2,100);
as = asinh(x);
at = tanh(0.8*x);
plot(x,x,'-k');hold on;
plot(x,as,'-r');
plot(x,at,'-b');