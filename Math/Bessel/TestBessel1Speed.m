n = (0:2e3).';
x = linspace(10,100+1i,2e3);

tic
for i = 1:length(n)
    besselj(n(i),x);
end
toc

tic
for i = 1:length(x)
    besselj(n, x(i));
end
toc;