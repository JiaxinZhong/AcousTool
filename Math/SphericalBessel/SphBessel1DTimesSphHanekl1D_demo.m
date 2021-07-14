max_degree = 2;

zj = 1e1+1e-4*1i;
zh = linspace(1e2,1e3,200)+1e-4*1i;

res = SphBessel1DTimesSphHankel1D(max_degree, zj, zh);
res_exact = SphBessel1DTimesSphHankel1D_Exact(max_degree, zj, zh);

fig = Figure;
for i = 1:max_degree+1
    subplot(max_degree+1, 1, i);
    plot(abs(zh), abs(res_exact(i,:)));
    hold on
    plot(abs(zh), abs(res(i,:)), 'o');
end
