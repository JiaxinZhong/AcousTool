max_order = 2;
z = linspace(.1, 100, 1e4) + 1e-4*1i;

jd = SphBessel1D(max_order, z);
jd_exact = SphBessel1D_Exact(max_order, z);

fig = Figure;
for i = 1:max_order+1
    subplot(max_order+1, 1, i);
    plot(abs(z), abs(jd_exact(i,:)));
    hold on
    plot(abs(z), abs(jd(i,:)),'-.');
end
