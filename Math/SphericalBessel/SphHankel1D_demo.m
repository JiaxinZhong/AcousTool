max_order = 2;
z = linspace(.1, 10, 1e4) + 1e-4*1i;

hd = SphHankel1D(max_order, z);
hd_exact = SphHankel1D_Exact(max_order, z);

fig = Figure;
for i = 1:max_order+1
    subplot(max_order+1, 1, i);
    plot(abs(z), abs(hd_exact(i,:)));
    hold on
    plot(abs(z), abs(hd(i,:)),'-.');
end
