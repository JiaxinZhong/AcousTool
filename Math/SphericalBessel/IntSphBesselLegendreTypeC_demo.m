max_order = 2;
kr = linspace(1,30)+1e-4*1i;
ka = 5+1e-5*1i;
kd = 50+1e-5*1i;
int_a = kd;
int_b = kd*1.1;

res = IntSphBesselLegendreTypeC(max_order, kr, ka, kd, ...
    int_a, int_b, 50);
res_exact = IntSphBesselLegendreTypeC_Exact(max_order, kr, ka, ...
    kd, int_a, int_b);

fig = Figure;
for i = 1:max_order+1
    subplot(max_order+1, 1, i);
    plot(real(kr), abs(res_exact(i,:)));
    hold on
    plot(real(kr), abs(res(i,:)), 'o');
end
