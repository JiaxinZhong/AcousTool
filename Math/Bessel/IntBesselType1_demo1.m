
lower_limit = 0;
upper_limit = 100+1i;

max_order = 3;
order = (0:max_order).';

int_num = (10:5e2);
int = 0 * order .* int_num;
tic
for i = 1:length(int_num)
    int(:,i) = IntBesselType1(max_order, lower_limit, upper_limit, 1, 'int_num', int_num(i));
end
toc


fig = Figure;
plot(int_num, abs(int));




