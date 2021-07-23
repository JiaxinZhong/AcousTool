
lower_limit = 1;
upper_limit = 2e3+1i;

max_order = 30;
order = (0:max_order).';

int_num = (10:5e2);
int = 0 * order .* int_num;
tic
for i = 1:length(int_num)
    int(:,i) = IntHankelType1(max_order, lower_limit, upper_limit, 1, 'int_num', int_num(i));
end
toc

int(:,end)

fig = Figure;
plot(int_num, abs(int));




