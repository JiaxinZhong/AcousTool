% function data = ParsePulseData(fn)

%     fn_num = size(fn);
    

clear all
fn = {'pls1.txt', 'pls2.txt'};

for i = 1:length(fn)
    data{i} = GetPulseAsciiFile(fn{i});
end

%% plot 
figure
plot(data{1}.freq, data{1}.spl)
hold on
plot(data{2}.freq, data{2}.spl);
legend('File 1', 'File 2');
xlabel('Frequency (Hz)');
ylabel('SPL (dB)');