%% Matrix solving - computational cost
clear
clc
close 'all'

exponents = 9:12;
time_req = nan(1,length(exponents));

for k = 1:length(exponents)
   N = 2^exponents(k);
   A = rand(N,N); b = rand(N,1);
   f = @() A\b;
   time_req(k) = timeit(f);
   
end

figure(1)
loglog(2.^exponents,time_req,'-b',...
    2.^exponents,1e-4*(2.^exponents).^3,'--r',...
    'linewidth',3)
title('Computational Time vs. N','fontsize',16,...
    'fontweight','bold');
xlabel('N','fontsize',14,'fontweight','bold');
ylabel('Time [s]','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on
