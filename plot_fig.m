%% plot figures 

result_data=cell(3,2);
gurobi_ip=cell(3,1);
%prediction hoizon 
result_data{1,1}=[8 10 30 45;0.1954 0.2042 0.6179 0.8864;0.0419 0.0590 0.1318 0.2052]';
result_data{1,2}=[8 10 30 45;1.1501 1.572 3.6162 5.2541;0.1887 0.221 0.4452 0.7612]';
gurobi_ip{1,1}=[8 10;1.9436 2.7031;2.2394 3.1234]';

figure(1)
subplot(2,2,[1 3])
plot(result_data{1,1}(:,1),result_data{1,1}(:,2),'-s');
hold all;
plot(result_data{1,1}(:,1),result_data{1,1}(:,3),'-s');
plot(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,2),'-s');
legend('GPAD-5.e-3','GPAD-5.e-2','Gurobi-IP')
xlabel('Prediction horizion');
ylabel('Average time (sec)');
subplot(2,2,[2 4])
plot(result_data{1,2}(:,1),result_data{1,2}(:,2),'-s');
hold all;
plot(result_data{1,2}(:,1),result_data{1,2}(:,3),'-s');
plot(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,3),'-s');

legend('GPAD (\epsilon=5\cdot 10^{-3}','GPAD(\epsilon=5\cdot 10^{-3}','Gurobi-IP')
xlabel('Prediction horizion');
ylabel('Maximum time (sec)');
%scenarios
result_data{2,1}=[60 120 240 480;0.1408 0.2042 0.2578 0.4619;0.0331 0.0590 0.0792 0.1424]';
result_data{2,2}=[60 120 240 480;0.6304 1.572 1.840 1.905;0.0975 0.221 0.2312 0.2860]';
gurobi_ip{2,1}=[60 120;1.1278 2.7031;1.2917 3.1234]';
figure(2)
subplot(2,2,[1 3])
plot(result_data{2,1}(:,1),result_data{2,1}(:,2),'-s');
hold all;
plot(result_data{2,1}(:,1),result_data{2,1}(:,3),'-s');
plot(gurobi_ip{2,1}(:,1),gurobi_ip{2,1}(:,2),'-s');
legend('GPAD (\epsilon_g=5\cdot 10^{-3})','GPAD(\epsilon=5\cdot 10^{-3}','Gurobi-IP')
xlabel('Scenarios');
ylabel('Average time (sec)');
subplot(2,2,[2 4])
plot(result_data{2,2}(:,1),result_data{2,2}(:,2),'-s');
hold all;
plot(result_data{2,2}(:,1),result_data{2,2}(:,3),'-s');
plot(gurobi_ip{2,1}(:,1),gurobi_ip{2,1}(:,3),'-s');

legend('GPAD-5.e-3','GPAD-5.e-2','Gurobi-IP')
xlabel('Scenarios');
ylabel('Maximum time(sec)');
%system_dimensions 
result_data{3,1}=[5 10 30 40;0.0420 0.2042 3.2187 6.2842;0.0154 0.0590 0.6685 1.3245]';
result_data{3,2}=[5 10 30 40;0.2232 1.572 10.2429 22.5273;0.0602 0.221 2.5011 4.2920]';
gurobi_ip{3,1}=[5 10;0.8149 2.7031;1.2019 3.1234]';


figure(3)
subplot(2,2,[1 3])
semilogy(result_data{3,1}(:,1),result_data{3,1}(:,2),'-s');
hold all;
semilogy(result_data{3,1}(:,1),result_data{3,1}(:,3),'-s');
semilogy(gurobi_ip{3,1}(:,1),gurobi_ip{3,1}(:,2),'-s');
legend('GPAD-5.e-3','GPAD-5.e-2','Gurobi-IP')
xlabel('Masses');
ylabel('Average time (sec)');
subplot(2,2,[2 4])
semilogy(result_data{3,2}(:,1),result_data{3,2}(:,2),'-s');
hold all;
semilogy(result_data{3,2}(:,1),result_data{3,2}(:,3),'-s');
semilogy(gurobi_ip{3,1}(:,1),gurobi_ip{3,1}(:,3),'-s');

legend('GPAD-5.e-3','GPAD-5.e-2','Gurobi-IP')
xlabel('Masses');
ylabel('Maximum time(sec)');

%% Not Normalised data 
result_data=cell(3,2);
gurobi_ip=cell(3,1);
%prediction hoizon 
result_data{1,1}=[128 256 512 1024 2048 4096;...
    0.0806 0.1098 0.2595 0.2449 0.7626 0.9803;...
    0.0499 0.0693 0.1268 0.1492 0.2429 0.5551;...
    0.0256 0.0343 0.0471 0.0753 0.1218 0.2284]';
result_data{1,2}=[128 256 512 1024 2048 4096;
    0.6098 1.0187 1.2957 2.5806 5.5187 7.0746;...
    0.3124 0.7586 0.5479 1.1991 1.0136 2.8822;...
    0.1004 0.2408 0.2625 0.3503 0.3732 0.6753]';
gurobi_ip{1,1}=[128 256 512 1024 2048 4096;...
    0.4661 0.9886 1.6445 3.3629 6.7626 14.015;...
    0.7297 1.8533 2.8263 6.004 9.7667 21.4]';

speedup=[51.4 51.8 68.86 73.07 83.146 78.68];
figure(1)
subplot(2,2,[1 3])
semilogx(result_data{1,1}(:,1),result_data{1,1}(:,2),'-s');
hold all;
semilogx(result_data{1,1}(:,1),result_data{1,1}(:,3),'-s');
semilogx(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s');
semilogx(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,2),'-s');
legend('GPAD-5.e-3','GPAD-1.e-2','GPAD-5.e-2','Gurobi-IP')
xlabel('Binary tree scenarions');
ylabel('Average time (sec)');
subplot(2,2,[2 4])
semilogx(result_data{1,2}(:,1),result_data{1,2}(:,2),'-s');
hold all;
semilogx(result_data{1,2}(:,1),result_data{1,2}(:,3),'-s');
semilogx(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s');
plot(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,3),'-s');

legend('GPAD-5.e-3','GPAD-1.e-2','GPAD-5.e-2','Gurobi-IP')
xlabel('Binary tree scenarions');
ylabel('Maximum time (sec)');

figure(2)
plot(result_data{1,1}(:,1),speedup);
xlabel('Binary tree scenarios');
ylabel('Speedup')

%% Normalised data 
result_data=cell(3,2);
gurobi_ip=cell(3,1);
%prediction horizon
result_data{1,1}=[10 20 30 40 50 60;...
    0.2489 0.4140 0.6913 0.8435 1.3154 1.5154 ;...
    0.1568 0.2616 0.4446 0.5568 0.8210 0.9710 ;...
    0.0764 0.1611 0.2769 0.3571 0.4601 0.5730 ]';
result_data{1,2}=[10 20 30 40 50 60;...
    1.9333 2.723 3.867 5.2063 7.5283 5.2695;...
    1.2681 1.054 2.565 2.7959 4.0952 3.8943;...
    0.1745 0.503 1.011 0.8983 1.3693 1.2693]';
gurobi_ip{1,1}=[10 20 30 40 50 60;...
    5.5781 10.8177 15.1898 20.6732 22.3625 26.6361;...
    7.4076 13.3254 20.9680 24.6445 28.5038 33.8569]';

figure(1)
semilogy(result_data{1,1}(:,1),result_data{1,1}(:,2),'-s','Linewidth',2,'MarkerSize',10);
hold all;
semilogy(result_data{1,1}(:,1),result_data{1,1}(:,3),'-s','Linewidth',2,'MarkerSize',10);
semilogy(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s','Linewidth',2,'MarkerSize',10);
semilogy(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,2),'-s','Linewidth',2,'MarkerSize',10);
legend('APG-0.005','APG-0.01','APG-0.05','Gurobi-IP')
xlabel('prediction horizon','FontSize',20);
ylabel('average time (sec)','FontSize',20);
axis tight;
grid on;

figure(2)
semilogy(result_data{1,2}(:,1),result_data{1,2}(:,2),'-s','Linewidth',2,'MarkerSize',10);
hold all;
semilogy(result_data{1,2}(:,1),result_data{1,2}(:,3),'-s','Linewidth',2,'MarkerSize',10);
semilogy(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s','Linewidth',2,'MarkerSize',10);
semilogy(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,3),'-s');
legend('APG-0.005','APG-0.01','APG-0.05','Gurobi-IP')
xlabel('prediction horizon','FontSize',20);
ylabel('average time (sec)','FontSize',20);
axis tight;
grid on;

%% Normalised data 
result_data=cell(3,2);
gurobi_ip=cell(3,1);
%Scenarios   
result_data{1,1}=[128 256 512 1024 2048 4096 8192;...
    0.0433 0.0805 0.0795 0.1350 0.2419 0.6822 1.3135;...
    0.0316 0.0512 0.0573 0.1020 0.1707 0.4483 0.8109;...
    0.0225 0.0249 0.0325 0.0621 0.1067 0.2453 0.3970]';
result_data{1,2}=[128 256 512 1024 2048 4096 8192;
    0.3148 0.4627 0.6497 0.7482 1.0050 3.8080 5.9283;...
    0.1543 0.4251 0.4234 0.4850 0.8289 2.850  3.7551;...
    0.1131 0.1408 0.0993 0.1623 0.2934 0.7938 1.7386]';
gurobi_ip{1,1}=[128 256 512 1024 2048 4096 8192;...
    0.3636 0.9884 1.4969 3.1814 6.9752 13.6275 32.8611;...
    0.4317 1.5497 2.0914 4.1012 10.035 20.8478 46.7915]';

speedup=[42 53.6 62.86 68.55 81.23 78.53 83.16];
figure(1)
subplot(2,2,[1 3])
loglog(result_data{1,1}(:,1),result_data{1,1}(:,2),'-s','Linewidth',2,'MarkerSize',10);
hold all;
loglog(result_data{1,1}(:,1),result_data{1,1}(:,3),'-s','Linewidth',2,'MarkerSize',10);
loglog(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s','Linewidth',2,'MarkerSize',10);
loglog(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,2),'-s','Linewidth',2,'MarkerSize',10);
legend('GPAD-5.e-3','GPAD-1.e-2','GPAD-5.e-2','Gurobi-IP')
xlabel('Binary tree scenarions','FontSize',20);
ylabel('Average time (sec)','FontSize',20);
axis tight;
subplot(2,2,[2 4])
loglog(result_data{1,2}(:,1),result_data{1,2}(:,2),'-s','Linewidth',2,'MarkerSize',10);
hold all;
loglog(result_data{1,2}(:,1),result_data{1,2}(:,3),'-s','Linewidth',2,'MarkerSize',10);
loglog(result_data{1,1}(:,1),result_data{1,1}(:,4),'-s','Linewidth',2,'MarkerSize',10);
loglog(gurobi_ip{1,1}(:,1),gurobi_ip{1,1}(:,3),'-s');

legend('GPAD-5.e-3','GPAD-1.e-2','GPAD-5.e-2','Gurobi-IP')
xlabel('Binary tree scenarions','FontSize',20);
ylabel('Maximum time (sec)','FontSize',20);
axis tight;
figure(2)
plot(result_data{1,1}(:,1),speedup,'Linewidth',2,'MarkerSize',10);
xlabel('Binary tree scenarios','FontSize',20);
ylabel('Speedup')
axis tight;