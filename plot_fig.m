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