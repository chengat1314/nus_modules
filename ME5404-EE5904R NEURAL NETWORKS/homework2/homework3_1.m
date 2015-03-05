function homework3_1()
%% Facial Image Based Glasses Wearing Recognition 
% Apply the Rosenblatt’s perceptron (single layer perceptron) 
clc,clear;close all
times =50;
record_error = zeros(times,2); 
for i = 1:times
    error =  single_perceptron();
    record_error(i,1) = error(1);
    record_error(i,2) = error(2);
end
plot(record_error(:,1),'*-');hold on
plot(record_error(:,2),'o-');legend('train error','test error');
xlabel('n');ylabel('mse');title('training error plot');
record_error
mean(record_error)
