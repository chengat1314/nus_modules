function homework3_2()
%% Facial Image Based Glasses Wearing Recognition 
% Apply the Multi-layer perceptron to the glass-wearing recognition problem 
% using batch mode of learning.
clc,clear;close all
n=20;times = 50;
train_error = zeros(times,n);
test_error = zeros(times,n);
for i = 1:50
    error = Multi_layer_perceptron(n);
    train_error(i,:) = error(1,:);
    test_error(i,:) = error(2,:);
end

plot(mean(train_error),'*-');hold on
plot(mean(test_error),'o-');legend('train error','test error');
xlabel('n');ylabel('mse');title('training error plot');

