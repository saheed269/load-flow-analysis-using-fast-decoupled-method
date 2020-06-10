clear all;
close all;
clc;
ch = input('Enter the bus system no.: (6 or 9 or 14 or 26 or 30): ');
while ch ~= 6 && ch ~= 14 && ch ~= 26 && ch ~= 30 && ch ~= 9
    fprintf('Invalid Input try again\n');
    ch = input('Enter the bus system no.: (6 or 9 or 14 or 26 or 30): ');
end
switch ch
    case 6
        data6
    case 9
        data9
    case 14
        data14
    case 26
        data26
    case 30
        data30
end

basemva = 100; 
accuracy = 0.001;  %Maximum Error
maxiter = 50;   %Maximum iteration
ybus  %form the bus admittance matrix
tic
decouple1 %Load flow solution by fast decoupled method(Algorithm)
toc
busout   %Prints the power flow solution on the screen
lineflow %Computes and displays the line flow and losses

