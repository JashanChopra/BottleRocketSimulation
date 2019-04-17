%% Housekeeping
clc; 
close all;
clear;
%% Reading in data
AP = readtable('yeet5 - Sheet1.csv');
AP = sortrows(AP,23); %sort the data by airspeed
AP = table2array(AP); %convert table to numbers
AP(~any(AP,2),:) = [];
for i = 1:length(AP)
if (AP(i,4) > 7) && (AP(i,4)<9)
    AP1(i,:) = AP(i,:);
elseif (AP(i,4) > 14) && (AP(i,4)<16)
    AP2(i,:) = AP(i,:);
elseif (AP(i,4) > 19) && (AP(i,4)<21)
    AP3(i,:) = AP(i,:);
end
end
%get rid of empty spaces in the matrices
AP1(~any(AP1,2),:) = [];
AP2(~any(AP2,2),:) = [];
AP3(~any(AP3,2),:) = [];

[C_d1,attack_angle1] = drag(AP1);
[C_d2,attack_angle2] = drag(AP2);
[C_d3,attack_angle3] = drag(AP3);
C_d = (C_d1+C_d2+C_d3)/3;
% figure
% plot(attack_angle1,C_d1,'Linewidth',2) %plot C_d for each speed
% hold on
% plot(attack_angle2,C_d2,'Linewidth',2)
% plot(attack_angle3,C_d3,'Linewidth',2)
% xlabel('Angle of Attack, radians')
% ylabel('Coeffiecient of Drag')
% title('C_d vs Angle of Attack')
% legend('14 m/s','20 m/s','25 m/s','Location','northwest')
