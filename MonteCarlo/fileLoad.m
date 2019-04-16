function [data] = fileLoad(filename)
%{

    Authors: Jashan Chopra (107689146)
    Date Created: April 4th, 2019

Script Purposes and goals: This simple function will load a file of Static
test stand data, specifically for usage in the ASEN2004 Bottle Rocket Lab

%}

data = load(filename);      % first we load the fileLoad
data = data(:,3);           % we want just the summed load column
data = data * 4.44822;      % convert summed force to newtons
