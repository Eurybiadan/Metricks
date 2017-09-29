function [ crossings ] = find_zero_crossings( signal )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

crossings = [];

for i=2:length(signal)
    
    if (sign(signal(i-1)) == 1 && ( (sign(signal(i)) == -1) ) ) || ...
       (sign(signal(i-1)) == -1 && ( (sign(signal(i)) == 1) ) )
    
        crossings = [crossings i];
    
    end

end

end

