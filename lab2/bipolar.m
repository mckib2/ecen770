% Get the bipolar representation
function [ R ] = bipolar(r)
    R = (-1*ones(size(r))).^r;
end