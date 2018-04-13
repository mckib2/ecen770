classdef ppath < handle
    properties
        s1 % from state
        s2 % to state
        m  % path metrix
    end
    
    methods
        function [ obj ] = ppath(s1,s2)
            obj.s1 = s1;
            obj.s2 = s2;
        end
    end
    
end