function [ num ] = bin2frac(code)
    num = 0;
    for ii = 1:numel(code)
        num = num + code(ii)*.5^(ii);
    end
end