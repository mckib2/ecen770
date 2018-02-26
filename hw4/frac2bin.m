function [ b ] = frac2bin(num,N)
    b = zeros(1,N);
    for ii = 1:N
        if (num - 2^(-ii)) >= 0
            b(ii) = 1;
            num = num - 2^(-ii);
        end
    end
end