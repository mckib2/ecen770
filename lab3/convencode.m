function [ c ] = convencode(m,g)

    % Convolve for each g(j)
    for ii = 1:size(g,1)
        c(ii,:) = mod(conv(m,g(ii,:)),2);
    end
    
    % Interleave the rows
    c = reshape(c,1,[]);
end