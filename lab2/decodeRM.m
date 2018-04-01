% Let's follow Algorithm 8.1, because we have (1,m) codes
function [ c ] = decodeRM(r,G)
    % (4) Find the bipolar representation of r
    R = bipolar(r);
    
    % (5) Compute the Hadamard transform
    T = R*hadamard(size(G,2));
    
    % (6) Find the coordinate ti with the largest magnitude
    [ ~,idx ] = max(abs(T));
    
    % (7) Let i have the binary expansion (i_m,..i_1)
    % Make sure to 0 index and convert from string!
    ii = dec2bin(idx - 1,log2(numel(r))) - '0';
    
    % (8) Build c, assume 1 is not sent
    c = zeros(1,size(G,2));
    for jj = 1:length(ii)
        % Recognize that rows of G are v's
        c = c + ii(numel(ii) - jj + 1)*G(end - jj + 1,:);
    end
    c = mod(c,2); % binary addition
    
    % % (10) 1 is sent - complement all the bits
    if T(idx) < 0
        c = 1 - c;
    end
end