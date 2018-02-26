function [ a ] = arithdeco1(code,p,chi,bnds)
    % get code as a fraction
    in = bin2frac(code);

    % define boundaries
    if nargin == 3
        bnds = [ 0 1 ];
    end
    % find our cummulative  sum
    cp = cumsum_p(p,bnds);
    
    % find the bin
    [ ~,idx ] = min((cp - in).^2);
    if in < cp(idx) % find the left index
        idx = idx - 1;
    end
    
    % check ending condtion
    if chi(idx) == '!'
        % fprintf('I found the end!\n');
        a = chi(idx);
    else
        % forge ahead, give new bounds
        bnds = [ cp(idx) cp(idx+1) ];
        a0 = chi(idx);
        a = [ a0 arithdeco1(code,p,chi,bnds) ];
    end
end

