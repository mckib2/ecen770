function [ code,l ] = arithenco1(symbols,chi,p,varargin)
    % which symbol are we encoding?
    if length(varargin) < 1
        ii = 1;
    else
        ii = varargin{1};
    end
    % define boundaries
    if length(varargin) < 2
        bnds = [ 0 1 ];
    else
        bnds = varargin{2};
    end
    % find our cummulative sum
    cp = cumsum_p(p,bnds);
    
    % Find a number for the current symbol
    s = symbols(ii);  % curr symbol

    % find the bin
    idx = find(chi == s);

    % check recursion conditions
    if s == chi(end)
        % we've hit the bang operator
        code = (cp(idx) + cp(idx+1))/2;

        % length of the boundary
        l = cp(idx+1) - cp(idx);
    else
        % find the number in the middle and take it
        bnds = [ cp(idx) cp(idx+1) ];
        [ code,l ] = arithenco1(symbols,chi,p,ii+1,bnds);
    end
end