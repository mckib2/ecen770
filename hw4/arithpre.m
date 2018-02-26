% There is a much better way to do this. I was trying
% to do something more elegant, but I ran out of time,
% so I did the easiest thing.  It works.  Just not
% very elegant.

function [ out ] = arithpre(N,symbols,p,chi,varargin)
    delim = ' ';
    out = [];
    ii = 1; start = 1;
    while ii <= numel(symbols)
        newsubstring = [ symbols(start:ii) chi(end) ];
        [ ~,l ] = arithenco1(newsubstring,chi,p);
        n = ceil(-log2(l)) + 1;
        
        if n == N
            %fprintf('Found next substring!\n');
            out = [ out newsubstring delim ];
            start = ii + 1;
        elseif n > N
            %fprintf('Found next substring!\n');
            out = [ out symbols(start:(ii-1)) chi(end) delim ];
            ii = ii - 1;
            start = ii + 1;
        elseif ii == numel(symbols)
            %fprintf('Found last substring!\n');
            out = [ out newsubstring ];
            break;
        end
        
        ii = ii + 1;
    end
    
    % Do some string-y stuff - spit out a cell matrix
    fprintf('With ''!'' symbols inserted: %s\n',out);
    out = split(out,delim);
end