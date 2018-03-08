%% Helper Functions
function [ out ] = Q(x)
    out = zeros(1,numel(x));
    for ii = 1:numel(x)
       out(ii) = 1/sqrt(2*pi)*integral(@probfun,x(ii),Inf);
    end
end

function [ out ] = probfun(n)
    out = exp(-n.^2/2);
end