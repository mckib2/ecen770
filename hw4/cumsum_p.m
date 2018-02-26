function [ cp ] = cumsum_p(p,bnds)
    cp = zeros(size(p));

    span = bnds(2) - bnds(1);
    cp(1) = bnds(1); cp(end) = bnds(2);

    for ii = 2:(numel(cp)-1)
        cp(ii) = bnds(1) + span*sum(p(1:ii));
    end
end