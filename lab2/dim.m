function [ k ] = dim(r,m)
    k = 1;
    for ii = 1:r
        k = k + nchoosek(m,r);
    end
end