%% Helper Functions
% function [ out ] = Q(x)
%     out = zeros(1,numel(x));
%     for ii = 1:numel(x)
%        out(ii) = 1/sqrt(2*pi)*integral(@probfun,x(ii),Inf);
%     end
% end
% 
% function [ out ] = probfun(n)
%     out = exp(-n.^2/2);
% end
function [ p ] = Q(xlist)
    p = zeros(size(xlist));
    i = 0;
    for x = xlist
        i = i+1;
        if(x < 0)
            p(i) = 1- 0.5*erfc(-x/sqrt(2));
        else
            p(i) = 0.5*erfc(x/sqrt(2));
        end
    end
end