%% Lab 3 - Viterbi
% Nicholas McKibben
% ECEn 770
% 2018-04-19

clear;
close all;

%% (1) Encoder
% Write a computer program that performs the encoding operation for a
% convolutional code with transfer function matrix
g1 = [ 1 0 1 ];
g2 = [ 1 1 1 ];

% Test the encoder
N = 10000;
m = randi([ 0 1 ],[ 1 N ]);
c = convencode(m,[ g1; g2 ]);

rng('default');

% Received signal here
r = c;

% From binary to decimal:
R = bin2dec(join(string(reshape(r,2,[]).'),''));

% Set the initial path costs, assume we start in state 0
nu = 4;
M(1) = 0;
M(2:8) = inf;

% Set P = zero set for initial paths
P = cell(8,1);
for ii = 1:8
    P{ii} =  [];
end

% set t = 0
t = 1;

% Define a transition,output matrices
% prev  this   next  output
% 0   ->  0  ->  0     00
% 2              1     11
%
% 0   ->  1   -> 2     01
% 2              3     10
%
% 1   ->  2   -> 0     11
% 3              1     00
%
% 1   ->  3   -> 2     10
% 3              3     01
s = zeros(nu,2);
s(1,:) = [ 1 2 ]; % 0 -> 0,1
s(2,:) = [ 3 4 ]; % 1 -> 2,3
s(3,:) = [ 1 2 ]; % 2 -> 0,1
s(4,:) = [ 3 4 ]; % 3 -> 2,3

% state -> prev states
sp = zeros(nu,2);
sp(1,:) = [ 1 3 ]; % 0 -> 0,2
sp(2,:) = [ 1 3 ]; % 1 -> 0,2
sp(3,:) = [ 2 4 ]; % 2 -> 1,3
sp(4,:) = [ 2 4 ]; % 3 -> 1,3

% Does state transition exist
se = zeros(nu);
se(1,[ 1 2 ]) = 1; % 0 -> 0,1
se(2,[ 3 4 ]) = 1; % 1 -> 2,3
se(3,[ 1 2 ]) = 1; % 2 -> 0,1
se(4,[ 3 4 ]) = 1; % 3 -> 2,3

% Output
o = ones(nu)*NaN;
o(1,1) = 0; o(1,2) = 3;
o(2,3) = 1; o(2,4) = 2;
o(3,1) = 3; o(3,2) = 0;
o(4,3) = 2; o(4,4) = 1;

% path decoder
pd(1,1) = 1;
pd(1,2) = 2;
pd(2,3) = 3;
pd(2,4) = 4;
pd(3,1) = 5;
pd(3,2) = 6;
pd(4,3) = 7;
pd(4,4) = 8;

first = true;
while (t <= numel(R))

    %fprintf('-----------------t = %d ----------------\n',t);
    
    % Construct the end of each existing path
    pii = [];
    for ii = 1:numel(P)
        if ~isempty(P{ii})
            p = P{ii};
            pii = [ pii p(end) ];
        else
            pii = [ pii NaN ];
        end
    end
    
    % For each next state q at time t + 1
    for ii = 1:numel(P)
        if ~isnan(pii(ii)) || first
            ps = pii(ii);
            if first
                % Start in state 0
                ps = 1;
            end
            idx = ii;
            prepend = []; inheritance = 0;
            for q = s(ps,:)
                if isinf(M(idx))
                    M(idx) = 0;
                end
                val = M(idx);
                M(idx) = inheritance + M(idx) + sum(abs(dec2bin(R(t),2) - dec2bin(o(ps,q),2)));
                
                %fprintf('%d -> %d, R: %d, out: %d, metric: %d\n',ps,q,R(t),o(ps,q),M(idx));
                P{idx} = [ P{idx} prepend q ];
                
                % After we continue the path, make a new path
                idx = find(isnan(pii),ii);
                idx = idx(end);
                if first
                    idx = idx + 1;
                    first = false;
                end
                prepend = P{ii};
                prepend = prepend(1:end-1);
                inheritance = val;
            end
        end
    end
    
    t = t + 1;
    
    if t >= 4
        % All paths have been allocated - meaning there are now two paths
        % going into every node.  We need to select the ones with minimum
        % metric.
        
        % Find all the ending nodes
        pii = [];
        for ii = 1:numel(P)
            p = P{ii};
            pii = [ pii p(end) ];
        end
        
        % Now find the node groupings
        qs(1,:) = find(pii == 1);
        qs(2,:) = find(pii == 2);
        qs(3,:) = find(pii == 3);
        qs(4,:) = find(pii == 4);
        
        % Now remove the paths with the best metric
        Pnew = cell(size(P));
        Mnew = ones(size(M))*inf;
        for ii = 1:nu
            [ ~,goodNode ] = min(M(qs(ii,:)));
            goodNode = qs(ii,goodNode);
            Pnew{ii} = P{goodNode};
            Mnew(ii) = M(goodNode);
        end
        P = Pnew;
        M = Mnew;
    end
end

% Choose the path with least metric at the end
[ ~,idx ] = min(M);
vitPath = P{idx};
vitPath = [ 1 vitPath ]; % assumed state 0 at beginning
for ii = 2:numel(vitPath)
    strbin = dec2bin(o(vitPath(ii-1),vitPath(ii)),2);
    r_hat(ii-1,:) = strbin - '0';
end

% Now find the message
m_hat = circshift(mod(r_hat(:,1) + r_hat(:,2),2),-1);

% Accumulate the error
e = sum(mod(m_hat(1:end-1).' + m(1:end-1),2));


function [ c ] = convencode(m,g)

    % Convolve for each g(j)
    for ii = 1:size(g,1)
        c(ii,:) = mod(conv(m,g(ii,:)),2);
    end
    
    % Interleave the rows
    c = reshape(c,1,[]);
    c = c(1:(numel(m)*size(g,1)));
end
