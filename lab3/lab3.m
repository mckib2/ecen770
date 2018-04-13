%% Lab 3 - Viterbi
% Nicholas McKibben
% ECEn 770
% 2018-04-12

clear;
close all;

%% (1) Encoder
% Write a computer program that performs the encoding operation for a
% convolutional code with transfer function matrix
g1 = [ 1 0 1 ];
g2 = [ 1 1 1 ];

% Test the encoder
m = [ 1 1 0 0 1 0 1 ];
c = convencode(m,[ g1; g2 ]);
% disp(table(c.'));

%% (2) Hard decision Viterbi
% Program and test the hard-decision Viterbi decoder for this encoder. Your
% final result should be a BER performance plot showing the performance of
% your code as a function of Eb/N0. Include a curve that shows the
% theoretical performance for uncoded BPSK. You may want to review Lab 1
% material. Do not forget to incorporate the rate of the code when
% generating noise at specific SNR levels as you (hopefully) did in Lab 1.
% Make the y-axis log-scale

% Get input
r = [ 1 1 1 0 0 0 1 0 1 1 0 1 0 0 0 1 ];
r = [ 0 0 r ]; % put a couple zeros out front
l = 2; % length of sequence

% Initializations
si{1,:} = 1; % States - start in state 1
paths = zeros(1,4); % Paths, TODO: make into a class

% Set up state - next state table
s(1,:) = [ 1 2 ];
s(2,:) = [ 3 4 ];
s(3,:) = [ 1 2 ];
s(4,:) = [ 3 4 ];

% Set up state - prev state table
ps(1,:) = [ 1 3 ];
ps(2,:) = [ 1 3 ];
ps(3,:) = [ 2 4 ];
ps(4,:) = [ 2 4 ];

% Get the outputs
%          in=0   in=1
%    s     c1 c2  c1 c2
sout(1,:) = [ 0  0   1  1 ];
sout(2,:) = [ 0  1   1  0 ];
sout(3,:) = [ 1  1   0  0 ];
sout(4,:) = [ 1  0   0  1 ];

% For each recieved sequence
t = 1; % time is 1 indexed here to match MATLAB
for ii = (1+l):l:numel(r)
    t = t + 1;
    
    % Get revieved sequence rt, 2 bits
    rt = r(ii:(ii+l-1));
    
    % Evaluate all paths from s_{t-1} to st
    sprev = si{t-1};

    % Find the next states from sprev
    snext = s(sprev,:);
    snext = sort(snext(:)).';

    % Update next state
    si{t,:} = unique(snext);

    % Evaluate all paths
    for kk = 1:2 % for input 0 and 1
        for nn = 1:numel(sprev)
            % Get the Hamming distance for each path
            paths(t-1,sprev(nn),s(sprev(nn),kk),:) = pdist2(rt,sout(sprev(nn),(kk*l-1):(kk*l)),'Hamming')*l;

            val = paths{t-1,sprev(nn),s(sprev(nn),kk)};
            fprintf('%d -> %d, pdist = %d\n',sprev(nn),s(sprev(nn),kk),val{1});
            
            % If any snext has more than 1 incoming path, choose the path
            % with minimum path metric
            if size(paths,2) >= max(ps(s(sprev(nn),kk),:)) && size(paths,2) >= max(s(sprev(nn),kk))
                if numel(paths{t-1,ps(s(sprev(nn),kk),:),s(sprev(nn),kk)}) > 1
                    disp('We have 2 coming in!');
                    
                    % Choose the one that has a smaller path metric
                    [ val,idx ] = mink(paths{t-1,ps(s(sprev(nn),kk),:),s(sprev(nn),kk)},1);
                    
                    prev = ps(s(sprev(nn),kk),:);
                    next = s(sprev(nn),kk);
                    
                    % Remove path, denote this by making it -1
                    paths(t-1,prev(prev ~= prev(idx)),next,:) = -1;
                    
                    % Choose the correct path
                    paths(t-1,prev(idx),next,kk) = val;
                    
                    fprintf('Chose %d -> %d, pdist = %d\n',sprev(nn),s(sprev(nn),kk),paths(t-1,prev(idx),next,kk));
                end
            end
        end
    end

end
