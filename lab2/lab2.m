%% Lab 2: Reed-Muller Codes
% Nicholas McKibben
% ECEn 770
% 2018-03-31

clear;
close all;

%% (1) Generator Matrices
% Find the generator matrices for three different Reed-Muller codes: (1)
% RM(1, 3), (2) RM(1, 4), (3) RM(2, 4).

G13 = [ 1 1 1 1 1 1 1 1;
        0 0 0 0 1 1 1 1;
        0 0 1 1 0 0 1 1;
        0 1 0 1 0 1 0 1 ];
    
G14 = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
        0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
        0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
        0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
        0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 ];

G24 = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
        0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
        0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
        0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
        0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1;
        0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1;
        0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1;
        0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1;
        0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1;
        0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 ];
    

%% (2) Interesting Quantities
% Give the dimension, blocklength, rate, and minimum distance for each of
% the three codes from Part 1.

% For (1,3)
k13 = dim(1,3);
dmin13 = dmin(1,3);
n13 = bl(1,3);
r13 = rate(1,3);

% For (1,4)
k14 = dim(1,4);
dmin14 = dmin(1,4);
n14 = bl(1,4);
r14 = rate(1,4);

% For (2,4)
k24 = dim(2,4);
dmin24 = dmin(2,4);
n24 = bl(2,4);
r24 = rate(2,4);

%% (3) Encode/Decode
% Program and test encoders and decoders for only the first two codes.

% Encode - the easy part
enc = @(m,G) m*G;

r = [ 1 0 0 1 0 0 1 0 ];
dec(r,G13)

%% Helper Functions
function [ k ] = dim(r,m)
    k = 1;
    for ii = 1:r
        k = k + nchoosek(m,r);
    end
end

function [ val ] = dmin(r,m)
    val = 2^(m-r);
end

function [ n ] = bl(~,m)
    n = 2^m;
end

function [ r ] = rate(r,m)
    r = dim(r,m)/bl(r,m);
end

% Get the bipolar representation
function [ R ] = bipolar(r)
    R = (-1*ones(size(r))).^r;
end

% Let's follow Algorithm 8.1, because we have (1,m) codes
function [ c ] = dec(r,G)
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