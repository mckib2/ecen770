%% Homework 4
% Nicholas McKibben
% ECEn 770
% 2018-02-26

clear;
close all;

% Use any programming environment/language to program an
% encoder/decoder pair for the arithmetic code algorithm
% presented in class with symbols and associating probabilities
% listed in the order given. Let the source alphabet be:
chi    = [ 'A' 'B' 'V' 'S' 'E' 'R' '!' ];
p      = [ .15 .08 .02 .2  .29 .18 .08 ];
N = 16; % how many bits we talkin bout

% Encode the source output sequence:
%     ASBVRAERBSESEEERARSAEREESS
% using your program.
% Limit the length of codewords to 16 bits, and continue to
% add symbols to codewords until you can no longer add the next
% symbol and stay within the 16 bit codeword requirement.
% Remember, each codeword must end with ‘!’.

%% (a)
% Write the sequence of source outputs with ‘!’ symbols inserted where
% they must go to do the encoding as specified.

symbols = arithpre(N,'ASBVRAERBSESEEERARSAEREESS',[ 0 p ],chi);

%% (b)
% Provide a list of codewords in decimal and binary for the sequence
% above (fill in 0’s at the end to make 16-bit codewords if needed).

b = zeros(numel(symbols),1); l = b;
b_code = zeros(numel(symbols),N);
for ii = 1:numel(symbols)
    [ b(ii),l(ii) ] = arithenco1(symbols{ii},chi,[ 0 p ]);
    tmp = frac2bin(b(ii),ceil(-log2(l(ii)))+1);
    padnum = N - numel(tmp); % zero pad ending if need be
    b_code(ii,:) = [ tmp zeros(1,padnum) ];
    recieved{ii,:} = arithdeco1(b_code(ii,:),[ 0 p ],chi);
end

if ~isequal(recieved,symbols)
    fprintf('Something has gone awry, Batman!\n');
end

display(table(symbols,recieved)); % demonstrate that it works
codewords = cell(numel(symbols),1);
for ii = 1:size(b_code,1)
    codewords{ii} = strjoin(string(b_code(ii,:)),'');
end
display(table(symbols,codewords,b))

%% (c)
% Decode the following binary list of codewords:
% (001101100001010,001111011101011)

s1 = [ 0 0 1 1 0 1 1 0 0 0 0 1 0 1 0 0 ];
s2 = [ 0 0 1 1 1 1 0 1 1 1 0 1 0 1 1 0 ];

a1 = arithdeco1(s1,[ 0 p ],chi);
a2 = arithdeco1(s2,[ 0 p ],chi);
a = replace(strcat(a1,a2),chi(end),'');
fprintf('Decoded string: %s\n',a);

