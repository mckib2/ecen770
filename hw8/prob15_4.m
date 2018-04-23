%% Problem 15.4
% Nicholas McKibben
% ECEn 770
% 2018-08-16

clear;
close all;

A = [ 1 1 0 1 0 0 0;
      0 1 1 0 1 0 0;
      0 0 1 1 0 1 0;
      0 0 0 1 1 0 1 ];

% Find H and then grab A2 from it
H = mod(rref(A),2);
A2 = H(:,5:7);

% Construct systematic G
G = [ A2; eye(3) ];
disp(G);

%% Verify we get the zero matrix
mod(H*G,2)
mod(A*G,2)