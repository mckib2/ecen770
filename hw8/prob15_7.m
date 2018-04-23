%% Probelm 15.7
% Nicholas McKibben
% ECEn 770
% 2018-04-16

clear;
close all;

A = [ 1 1 1 0 0 1 1 0 0 1;
      1 0 1 0 1 1 0 1 1 0;
      0 0 1 1 1 0 1 0 1 1;
      0 1 0 1 1 1 0 1 0 1;
      1 1 0 1 0 0 1 1 1 0 ];

  %% (a)
r = [ 1 0 1 0 1 0 0 1 0 1 ];

c_hat = r;

z = mod(c_hat*A.',2);
disp(table(z));

v = z*A;
disp(v);

%% (b)
r = [ 1 0 1 0 1 1 0 1 0 1 ];

c_hat = r;
z = mod(c_hat*A.',2);
disp(table(z));

v = z*A;
fprintf('v\n\t');
disp(v);

[ ~,idx ] = maxk(v,3);
c_hat(idx) = ~c_hat(idx);

z = mod(c_hat*A.',2);
disp(table(z));
v = z*A;
disp(v);

[ ~,idx ] = maxk(v,1);
c_hat(idx) = ~c_hat(idx);

z = mod(c_hat*A.',2);
disp(table(z));
