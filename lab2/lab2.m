%% Lab 2: Reed-Muller Codes
% Nicholas McKibben
% ECEn 770
% 2018-03-31

if exist('h','var')
    delete(h);
end
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

% Test out Book Example 8.9
r_example = [ 1 0 0 1 0 0 1 0 ];
c_example = decodeRM(r_example,G13);

disp(r_example);
disp(c_example);

% Try G14
% r_try = [ 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 ];
r_try = mod(G14(1,:) + G14(2,:),2);
c_try = decodeRM(r_try,G14);

%% (4) Performace of RM(1,3)
% Display the performance curves as functions of Eb/N0 for the first two
% codes on the same plot, along with a theoretical curve for uncoded BPSK
% (you may want to review Coding Project 1). Identify the coding gain as
% measured at probability of bit error 10e−4 . Simulate proper ranges of
% Eb/N0 to give good estimates for the coding gains.

% Following Algorithm 1.2...
% (1) Initialization
rng('default');
Ec = 1;
Eb = (n13/k13)*Ec;
% gammas = unique([ linspace(1,3,50) linspace(4,5,5) ]); % signal-to-noise ratio
% gammas = linspace(0,8,9);
gammas = 8.4;
N = 150*ones(size(gammas));
sigma2 = zeros(1,numel(gammas));
N0 = sigma2;
Pe13 = sigma2;

% Probability of flipping
crossprob = @(N00) qfunc(sqrt(2*Ec/N00));

% Let's see how long this is going to take...
h = waitbar(0,'Simulating...');
steps = numel(gammas);

tic;
% For each signal-to-noise ratio gamma
for ii = 1:numel(gammas)

    % Compute N0, sigma2
    N0(ii) = Ec/(r13*gammas(ii));
    sigma2(ii) = N0(ii)/2;
    p = crossprob(N0(ii));
    
    nn = 0;
    nbits = 0;
    while nn < N(ii)
        % Generate recieved codeword - assume random codeword to deal with
        % max situation in MATLAB...
        m = randi(2,1,k13) - 1;
        c = mod(m*G13,2);
        
        % Generate noise
        n = double(rand(1,n13) < p);
        r = mod(c + n,2);
        
        % Add some bits
        nbits = nbits + k13;

        % Decode the recieved codeword
        %c_hat = decodeRM(r,G13_systematic);
        c_hat = mod(decodeRM(r,G13),2);
        
        % Now find the message
        %m_hat = c_hat(1:k13);
        m_hat = mod([ c_hat(1) (c_hat(5) + c_hat(1)) (c_hat(3) + c_hat(1)) (c_hat(2) + c_hat(1)) ],2);
        %disp(m_hat);
        
        
        % Accumulate error
        nn = nn + sum(mod(m_hat + m,2));
        
        
        %waitbar((ii-1)/steps,h,sprintf('nn = %d and %2.2f%%',nn,100*(ii-1)/steps));
    end

    Pe13(ii) = nn/nbits;

    % Update waitbar
    waitbar(ii/steps,h,sprintf('%2.2f%%',100*ii/steps));
end
sim_time = toc;
fprintf('Simulation took %f seconds to run.\n',sim_time);
delete(h); % remove wait bar


%% Plots
gammast = [ gammas linspace(max(gammas),max(gammas)*2,15) ];
N0t = Eb./(gammast);
Pe_theoretical = qfunc(sqrt(2*Eb./N0t));

% x13 = 10*log10(Eb./(r13*N0));
x13 = 10*log10(Eb./N0);
x2 = 10*log10(Eb./N0t);
figure(1);
semilogy(10*log10(gammas),Pe13,'k-','DisplayName','Simulated RM(1,3)'); hold on; grid on;
semilogy(x2,Pe_theoretical,'k--','DisplayName','Theoretical');
title('Probability of Error');
xlabel('E_b/N_0 (dB)');
ylabel('P_b');
legend(gca,'show');

%% (4) Performace of RM(1,4)
% Display the performance curves as functions of Eb/N0 for the first two
% codes on the same plot, along with a theoretical curve for uncoded BPSK
% (you may want to review Coding Project 1). Identify the coding gain as
% measured at probability of bit error 10e−4 . Simulate proper ranges of
% Eb/N0 to give good estimates for the coding gains.

% Following Algorithm 1.2...
% (1) Initialization
rng('default');
Ec = 1;
Eb = (n14/k14)*Ec;
sigma2 = zeros(1,numel(gammas));
N0 = sigma2;
Pe14 = sigma2;

% Let's see how long this is going to take...
h = waitbar(0,'Simulating...');
steps = numel(gammas);

tic;
% For each signal-to-noise ratio gamma
for ii = 1:numel(gammas)

    % Compute N0, sigma2
    N0(ii) = Ec/(r14*gammas(ii));
    sigma2(ii) = N0(ii)/2;
    p = crossprob(N0(ii));
    
    nn = 0;
    nbits = 0;
    while nn < N(ii)
        % Generate recieved codeword - same deal as before
        m = randi(2,1,k14) - 1;
        c = mod(m*G14,2);
        
        % Generate noise
        n = double(rand(1,n14) < p);
        r = mod(c + n,2);
        
        % Add some bits
        nbits = nbits + k14;

        % Decode the recieved codeword
        %c_hat = decodeRM(r,G14_systematic);
        c_hat = mod(decodeRM(r,G14),2);
        
        % Now find the message
        %m_hat = c_hat(1:k14);
        m_hat = mod([ c_hat(1) (c_hat(9) + c_hat(1)) (c_hat(5) + c_hat(1)) (c_hat(3) + c_hat(1)) (c_hat(2) + c_hat(1)) ],2);

        % Accumulate error
        nn = nn + sum(mod(m_hat + m,2));
    end

    Pe14(ii) = nn/nbits;

    % Update waitbar
    waitbar(ii/steps,h,sprintf('%2.2f%%',100*ii/steps));
end
sim_time = toc;
fprintf('Simulation took %f seconds to run.\n',sim_time);
delete(h); % remove wait bar

%% Plots
x14 = 10*log10(Eb./(N0));
% x14 = 10*log10(Eb./(r14*N0));
semilogy(x14,Pe14,'k-.','DisplayName','Simulated RM(1,4)'); hold on; grid on;
legend(gca,'show');
