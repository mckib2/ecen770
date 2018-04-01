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

% Try G14
% r_try = [ 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 ];
r_try = mod(G14(1,:) + G14(2,:),2);
c_try = decodeRM(r_try,G14);

%% (4) Performace
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
gammas = unique([ linspace(1,3,50) linspace(4,5,5) ]); % signal-to-noise ratio
N = [ 500*ones(1,45) 100*ones(1,10) ]; % num of errors
sigma2 = zeros(1,numel(gammas));
N0 = sigma2;
Pe = sigma2;

% Is this the prob?
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

    nn = 0;
    nbits = 0;
    while nn < N(ii)
        % Generate recieved codeword - assume all 0 codeword transmitted
        c = zeros(1,n13);
        
        % Generate noise
        p = crossprob(N0(ii));
        n = double(rand(1,n13) < p);
        r = c + n;
        
        % Add some bits
        nbits = nbits + k13;

        % Decode the recieved codeword
        c_hat = decodeRM(r,G13);
        
        % Now find the message - I think this is how you do it?
        m_hat = sign(c_hat/G13);
        m_hat(m_hat < 0) = 0;

        % Accumulate error
        nn = nn + sum(m_hat);
    end

    Pe(ii) = nn/nbits;

    % Update waitbar
    waitbar(ii/steps,h,sprintf('%2.2f%%',100*ii/steps));
end
sim_time = toc;
fprintf('Simulation took %f seconds to run.\n',sim_time);
delete(h); % remove wait bar

%% Plots
gammast = [ gammas linspace(5,15,15) ];
N0t = Eb./(r13*gammast);
Pe_theoretical = qfunc(sqrt(2*Eb./N0t));

x1 = 10*log10(Eb./N0);
x2 = 10*log10(Eb./N0t);
figure(1);
semilogy(x1,Pe,'k-','DisplayName','Simulated'); hold on; grid on;
semilogy(x2,Pe_theoretical,'k--','DisplayName','Theoretical');
title('Probability of Error');
xlabel('E_b/N_0 (dB)');
ylabel('P_b');
legend(gca,'show');
xlim([ 0 9 ]);

%% Coding Gain
fSNR_sim = fit(Pe',x1','poly2');
SNR_sim = fSNR_sim(10e-4);
fSNR_the = fit(Pe_theoretical',x2','poly2');
SNR_the = fSNR_the(10e-4);

coding_gain = abs(SNR_sim - SNR_the);
fprintf('The coding gain for Pb = 10^-4 is %f dB\n',coding_gain);