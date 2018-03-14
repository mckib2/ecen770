%% Lab 1 - 8-PSK Portion
% Nicholas McKibben
% ECEn 770
% 2018-03-07

clear;
close all;

%% (1) 8-PSK Simulation
% Write a program that will simulate an 8-PSK communication system with
% equal prior bit probabilities. Use a signal constellation in which the
% points are numbered in Gray code order. Make your program so that you can
% estimate both the symbol error probability and the bit error probability.
% Decide on an appropriate value of N.

% Follow algorithm 1.2...
% (1) Initialization
rng('default');
M = 8;
Eb = 1;
Es = Eb*log2(M);
A = sqrt(Es);
N = 500; % num of errors we wait for
gammas = linspace(1,50,30); % signal-to-noise ratio
sigma2 = zeros(1,numel(gammas));
N0 = sigma2;
Ps = sigma2;

% Let's see how long this is going to take...
h = waitbar(0,'Simulating...');
steps = numel(Ps);

tic;
% (2) for each signal-to-noise ratio
for ii = 1:numel(gammas)
    
    % (3) Compute N0, sigma2
    N0(ii) = Es/gammas(ii);
    sigma2(ii) = N0(ii)/2;
    
    % (4) do
    nn = 0;
    nbits = 0;
    while nn < N
        % (5) Generate transmitted bits
        nbits = nbits + log2(M);
        btx = bin2dec(strjoin(string(randi(2,1,log2(M)) - 1)));
        
        
        % (6) Map the bits into the signal constellation, Gray encoding
        s = btx;
        s(s == 7) = 1;
        s(s == 6) = A*exp(1j*pi/4);
        s(s == 2) = A*exp(1j*pi/2);
        s(s == 3) = A*exp(1j*3*pi/4);
        s(s == 1) = A*exp(1j*pi);
        s(s == 0) = A*exp(1j*5*pi/4);
        s(s == 4) = A*exp(1j*3*pi/2);
        s(s == 5) = A*exp(1j*7*pi/4);
        idxs = [ 7 6 2 3 1 0 4 5 ];
        
        
        % (7) Generate noise
        n = randn(1)*sqrt(sigma2(ii)) + 1j*randn(1)*sqrt(sigma2(ii));
        
        % (8) Add noise to get r
        r = s + n;
        
        % (9) Perform signal detection
        con = A*exp(1j*[ 0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 ]).';
        [ ~,idx ] = min((con - r).^2);
        s_hat = con(idx);
        
        % (10) Determine brx
        brx = idxs(idx);
        
        % (11) Compare tx to rx
        if brx ~= btx
            % (12) Accumulate error
            brx = dec2bin(brx,log2(M));
            btx = dec2bin(btx,log2(M));
            bb = brx - btx;
            
            nn = nn + numel(bb(bb == 0));
        end
    end
    
    Ps(ii) = nn/nbits;
    
    % Update waitbar
    waitbar(ii/steps,h,sprintf('%2.2f%%',100*ii/steps));
end
sim_time = toc;
fprintf('Simulation took %f seconds to run.\n',sim_time);
delete(h); % remove wait bar

%% (2) Theoretical
% Prepare data from which to plot the bound on the probability of symbol
% error Ps using (1.26) and probability of bit error Pb using (1.27).

dmin = norm(A*[ 1 0 ] - A*[ 1/sqrt(2) 1/sqrt(2) ]);
Ps_theoretical_book = 2*Q(dmin./(2*sqrt(sigma2)));

% Ps_theoretical = erfc(sqrt(log2(M)*Eb./N0)*sin(pi/M));
Ps_theoretical = 2*Q(sqrt(2*Es./N0)*sin(pi/M));

figure(1);
subplot(2,1,1);
x = 10*log10(Eb./N0);
semilogy(x,Ps_theoretical/log2(M)); grid on; hold on;
semilogy(x,Ps_theoretical_book/log2(M));
xlabel('E_b/N_0');
ylabel('P_b');

subplot(2,1,2);
x = 10*log10(Es./N0);
semilogy(x,Ps_theoretical); grid on; hold on;
semilogy(x,Ps_theoretical_book);
xlabel('E_s/N_0');
ylabel('P_s');

%% (3) Plots
% Plot the simulated probability of symbol error and bit error on the same
% axes as the bounds on the probabilities of error.

figure(2);
semilogy(x,Ps,'DisplayName','Simulated P_s'); grid on; hold on;
semilogy(x,Ps_theoretical_book,'DisplayName','Theoretical P_s');
legend(gca,'show');

%% (4) Compare
% Compare the theoretical and simulated results. Comment on the accuracy of
% the bound compared to the simulation and the amount of time it took to
% run the simulation.