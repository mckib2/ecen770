%% Lab 1 - Simulating a Communications Channel
% Nicholas McKibben
% ECEn 770
% 2018-03-06

clear;
close all;

%% (1) BPSK Simulation
% Write a program that will simulate a BPSK communication system with
% unequal prior bit probabilities.  Using your program, create data from
% which to plot the probability of bit error obtained from your simulation
% for SNRs in the range from 0 to 10 dB, for the three cases that P0 = .5,
% P0 = .25, and P0 = .1.  Decide on an approproate value of N.

% Following Algorithm 1.2...
% (1) Initialization
rng('default');
Eb = 1;
N = 5; % num of errors we wait for
P0 = [ .5 .25 .1 ];
gammas = linspace(1,10,15); % signal-to-noise ratio
sigma2 = zeros(numel(P0),numel(gammas));
N0 = sigma2;
Pe = sigma2;

% Let's see how long this is going to take...
h = waitbar(0,'Simulating...');
steps = numel(Pe);

% (2) for each signal-to-noise ratio gamma
for pp = 1:numel(P0)
    for ii = 1:numel(gammas)

        % (3) Compute N0, sigma2
        N0(pp,ii) = Eb/gammas(ii);
        sigma2(pp,ii) = N0(pp,ii)/2;

        % do
        nn = 0;
        nbits = 0;
        while nn < N
            % (5) Generate transmitted bits (for all P0s)
            nbits = nbits + 1;
            btx = double(rand(1) > P0(pp));

            % (6) Map the bits into the signal constellation
            s = btx;
            s(s == 0) = -sqrt(Eb);
            s(s == 1) = sqrt(Eb);

            % (7) Generate noise
            n = randn(1)*sqrt(sigma2(pp,ii));

            % (8) Add noise to get r
            r = s + n;

            % (9) Perform signal detection
            con = [ -sqrt(Eb) sqrt(Eb) ].';
            [ ~,idx ] = min((con - r).^2);
            s_hat = con(idx);

            % (10) Determine brx
            brx = s_hat;
            brx(brx > 0) = 1;
            brx(brx < 1) = 0;

            % (11) Compare tx to rx
            if brx ~= btx
                % (12) Accumulate error
                nn = nn + 1;
            end
        end

        Pe(pp,ii) = nn/nbits;
        
        % Update waitbar
        currstep = sub2ind(size(Pe.'),ii,pp);
        waitbar(currstep/steps,h,sprintf('%2.2f%%',100*currstep/steps));
    end
end
close(h);

%% (2) Theoretical
% Prepare data from which to plot the theoretical probability of error
% (1.24) for the same three values of P0.
tau = 0;
Pe_theoretical = zeros(size(Pe));
for pp = 1:numel(P0)
    Pe_theoretical(pp,:) = Q((tau + sqrt(Eb))./sqrt(sigma2(pp,:)))*P0(pp) ...
        + Q((sqrt(Eb) - tau)./sqrt(sigma2(pp,:))).*(1 - P0(pp));
end

%% (3) Plots
% Plot the simulated probability of error on the same axes as the
% theoretical probability of error.  The plots should have Eb/N0 in dB as
% the horizontal axis and the probability as the vertical axis, plotted on
% a logarithmic scale.

figure(1);
for pp = 1:numel(P0)
    semilogy(10*log10(Eb./N0(pp,:)),Pe_theoretical(pp,:), ...
        'DisplayName',sprintf('P0_%d',pp));
    hold on;
end
title('Probability of error for BPSK signalling');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of error');
legend('show');

%% Helper Functions
function [ out ] = Q(x)
    out = zeros(1,numel(x));
    for ii = 1:numel(x)
       out(ii) = 1/sqrt(2*pi)*integral(@probfun,x(ii),Inf);
    end
end

function [ out ] = probfun(n)
    out = exp(-n.^2/2);
end