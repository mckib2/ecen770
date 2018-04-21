%% Lab 3 - Viterbi
% Nicholas McKibben
% ECEn 770
% 2018-04-19

if exist('h','var')
    delete(h);
end
clear;
close all;

load('hard.mat');
Pehard = Pe;

mbits = 1e4;
% [ r,r_hat,m,m_hat,e ] = viterbialg(N);

% Simulation:
Ec = 1;
Eb = 2*Ec;
gammas = linspace(1,2,5);
N = 2500;
sigma2 = zeros(1,numel(gammas));
N0 = sigma2;
Pe = sigma2;

% Probability of flipping
rng('default');
crossprob = @(N00) Q(sqrt(2*Ec/N00));

% Let's see how long this is going to take...
h = waitbar(0,'Simulating...');
steps = numel(gammas);

tic;
% For each signal-to-noise ratio gamma
for ii = 1:numel(gammas)

    % Compute N0, sigma2
    N0(ii) = Ec/(.5*gammas(ii));
    sigma2(ii) = N0(ii)/2;
    p = crossprob(N0(ii));
    
    nn = 0;
    nbits = 0;
    while nn < N
        % Generate recieved codeword - assume random codeword
        [ ~,~,~,~,e ] = viterbialg(mbits,p,N0(ii)/2,'soft');
        
        % Add some bits
        nbits = nbits + mbits - 1;

        % Accumulate error
        nn = nn + e;
        
        if nn
            fprintf('nn: %d\n',nn);
        end
    end

    Pe(ii) = nn/nbits;

    % Update waitbar
    waitbar(ii/steps,h,sprintf('%2.2f%%',100*ii/steps));
end
sim_time = toc;
fprintf('Simulation took %f seconds to run.\n',sim_time);
delete(h); % remove wait bar


%% Plots
gammast = 0:10;
N0t = Eb./(gammast);
Pe_theoretical = Q(sqrt(2*Eb./N0t));

x2 = 10*log10(Eb./N0t);
figure(1);
semilogy(10*log10(linspace(1,3.5,5)),Pehard,'k--','DisplayName','Hard Decision');
hold on; grid on;
semilogy(10*log10(gammas),Pe,'k-.','DisplayName','Soft Decision');
semilogy(x2,Pe_theoretical,'k-','DisplayName','Theoretical');
title('Probability of Error');
xlabel('E_b/N_0 (dB)');
ylabel('P_b');
legend(gca,'show');

