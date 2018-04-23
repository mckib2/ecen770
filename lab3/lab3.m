%% Lab 3 - Viterbi
% Nicholas McKibben
% ECEn 770
% 2018-04-19

if exist('h','var')
    delete(h);
end
clear;
close all;

method = 'soft';

if strcmp(method,'hard')
    load('soft.mat');
    Pesoft = Pe;
    gammas_soft = gammas;
else
    load('hard.mat');
    Pehard = Pehard;
    gammas_hard = gammas;
end
    
mbits = 1e4;
% [ r,r_hat,m,m_hat,e ] = viterbialg(N);

% Simulation:
Ec = 1;
Eb = 2*Ec;
gammas = linspace(1,5,4);
N = [ 3000 100 5 2 ];
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
    while nn < N(ii)
        % Generate recieved codeword - assume random codeword
        [ ~,~,~,~,e ] = viterbialg(mbits,p,sqrt(sigma2(ii)),method);
        
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
if strcmp(method,'hard')
    Pehard = Pe;
    gammas_hard = gammas;
else
    Pesoft = Pe;
    gammas_soft = gammas;
end

gammast = 0:10;
N0t = Eb./(gammast);
Pe_theoretical = Q(sqrt(2*Eb./N0t));

x2 = 10*log10(Eb./N0t);
figure(1);
semilogy(10*log10(gammas_hard),Pehard,'k--','DisplayName','Hard Decision');
hold on; grid on;
semilogy(10*log10(gammas_soft),Pesoft,'k-.','DisplayName','Soft Decision');
semilogy(x2,Pe_theoretical,'k-','DisplayName','Theoretical');
title('Probability of Error');
xlabel('E_b/N_0 (dB)');
ylabel('P_b');
legend(gca,'show');

