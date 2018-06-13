%Setting up the parameters

S0 = 999; %Initial number of survivors
K = 1000; %Total population
Beta_true = 0.005; 
StandDev = 10;

% Generating Data

n = 50; %Number of data points generated
Y = zeros(1,n); %Vector of the population at time t
T = zeros(1,n); %Vector of times t

for i = 1:n
    t = i*(10/n);
    Y(i) = min(999, max(1, S(t, Beta_true) + random('norm',0,StandDev)));
    T(i) = t;
end


%Metropolis Hastings

BETAS = zeros(1, 1); %This is the vector of accepted Betas in the MCMC
BetaI = random('uniform', 0, Beta_true*10); %BetaI is the parameter added to the BETAS vector at every iteration

for mcmc_iteration = 1:100000

    BETAS(mcmc_iteration) = BetaI;
    
    Betac = BetaI; %Betac is the current Beta (last value of beta that has been accepted)
    Betan = random('norm',Betac,0.0001); %Betan is the new beta being tested
    
    %Alpha is the ratio of the posterior of Betan by the posterior of Betac
    %By the uniform Prior distribution of Beta, if the random Beta isn't in
    %[0,1], then we set Alpha to be zero, making this new Beta
    %automatically rejected.
    
    if Betan > 1  
        Alpha=0;
    else
        if Betan < 0 
            Alpha = 0; 
        else
            Alpha = exp(LogLikelihood(Betan, Y, T, StandDev, n) - LogLikelihood(Betac, Y, T, StandDev, n));
        end
    end
        
            
    u = random('uniform', 0, 1); %Generating a random number uniformly chosen between 0 and 1
    if u > Alpha
        BetaI = BetaI; %Rejecting the tested Beta if the random number u is bigger than Alpha
    else
        BetaI = Betan; %Accepting it otherwise
    end
    
    
end


histogram(BETAS,10000)
xlim([Beta_true*0.9,Beta_true*1.1])


% % Defining symbolic functions for plots
% 
% syms t
% S_plot001(t) = (S0*K) / (S0 + (K - S0)*exp(0.001*K*t));
% Z_plot001(t) = K - S_plot001(t);
% 
% S_plot002(t) = (S0*K) / (S0 + (K - S0)*exp(0.002*K*t));
% Z_plot002(t) = K - S_plot002(t);
% 
% % Plotting figures 1&2
% 
% figure
% 
% subplot(2,1,1)
% hold on
% fplot(S_plot001,[0 10], 'k')
% fplot(Z_plot001,[0 10], '--k')
% hold off
% xlabel('Days')
% ylabel('Population')
% legend('S(t)','Z(t)')
% title('Beta = 0.001')
% 
% subplot(2,1,2)
% hold on
% fplot(S_plot002,[0 10], 'k')
% fplot(Z_plot002,[0 10], '--k')
% hold off
% xlabel('Days')
% ylabel('Population')
% legend('S(t)','Z(t)')
% title('Beta = 0.002')


% % Plotting figure 3

% hold on
% fplot(S_plot001,[0 10])
% scatter(T,Y)
% xlabel('Days')
% ylabel('Population')
% hold off

%Defining the loglikelihood

function l = LogLikelihood(beta, Y, T, StandDev, n)
    l=0;
    for i = 1:n
        l = l -((Y(i)-S(T(i), beta))^2)/(2*StandDev^2)/sqrt(2*pi*StandDev^2);
    end
end


%Defining the solution of the ODE

function s = S(t, beta)
    S0 = 999;
    K = 1000;
    s = (S0*K) / (S0 + (K - S0)*exp(beta*K*t));
end

