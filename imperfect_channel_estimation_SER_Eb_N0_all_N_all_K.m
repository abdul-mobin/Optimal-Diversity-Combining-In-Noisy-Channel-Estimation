 % Objective: Investigating the performance of optimum diversity combining in
 %           presence of imperfect channel estimation. 
 % 
 %           by Mohammad Abdul Mobin
 %                original:4.6.2023
 %                edited:  7.27.2024
 % dependencies: N/A

%for all N & particular M, plot for value of K=[0,5,10]

clc
clear all
close all

%inputs
N=[1,4];                  %Diversity (Number of Receiver Antennas)
M=2;                      %M-ary PSK
error=0;                  %Error Count
total_frames=1e5;         %Each frame containing 101 symbols

for Idx_N=1:2

%Initialization
source_data= zeros(101,1);
symbol     = zeros(101,1);
u          = zeros(N(Idx_N),1);
h          = zeros(N(Idx_N),1);
y          = zeros(101,N(Idx_N));
yn         = zeros(101,N(Idx_N));
noise      = zeros(101,N(Idx_N));
h_estimated= zeros(1,N(Idx_N));
alpha      = zeros(1,101);
detected_symbol= zeros(101,1);
SER        =zeros(3,11);
SER_analytical  =zeros(3,11);

Eb_N0=[0:2:20];             
K_db=[0,5,10];               %Rician Factor in dB
    
for Idx_K=1:3

k=10^(K_db(Idx_K)/10);   %Rician Factor in pure number

for Idx_Eb_N0=1:11
    error=0;
    
    for q=1:total_frames                      %Sending frames each with 101 symbols, 1st symbol of each frame is pilot
        %Source
        source_data=randi([1,M],101,1);
        symbol = exp(j*(source_data-1)*((2*pi)/M));
        
        %channel
        u= sqrt(k/2)*(1+j)*ones(N(Idx_N),1);  % Mean of the channel, if Rayleigh channel is used, mean is zero
        h= u + (1/sqrt(2))*(randn(N(Idx_N),1) +j* randn(N(Idx_N),1));
        
        y=symbol*h.';                         % Signal Trasmitted, each symbol is being multiplied by N channel gains
        
        %noise
        N0=1/((10^(Eb_N0(Idx_Eb_N0)/10))*log2(M));
        noise=sqrt(N0/2)*(randn(101,N(Idx_N)) + j*randn(101,N(Idx_N)));
        
        yn=y+noise;                           % Signal Received by the Receiver
        
        %channel estimation
        pilot=symbol(1,1);                    % First symbol is being used as pilot
        h_estimated= yn(1,:)/pilot;
        
        %detection 
        ro_square=1/(1+N0);                   %covariance co-efficient
        alpha = (ro_square*h_estimated.' + (1-ro_square)*u)'*yn.';      % Decision variable, this variable is derived by Maximum a posteriori(MAP) rule. Ref: eqn 16 
        
                                               % calculating eucledian distance between all possible symbols and decision variabls
                                        
        for i=2:101
            min=1e5;
            for w=1:M
                a= ( abs( alpha(i) - exp(j*(w-1)*((2*pi)/M))))^2;
                if a<min
                    min=a;
                    symb=w;
                end
            end
            detected_symbol(i)=symb;
        end
        
        %error calculation
        for i=2:101
            if( source_data(i) ~= detected_symbol(i) )
                error=error+1;
            end
        end
    end
    
    % Analytical error calculation
    gamma=(k+1)/N0;
    gamma_rice1=(ro_square*gamma)/(gamma*(1-ro_square)+1+k);  %gamma of rician channel derived in paper 
    
    fun1 = @(x) exp(((N(Idx_N)*k)/ro_square)./(1+((gamma_rice1*(sin(pi/M))^2)./(sin(x)).^2)))./(1+((gamma_rice1*(sin(pi/M))^2)./(sin(x)).^2)).^N(Idx_N);   
    
    SER(Idx_K,Idx_Eb_N0)=error/(total_frames*100);
    SER_analytical(Idx_K,Idx_Eb_N0) = (exp(-(N(Idx_N)*k)/ro_square)/pi)*integral(fun1,0,(pi-(pi/M)));
 
end
end


%plotting
semilogy(Eb_N0,SER(1,:),'^')
xlabel('EB/N0(db)');
ylabel('log10(SER)');
grid on
hold on
semilogy(Eb_N0,SER_analytical(1,:),'-')
hold on
semilogy(Eb_N0,SER(2,:),'*')
hold on
semilogy(Eb_N0,SER_analytical(2,:),'.-')
hold on
semilogy(Eb_N0,SER(3,:),'o')
hold on
semilogy(Eb_N0,SER_analytical(3,:),'--')
hold on

end
legend('k=0, simulation', 'k=0, analytical', 'k=5db, simulation', 'k=5db, analytical', 'k=10db, simulatio', 'k=10db, analytical');