                   % ====== Experimental ====== % 
                   % === For Uplink System ==== %
clear all; clc ; warning off ;
N_frame = 15;      % Number of frames/packet
N_packet= 1200;   % Number of packets
IT = N_packet;
Eb=1;
MS=N_frame;
NR     = 3; % N  % Number of B.S.'s antenna for each user 
N_user = 2; % K  % Number of user
NT     = 3; % M  % Number of user's antenna

N=NT;    K=N_user;   M=NR;     
b=2;                        % Number of bits per QPSK symbol
N_pbits = N_frame*K*b;      % Number of bits per packet
N_tbits = N_pbits*N_packet; % Number of total bits

SNRdBs = [-5:2:25]; 
K_dB = 30 ; % for Rician channel = [ 0 ... 30 ]
    B= 0.85  ; % beta -> correction parameter =[ 0.1 ... 0.9 ]
  

    
                     hand=waitbar(0,'Please Wait....');  
                     
 parfor SNR=1:length(SNRdBs);
                
    SNRdB=SNRdBs(SNR);
    scf=sqrt(0.5*10^(-SNRdB*0.1));sq=1/sqrt(2);
   N_ebits=0; rand('seed',1); randn('seed',1);
   
    
   
   
   for count=1:N_packet; 
                 waitbar((((SNR-1)*N_packet)+count-1)/(15*N_packet));
        %% =============== Base Station (B.S.) ==============
                 %------------symbols-----------
        msg_bit=randint(N_pbits,1);  % (N_frame*K*b X 1) 40x1
        symbol = QPSK_mapper(msg_bit);
        s=reshape(symbol,K,N_frame);
      %------------ channels between B.S.& Users ------------
      
        H1=(randn(NR,NT)+1i*randn(NR,NT)); %/sqrt(2); %Rayleigh channel between B.S. & user 1) 2x2
         %  Rician_ch1 = Ric2deig_model(K_dB,M,N); H1 = Rician_ch1; 
          
        H2=(randn(NR,NT)+1i*randn(NR,NT)); %/sqrt(2); % Rayleigh channel between B.S. & user 2
          % Rician_ch2 = Ric2deig_model(K_dB,M,N); H2 = Rician_ch2;
          
        % hsr=(randn(1,1)+1i*randn(1,1)); %/sqrt(2); % Rayleigh channel betwe S & R 
          Rician_chsr = Ric2deig_model(K_dB,1,1);hsr = Rician_chsr;
          
            %------- extended matrices ( interference channels ) --------
        Hi1= (1/sqrt(2)).*(randn(NR,NT)+1i*randn(NR,NT)); %for user 1  2x2
        Hi2= (1/sqrt(2)).*(randn(NR,NT)+1i*randn(NR,NT)); %for user 2
       
        %--------- pre-coding (w) -----------
        
        [V1,D1] = eig(H1'*H1,Hi1'*Hi1);% v1=V1;d=D1;
         w1=V1(:,N);%/sqrt(V1(:,N)'*V1(:,N));% 2*1
        
        [V2, D2] = eig(H2'*H2,Hi2'*Hi2); %v2=V2;d2=D2;
         w2=V2(:,N);%/sqrt(V2(:,N)'*V2(:,N));% 2x1
        
        
        [V12,D12] = eig(Hi1'*Hi1,H1'*H1);% v1=V1;d=D1;
         w12=V12(:,N);%/sqrt(V1(:,N)'*V1(:,N));% 2*1
         
        
         [V21,D21] = eig(Hi2'*Hi2,H2'*H2);% v1=V1;d=D1;
          w21=V21(:,N);%/sqrt(V1(:,N)'*V1(:,N));% 2*1
        
        
        %--------- transmit vector x ---------------
       
        x1=w1*s(1,:);%2*1
        x2=w2*s(2,:);
          %------------ noise vector ------------
        
        n1=scf*(randn(M,1)+1i*randn(M,1));% 2x10
        noise1=repmat(n1,1,N_frame);
        
        n2=scf*(randn(M,1)+1i*randn(M,1));
        noise2=repmat(n2,1,N_frame);
        
        noise12=scf*(randn(NR,1)+1i*randn(NR,1));NoiS1= repmat(noise12(1,1),[1 N_frame 1]);% [heq2 heq2 heq2 heq2 ...]
                                                 
        noise22=scf*(randn(NR,1)+1i*randn(NR,1));NoiS2= repmat(noise22(1,1),[1 N_frame 1]);
       %------------ ------------
        
      %  msr=real(hsr);
       % Ehsr=sum(msr.^2)/NT;% si = NT
       % Nosr=Eb*Ehsr/SNRdB ;
       % BT=sqrt(1/(Ehsr *Eb + Nosr));
         BT=0.85;
        %------------ Receive vector y in U1 from U2 
        s1co = hsr*s(2,:)+NoiS1; % in user 1
        yu1  = ((hsr')/(hsr'*hsr))* s1co ;
        %------------ Receive vector y in U2 from U1
        s2co = hsr*s(1,:)+NoiS2; % in user 2
        yu2  = ((hsr')/(hsr'*hsr))* s2co ;
                % ------------- GAIN --------------
        GAIN_2D= (w21'*Hi2')/(w21'*Hi2'*Hi2*w21);
        GAIN_1D= (w12'*Hi1')/(w12'*Hi1'*Hi1*w12);
        
        %------------ Receive vector y in B.S from user 1 & user 2 
        y1=  H1*x1+Hi2*x2+noise1; % from user 1
        y11= Hi2*w21*(BT*yu2)+ noise1;
        
        
         y2=  H2*x2+Hi1*x1+noise2; % from user 2
         y22= Hi1*w12*(BT*yu1)+ noise2;
      
        GAIN_1= (w1'*H1')/(w1'*H1'*H1*w1); % Equalizer for the 1st user (ML)   
        GAIN_2= (w2'*H2')/(w2'*H2'*H2*w2); % Equalizer for the 2nd user
      
        %--------------- SIC for user1 signal ------------
       S1_hat= GAIN_1 * y1 ;
       y1_hat= y1- H1*w1*S1_hat;
       S22_hat= ((w2'*Hi2')/(w2'*Hi2'*Hi2*w2))*y1_hat;
        SS1_hat = GAIN_2D* y11 ;
        
       S2_hat=GAIN_2* y2;
       y2_hat= y2-H2* w2* S2_hat;
       S11_hat= ((w1'*Hi1')/(w1'*Hi1'*Hi1*w1))* y2_hat;
       SS2_hat = GAIN_1D* y22;
       
                               %--------
          T_S1_hat = S1_hat + (SS1_hat + B* S11_hat); % NO inter-user channel with beta
          T_S2_hat = S2_hat + (SS2_hat + B* S22_hat);
        
        %---------------------
        S_T = [ T_S1_hat ; T_S2_hat ];
        
        symbol_hat = reshape(S_T,K*N_frame,1);
        
        symbol_sliced = QPSK_slicer(symbol_hat);
        demapped = QPSK_demapper(symbol_sliced);
       p = msg_bit;%(1:2,:); 
       q = demapped;
        N_ebits = N_ebits + sum(p'~=q);

 end
    BER(SNR) = N_ebits/N_tbits;
 
end
 close(hand);
figure
semilogy(SNRdBs,BER,'r-o'); grid on
hold on
axis([-5 25 10^-5 10^0]);
%colours = ['k-o';'k-*';'k-s';'k-+';'k-^';'k-h';'k-v';'k-p'];
%colours = ['b-o';'r-d';'g-s';'k-v';'m-^';'b-<';'r->';'g-p'];
%title(['B.S.a(M)=',num2str(M),', n.user(K)=',num2str(K),', user.a(N)=',num2str(N),', Rician(k)=',num2str(K_dB),',B=',num2str(B)]);
title([',B.S.a(M)=',num2str(M),', n.user(K)=',num2str(K),', user.a(N)=',num2str(N),', Rayleigh',', B=',num2str(B) ]);
legend(['Up.SLR.T=',num2str(N),'-R=',num2str(M),' IT=',num2str(N_packet),' MS=',num2str(N_frame)],['BF-TR=',num2str(N),'-RC=',num2str(M)]);
% the end