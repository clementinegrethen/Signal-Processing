%% Projet traitement du signal (grethen/Fernandez)

load donnees1.mat;
load donnees2.mat;
%définition des paramètres utiles:
fp1=0;%porteuse 2
fp2=46e3;%porteurse1
Fe=120e3;%fréquence échantillonage
Te=1/Fe;%période echantiollonage
T=40e-3;% durée d'un timeslot
N=length(bits_utilisateur1)% longueur message binaire
Ts = T/N;
Ns = Ts/Te;
% 3.1 implantation
m1=kron(2*bits_utilisateur1-1,ones(1,Ns));
m2=kron(2*bits_utilisateur2-1,ones(1,Ns));
t=(0:N*Ns-1)/Fe
figure(1)
plot(t,m1);
xlabel('secondes');
ylabel('Amplitude');
figure(2);
plot(t,m2);
xlabel('seconde');
ylabel('Amplitude');
%3.1: densité spectrale:
f=(0:N*Ns-1)*Fe
fft_m1=fftshift(fft(m1));
dsp_m1=(abs(fft_m1).^2)*1/N
figure(3)
plot(f,dsp_m1);
xlabel("Hz");
ylabel("Densité spectrale de m1");
fft_m2=fftshift(fft(m2));
dsp_m2=(abs(fft_m2).^2)*1/N
figure(4)
plot(f,dsp_m2)
xlabel("Hz");
ylabel("Densité spectrale de m2");

%3.2.2 construction du signal MF_TDMA
%(a)

signal_1=zeros(1,N*Ns*5);
signal_2=zeros(1,N*Ns*5);
signal_1(N*Ns+1:2*N*Ns) = m1;
signal_2(4*N*Ns+1:5*N*Ns) = m2;
figure(5);
plot([0:Te:Te*(5*N*Ns-1)],signal_1);
xlabel("temps en secondes");
ylabel("signal envoyé sur la porteuse pour m1");
figure(6)
plot([0:Te:Te*(5*N*Ns-1)],signal_2)
xlabel("temps en secondes");
ylabel("signal envoyé sur la porteuse  m2");
%(b)
x1=signal_1.*cos(2*pi*fp1*[0:Te:Te*(5*N*Ns-1)]);
x2=signal_2.*cos(2*pi*fp2*[0:Te:Te*(5*N*Ns-1)]);
%2
somme_signal=x1+x2;
SNR = 100;
Ps = mean(abs(somme_signal).^2);
Pb = Ps * 10^(-(SNR/10));
Bruit = (sqrt(Pb))*randn(1,5*N*Ns);
x = somme_signal + Bruit;
figure(7);
plot([0:Te:Te*(5*N*Ns-1)],x);
xlabel("fréquence en Hz");
ylabel("signal MF-TDMA");

%3:densité spectrale
dsp_x=(1/length(x))*abs(fft(x)).^2;
figure(8)
plot([0:Te:Te*(5*N*Ns-1)],fftshift(dsp_x));
fft_x=fft(x);
xlabel("fréquence en Hz");
ylabel("densité spectrale du signal MF-TDMA");

%4:Mise en place du récepteur MF-TDMA:
%4.1.1 synthèse du filtre passe-bas:

ordre=61;
fc=23e3;
I_61=[-Te/2*(ordre-1):Te:Te/2*(ordre-1)];
sinus_card_61=2*fc/Fe*sinc(I_61*2*fc);

signal_filtre=conv(x,sinus_card_61,'same');
figure(50)
plot([0:Te:Te*(5*N*Ns-1)],signal_filtre);
xlabel('Temps en secondes');
ylabel("signal filtré par le passe bas")

%tracé rep impulsionnelle du filtre passe-bas:
figure(9)
plot(I_61,sinus_card_61);
xlabel('Temps en secondes');
ylabel("réponse impulsionnelle du PB")
% réponse en fréquence du  filtre passe-bas en échelle logarithmique:

n = 2^nextpow2(length(x));
H= fftshift(fft(sinus_card_61,n));
figure(10)
semilogy(linspace(-Fe/2,Fe/2,length(H)),abs(H));
xlabel('fréquence en Hz');
ylabel("réponse en fréquence du PB")
%tracé de la densité spectrale de puissance 
figure(11)
semilogy(linspace(-Fe/2,Fe/2,length(fft_x)),fftshift(dsp_x));
hold on,semilogy(linspace(-Fe/2,Fe/2,length(H)),abs(H));
title('Figure 11 :  desniét spectrale de x et réponse du filtre');
legend('DSP_x','réponse du filtre');
%synthèse du filtre passe-haut:
H_passe_haut=-sinus_card_61;
H_passe_haut(1+(length(H_passe_haut)-1)/2)= 1 + H_passe_haut(1+(length(H_passe_haut)-1)/2);
signal_passe_haut = conv(x,H_passe_haut,'same');

figure(13)
plot([0:Te:Te*(5*N*Ns-1)],signal_passe_haut);
xlabel('Temps en secondes');
ylabel("signal filtré par le passe haut")
H_PH=fftshift(fft(H_passe_haut,n));
figure(14)
plot(I_61,H_passe_haut);
xlabel('Temps en secondes');
ylabel("Réponse impulsionnelle du passe haut")
figure(15)
semilogy(linspace(-Fe/2,Fe/2,length(H_PH)),abs(H_PH));
xlabel('Hz');
ylabel('réponse en fréquence du filtre passe haut')
figure(16);
figure,semilogy(linspace(-Fe/2,Fe/2,length(fft_x)),fftshift(dsp_x));
hold on,semilogy(linspace(-Fe/2,Fe/2,length(H_PH)),abs(H_PH));
legend('DSP_x','module de la réponse du filtre');

%4.1.3: filtrage:
%recopier les signaux au dessus
x2_filtre=signal_passe_haut;
x1_filtre=signal_filtre;

%4.2 retour en  bande de phase
x1_base=x1_filtre.*cos(2*pi*fp1*[0:Te:Te*(5*N*Ns-1)])
x2_base=x2_filtre.*cos(2*pi*fp2*[0:Te:Te*(5*N*Ns-1)])
%4.3: détection du slot utile:
p1 = reshape(x1_base,[length(x1_base)/5,5]);
p2 = reshape(x2_base,[length(x1_base)/5,5]);

%détecteur d'énergie:
E_1=sum(abs(p1).^2)
E_2 =sum((abs(p2)).^2)

%recherche du maximum:
max_1=find(E_1==max(E_1))
max_2=find(E_2==max(E_2))

%détection slot utile
MessageRetrouve1=p1(:,max_1);
MessageRetrouve2=p2(:,max_2);

%demodulation bande de base:
SignalFiltre1=filter(ones(1,Ns),1,MessageRetrouve1) ;
SignalEchantillonne1=SignalFiltre1(Ns :Ns :end) ;
BitsRecuperes1=(sign(SignalEchantillonne1)+1)/2 ;
texte1=bin2str(BitsRecuperes1);
disp(texte1);

SignalFiltre2=filter(ones(1,Ns),1,MessageRetrouve2) ;
SignalEchantillonne2=SignalFiltre2(Ns :Ns :end) ;
BitsRecuperes2=(sign(SignalEchantillonne2)+1)/2 ;
texte2=bin2str(BitsRecuperes2);
disp(texte2);










