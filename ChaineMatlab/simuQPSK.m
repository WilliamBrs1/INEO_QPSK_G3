
 % clear all
 % close all

%Pour Octave : installer éventuellement le package octave-signal
%isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
%if isOctave
%    pkg load signal
%end

constelQPSK = [1; 1i; -1; -1i] * exp(1i * pi / 4);

% génération filtre demi-Nyquist
rolloffNyquist = 0.2;
nbSymbHalfNyq = 40;
overSamplingFactor = 8;
filterHalfNyquist = CalculHalfNyquist(overSamplingFactor, rolloffNyquist, nbSymbHalfNyq);
Lfilter = length(filterHalfNyquist);

figure(1);
freqz(filterHalfNyquist, 1, 8192);

nbSymb = 10000;
SNRmin = -5;
SNRmax = 10;
%SNRmin = 20;
%SNRmax = 21;
SNRStep = 0.5;

TEB = [];


for snrdB = SNRmin:SNRStep:SNRmax

    %% génération des symboles QPSK
    randIdx = floor(rand(nbSymb,1) * length(constelQPSK)) + 1;
    symb = constelQPSK(randIdx);

    %% Suréchantillonnage d'un facteur 8 et filtrage 1/2 Nyquist
    LoverSamp = overSamplingFactor * nbSymb;
    symbOverSampled = zeros(LoverSamp, 1);
    symbOverSampled(1 : overSamplingFactor : end) = symb;
    symbOverSampled = filter(filterHalfNyquist, 1, symbOverSampled); 
    symbOverSampled = symbOverSampled(Lfilter : end);


    %% Ajout de bruit blanc Gaussien fonction du SNR
    snrLin = 10^(snrdB / 10);
    Lsig = length(symbOverSampled);
    noise =  sqrt(overSamplingFactor / (2 * snrLin)) * (randn(Lsig, 1) + 1i * randn(Lsig, 1));
    symbNoisy = symbOverSampled + noise;

    % TODO : ajout d'un perturbateur bande étroite
    symbPSK = create_NBI(48, 0.2, 40);
    requiredpadding = length(symbPSK)-length(symbNoisy);
    symbPadded=[symbNoisy' zeros(1, requiredpadding)];
    symbNoisy = symbPadded' + symbPSK;
    outputsig = filter(Hd, symbNoisy);

    %% Affichage spectre signal + bruit
    nfft = 4096;
    freq = (-nfft/2: (nfft/2-1));
    [symbSpectrum,~]      = hann_spectrogram(symbOverSampled, nfft);
    [symbNoiseSpectrum,~] = hann_spectrogram(symbNoisy, nfft);
    [symbPSKSpectrum,~] = hann_spectrogram(symbPSK, nfft);
    symbSpectrum_dB = 10*log10(symbSpectrum+1e-10);
    symbNoiseSpectrum_dB = 10*log10(symbNoiseSpectrum+1e-10);
    symbPSKSpectrum_dB = 10*log10(symbPSKSpectrum+1e-10);
    figure(2);
    plot(freq, fftshift(symbSpectrum_dB));
	hold on
    plot(freq, fftshift(symbNoiseSpectrum_dB));
    hold on
    plot(freq, fftshift(symbPSKSpectrum_dB));
    hold off
    title("densité spectrale de puissance signal BPSK");
    xlabel("fréquence");
    ylabel("DSP (dB)");
    legend({"signal QPSK seul";"signal QPSK + NBI";"NBI PSK"});

    %% filtrage 1/2 Nyquist à la réception
    symbFiltered = filter(filterHalfNyquist, 1, outputsig);
    symbFiltered  = symbFiltered(Lfilter:end);

    %% Synchronisation pour trouver meilleur instant d'échantillonnage
    powerPeigne = zeros(overSamplingFactor, 1);
    for n = 1:overSamplingFactor
        powerPeigne(n) = mean(abs(symbFiltered(n:overSamplingFactor:end)).^2);
    end
    [~, tSync] = max(powerPeigne);

    %% Décimation à l'instant d'échantillonnage idéal
    symbDecimated = symbFiltered(tSync:overSamplingFactor:end);
    figure(3);
    plot(symbDecimated,"x");
    title("constellation BPSK signal avec bruit")

%     % On corrèle avec les 1000 premiers symboles pour trouver la position de synchronisation
%    corrSync = abs(xcorr(symbDecimated, symb(1:1000))).^2;
%    [~, idxSync] = max(corrSync);
%    idxSync = idxSync - length(symbDecimated) + 1;
%    % if idxSync >= 0
%    %     symbSync = symbDecimated(idxSync:end);
%    %     symbRef = symb(1:length(symbSync));
%    % else
%    %     symbSync = symbDecimated(1:end);
%    %     symbRef = symb((2 - idxSync) : (length(symbSync) + 1 - idxSync));
%    % end
% 
%    if idxSync >= 0
%     symbSync = symbDecimated(idxSync:end);
%     startIdx = max(1, 1);
%     endIdx = min(length(symb), length(symbSync));
%     symbRef = symb(startIdx:endIdx);
%     else
%     symbSync = symbDecimated(1:end);
%     startIdx = max(1, (2 - idxSync));
%     endIdx = min(length(symb), (length(symbSync) + 1 - idxSync));
%     symbRef = symb(startIdx:endIdx);
%     end
% 
% 
%     % Ensure symbSync and symbRef are the same length
%     minLength = min(length(symbSync), length(symbRef));
%     symbSync = symbSync(1:minLength);
%     symbRef = symbRef(1:minLength);
% 
%     % Now calculate SNR output
%     snrOut = -10*log10(mean(abs(symbSync / overSamplingFactor - symbRef).^2));
% 
%    %  % Calcul SNR en sortie de démod
%    % snrOut = -10*log10(mean(abs(symbSync / overSamplingFactor - symbRef).^2))
% 
%     % calcul TEB en supposant que partie réelle et imaginaire portent chacun un bit
%    nbErrorsBits = sum((real(symbSync) .* real(symbRef) < 0) + (imag(symbSync) .* imag(symbRef) < 0));
%    TEB = [TEB, nbErrorsBits / (2 * length(symbRef))];
%     pause(0.1);
% 
end
% 
% figure(4);
% SNR_dB = SNRmin:SNRStep:SNRmax;
% SNR_lin = 10.^(SNR_dB/10);
% semilogy(SNR_dB, TEB);
% hold on
% EbOverN0 = SNR_lin/2;
% TEBtheorique = 0.5*erfc(sqrt(EbOverN0));
% semilogy(SNR_dB, TEBtheorique,"r");
% hold off
% grid
% grid minor
% title ("taux d'erreur binaire vs rapport signal-sur bruit");
% ylabel("TEB")
% xlabel("RSB (dB)")
% legend({"TEB mesuré"; "TEB théorique"})
