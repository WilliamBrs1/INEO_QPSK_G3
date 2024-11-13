function symbPSK = create_NBI(overSamplingFactor, rolloffNyquist, nbSymbHalfNyq)


%Pour Octave : installer éventuellement le package octave-signal
%isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
%if isOctave
%    pkg load signal
%end

constelQPSK = [1; 1i; -1; -1i] * exp(1i * pi / 4);

% génération filtre demi-Nyquist
filterHalfNyquist = CalculHalfNyquist(overSamplingFactor, rolloffNyquist, nbSymbHalfNyq);
Lfilter = length(filterHalfNyquist);

nbSymb = 10000;
snrdB = 10;

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
    symbPSK = symbOverSampled ;