
function [y_v, yhann_m] = hann_spectrogram(x_v, nfft)
% Calcul du spectrogramme et du spectre moyen d'un signal
% pour une taille de FFT spécifiée avec utilisation d'une fenêtre de Hann et overlap 50%
%
% INPUT:
%   x_v                   : signal en entrée de longueur L
%   nfft                  : taille de FFT utilisée pour spectrogramme et spectre
%
% OUTPUT:
%   y_v                   : spectre moyen (taille nfft X 1)
%   yhann_m               : spectrogramme, affichable sous forme de waterfall avec pcolor 
%                           (taille nfft X (2*nb_frames) et nb_frames = ceil(L / nfft))
%
% Auteur(s): Paul Gagneur
%%

    %% Division du signal en trames de longueur nfft
    nb_frames = ceil(length(x_v) / nfft);
    Ltot = nb_frames * nfft;
    xfull_v = zeros(Ltot, 1);
    xfull_v(1:length(x_v)) = x_v;
    xhalf_m = reshape(xfull_v, nfft, nb_frames);

    %% Creation des trames intermédiaires (overlap de 50%)
    xfull_m = zeros(nfft, 2*nb_frames);
    xfull_m(:,1:2:end) = xhalf_m;
    xfull_m(1:nfft/2,2:2:end) = xhalf_m((nfft/2+1):nfft, :);
    xfull_m(nfft/2+1:nfft,2:2:end-2) = xhalf_m(1:nfft/2, 2:end);

    %% multiplication des trames par une fenêtre de Hanning
    hann_win = 0.5 *(1-cos(2*pi*(0:nfft -1)/nfft))';
    %%    hann_win = hanning(nfft);
    hann_full_m = repmat(hann_win,1,2*nb_frames);
    xhann_m = xfull_m.*hann_full_m;

    %% calcul du spectrogramme
    yhann_m = abs(fft(xhann_m, nfft)).^2 / (nfft * sum(hann_win.^2));

    %% calcul du spectre moyen
    y_v = mean(yhann_m, 2);
end
