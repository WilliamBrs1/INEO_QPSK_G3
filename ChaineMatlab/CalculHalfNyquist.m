
function filterHalfNyquist = CalculHalfNyquist(NbEcht_Sym, alpha, nbSymbHalfNyq)
    % calcul de la fonction de mise en forme
    % NbEcht_Sym est le nombre d'�chantillons par p�riode symbole
    %
    L = nbSymbHalfNyq / 2; % nombre de symboles avant et apr�s sur lequel s'�tend le filtre
    % l'�tendue du filtre est donc 2*L.

    hd = zeros(2 * L * NbEcht_Sym + 1, 1);
    t = (-L * NbEcht_Sym : L * NbEcht_Sym) / NbEcht_Sym;

    %formule issue de https://en.wikipedia.org/wiki/Root-raised-cosine_filter
    hd =  (sin(pi*t*(1-alpha)) + 4*alpha*t .* cos(pi*t*(1+alpha))) ./ (pi*t.*(1-(4*alpha*t).^2));

    hd(L*NbEcht_Sym+1) = (1+alpha*(4/pi-1));
    nplus = round(NbEcht_Sym/(4*alpha));
    if abs(nplus -(NbEcht_Sym/(4*alpha))) < 1e-5
        val = alpha / sqrt(2)*((1 + 2/pi) * sin(pi/(4*alpha)) + (1 - 2/pi) * cos(pi/(4*alpha)));
        hd(L*NbEcht_Sym+1+nplus) = val;
        hd(L*NbEcht_Sym+1-nplus) = val;
    end
    %filterHalfNyquist = sqrt(NbEcht_Sym) * hd / sum(hd);
    filterHalfNyquist = NbEcht_Sym * hd / sum(hd);
    %filterHalfNyquist = hd / sum(hd);
end

