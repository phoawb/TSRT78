function E = getEnergy(s, domain)
    %getEnergy get the energy of the signal
    N = length(s);

    switch domain
        case "time"
            E = sum(s .^ 2);
        case "frequency"
            S = fft(s);
            E = (1 / N) * sum(abs(S) .^ 2);
        otherwise
            error("Domain must be 'time' or 'frequency'");
    end

end
