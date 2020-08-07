function I = compute_intensity(E)

I = gather(sum(E.*conj(E)));