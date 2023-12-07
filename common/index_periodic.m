function idx = index_periodic(i, N)
    idx = mod(i-1,N) + 1;
end