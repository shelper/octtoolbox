function win_func = hann_window(n)

win_func = (0.5 * (1 - cos((2 * pi * (0:(n - 1))) / (n - 1))));
