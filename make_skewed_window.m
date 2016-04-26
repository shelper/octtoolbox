function win = make_skewed_window(winfunc, winfunc2, win_length, win_center_pct)

if (nargin < 4)
    win_center_pct = win_length;
    win_length = winfunc2;
    winfunc2 = winfunc;
end

win_center = round(win_length * win_center_pct);
win = zeros(win_length, 1);

len1 = ((win_center * 2) - 1);
win1 = winfunc(len1);
win(1:win_center) = win1(1:win_center);

len2 = (((win_length - win_center) * 2) + 1);
win2 = winfunc2(len2);
win(win_center:win_length) = win2((end - (win_length - win_center)):end);
