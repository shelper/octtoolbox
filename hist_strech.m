function img_stretched = hist_strech(img, peak_l, peak_h, strch)
img = double(img);
peak = (peak_l + peak_h)/2;
width = peak_h - peak_l;
peak_ll = peak - width*strch/2;
peak_hh = peak + width*strch/2;

img_l = img .* (img< peak_l);
img_h = img .* (img> peak_h);
img_m = img -img_l-img_h;

img_m = img_m * strch - peak * (strch-1);
img_l = img_l*peak_ll/peak_l;
img_h = (img_h + (peak_hh-peak_h)) * (65535-peak_hh)/(65535-peak_h);

img_stretched = uint16(img_h)+uint16(img_l)+uint16(img_m);
