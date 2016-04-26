function display_range = GetDisplayRange(image_data, pct)

if (size(image_data, 2) > 1)
    image_data = sort(image_data(:));
else
    image_data = sort(image_data);
end

if (numel(pct) == 2)
    pct_lo = pct(1);
    pct_hi = pct(2);
else
    pct_lo = pct(1);
    pct_hi = pct(1);
end

offset_lo = round(pct_lo * length(image_data));
offset_hi = round(pct_hi * length(image_data));

display_range = [image_data(max(offset_lo, 1)) image_data(max(length(image_data) - offset_hi), 1)];
 