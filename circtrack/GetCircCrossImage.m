function circ_data = GetCircCrossImage(input_vol, vol_size, circ_dia, circ_anum, v_motion)
if nargin == 4
    vx = 0; vy = 0;
else
    vx = real(v_motion); vy = imag(v_motion);
end

circ_data = zeros(size(input_vol, 1), circ_anum);
scan_track = zeros(2,circ_anum);
for i = 1 : circ_anum
    cx = circ_dia * sin(i*2*pi/circ_anum)/2 + vol_size/2 + i *vx;
    cy = -circ_dia * cos(i*2*pi/circ_anum)/2 + vol_size/2 + i *vy;
    cx = cx * (size(input_vol,2)-1)/vol_size +1;
    cy = cy * (size(input_vol,3)-1)/vol_size +1;
    
    cx1 = floor(cx); cx2 = ceil(cx);
    cy1 = floor(cy); cy2 = ceil(cy);
    dx1= cx-cx1; dx2= cx2-cx;
    dy1= cy-cy1; dy2= cy2-cy;
    circ_data(:,i) = ...
        input_vol(:, cx1, cy1) * dx2 * dy2 + ...
        input_vol(:, cx1, cy2) * dx2 * dy1 + ...
        input_vol(:, cx2, cy1) * dx1 * dy2 + ...
        input_vol(:, cx2, cy2) * dx1 * dy1;
    scan_track(1,i) = cx; scan_track(2,i) = cy;
end
% figure;plot(scan_track(1,:), scan_track(2, :));




