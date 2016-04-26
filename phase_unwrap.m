function dop_img = phase_unwrap(dop_img, flow_cntr, flow_rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subflow_cntr is the center pixel position (x, y) of the subflow
% flow_rad is the radius of the subflow in pixel
% the function will only treat the pixels within the circle
% some margin for caclulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crop the subflow area from the input image
subflow = dop_img(flow_cntr(1) - flow_rad : flow_cntr(1) + flow_rad, flow_cntr(2) - flow_rad : flow_cntr(2) + flow_rad);
subflow_cntr = [flow_rad+1, flow_rad+1];
% subflow = double(subflow)/256.0*2*pi - pi;

% flt_win =  [0, 1, 0; 1, 1, 1; 0, 1, 0]; %disk shape filter
flt_win =  fspecial('gaussian', 3, 0.95); %disk shape filter
dx = [subflow(:,1),diff(subflow, 1,2)];
figure;imshow(dx)
dx = dx + (dx<-pi)*2*pi; dx = dx - (dx>pi)*2*pi;
figure;imshow(dx,[0,1])
std_dx = sqrt(imfilter(dx.^2, flt_win)  - imfilter(dx, flt_win).^2); %get quality in lateral
figure;imshow(std_dx,[0,1])

dy = [subflow(1,:);diff(subflow, 1,1)];
dy = dy + (dy<-pi)*2*pi; dy = dy - (dy>pi)*2*pi;
figure;imshow(dy,[0,1])
std_dy = abs(sqrt(imfilter(dy.^2, flt_win)  - imfilter(dy, flt_win).^2)); %get quality in axial
figure;imshow(std_dy,[0,1])
quality_map = -(std_dx + std_dy); %get quality map
figure;imshow(quality_map);

%% start quality_map guided unwrapping
% % select start pix from the center of the image
% cntr_size = [10,10];
% cntr_pos = round(size(quality_map)/2);
% cntr_map = quality_map(cntr_pos(1)-cntr_size(1):cntr_pos(1)+cntr_size(1), ...
%                          cntr_pos(2)-cntr_size(2):cntr_pos(2)+cntr_size(2));
% [~,ind] = max(cntr_map(:));
% [row,col] = ind2sub(size(cntr_map), ind);
% ind = sub2ind(size(quality_map), ...
%       row+cntr_pos(1)-cntr_size(1)-1,col+cntr_pos(2)-cntr_size(2)-1);
ind = sub2ind(size(quality_map), flow_rad+1, flow_rad+1);

[row_size,col_size] = size(subflow);
joint_index = [ind];
joint_quality = quality_map(joint_index);
wrapped = ones(size(subflow));

while ~isempty(joint_index)
    [~,i] = max(joint_quality);
    ind = joint_index(i);
%     [row, col]=ind2sub([row_size,col_size], ind);
%     if sqrt((row-subflow_cntr(1))^2 + (col-subflow_cntr(2))^2) < flow_rad %boundary test
    if ind+1 <= row_size * col_size 
        if wrapped(ind+1)
            subflow(ind+1) = subflow(ind) + dy(ind+1);
            joint_index(end+1)=ind+1;
            joint_quality(end+1)=quality_map(ind+1);
        end
    end
    if  ind+row_size <= row_size * col_size 
        if wrapped(ind+row_size)
            subflow(ind+row_size) = subflow(ind) + dx(ind+row_size);
            joint_index(end+1)=ind+row_size;
            joint_quality(end+1)=quality_map(ind+row_size);
        end
    end
    if  ind-1 >= 1 
        if wrapped(ind-1)
            subflow(ind-1) = subflow(ind) - dy(ind);
            joint_index(end+1)=ind-1;
            joint_quality(end+1)=quality_map(ind-1);
        end
    end
    if  ind-row_size >= 1 
        if  wrapped(ind-row_size)
            subflow(ind-row_size) = subflow(ind) - dx(ind);
            joint_index(end+1)=ind-row_size;
            joint_quality(end+1)=quality_map(ind-row_size);
        end
    end 
    % set current pixel from joint pixel as unwrapped
    wrapped(ind) = 0;
    % remove current pixel from joint pixel list
    joint_index(i)=[];
    joint_quality(i)=[];
%     set(h_map, 'CData', wrapped); pause;
end

%% find the surround pixel's subflow value as DC
% background = subflow(:,1:3);
% phase_dc = mean(background(:))
% subflow = subflow-phase_dc;

mean(subflow(:)), max(subflow(:)), min(subflow(:))
dop_img(flow_cntr(1) - flow_rad : flow_cntr(1) + flow_rad, flow_cntr(2) - flow_rad : flow_cntr(2) + flow_rad) =subflow;
% figure;imshow(dop_img,[]);colormap jet;
end





%% additional thoughts:
% 1. for Doppler OCT subflow, center can have higher quality weight
% 2. High intensity pixel can have higher quality weight
% 3. Parallel processing (region growth) from several points
