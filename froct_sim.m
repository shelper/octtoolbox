%% main function to export the following functions
function  funcs = get_oct_sim_funcs()
  funcs.gen_cplx_data=@gen_cplx_data;
  funcs.modulate_phase=@modulate_phase;
  funcs.update_alt_phase_amp=@update_alt_phase_amp;
  funcs.add_flow=@add_flow;
  funcs.phase_alt_deconj=@phase_alt_deconj;
  funcs.linear_mod_deconj=@linear_mod_deconj;
end

%% generate complex interferogram
function gen_cplx_data(pix_num, aline_num, sp_profile, depth_profile, xpos_profile,oversamp_rate)
    k = 0 : 1/(pix_num-1) : 1;
    magnitude = bxsfun(@times, randn(pix_num, aline_num),depth_profile);
    magnitude = bxsfun(@times, magnitude,xpos_profile);
    phase = randn(pix_num, aline_num)*2*pi;
    filter = fspecial('gaussian', oversamp_rate*2, 3);
    magnitude = imfilter(magnitude, filter);
    phase = imfilter(phase, filter);
    cplx_signal = zeros(pix_num, aline_num);
    for x_i = 1:aline_num
        if xpos_profile(x_i)
            for d_i = 1:pix_num
                if depth_profile(d_i)
                    temp = exp(1i * 2 * pi * signal_depth(d_i-pix_num/2) * k)
                    temp = temp * magnitude(d_i, x_i) *exp(1i * phase(d_i, x_i));
                    cplx_signal(:,x_i) = cplx_signal(:,x_i) + temp';
                end
            end
            cplx_signal(:,x_i) = cplx_signal(:,x_i) .* sp_profile;
            x_i
        end
    end
end

function data = modulate_phase(cplx_data, phase_amp, mod_freq)
    [pix_num, aline_num] = size(cplx_data);
    if size(phase_amp, 1)<size(phase_amp, 2)
        phase_amp = phase_amp';
    end
    if size(mod_amp, 1) == 1 || size(mod_amp, 1) == pix_num
        mod_phase = phase_amp * cos((1:aline_num)*2*pi*mod_freq);
    else
        fprintf('mod_amp size does not match pix_num\n')
        return
    end
    data = real(bxsfun(@times, cplx_data, exp(1i*mod_phase));
end

function data = update_alt_phase_amp(data, orig_phase_amp, new_phase_amp)
    [pix_num, aline_num] = size(data);
    H0 = GetFreqComp(data',[1, aline_num/4],3)';
    H0 = bsxfun(@times, H0, tan(orig_alt_phase)./tan(new_alt_phase));
    H1=GetFreqComp(data',[aline_num/4+1,aline_num/2+1],3)';
    data = H0 + H1;
end

function flow_data = add_flow(cplx_data, flow_radius, flow_rate, flow_depth, flow_xpos)
    [pix_num, aline_num] = size(cplx_data);
    flow_phase =fspecial('gaussian', flow_radius*2+1, flow_radius/1.5);
    flow_phase = flow_phase/max(flow_phase(:)) * flow_rate;
    flow_phase = cumsum(flow_phase, 2);
    aline_phase = zeros(pix_num, size(flow_phase, 2));
    aline_phase(flow_depth-flow_radius: flow_depth+flow_radius,:) = flow_phase;
    flow_data = fft(cplx_data);
    if flow_xpos-flow_radius >=1 && flow_xpos+flow_radius <= aline_num
        flow_data(:, flow_xpos-flow_radius:flow_xpos+flow_radius) = flow_data(:, i) .*exp(1i*aline_phase);
    else
        fprintf('flow size exceeds the boundary');
    end
    flow_data = real(ifft(flow_data));
end

function fr_data = phase_alt_deconj(data, phase_amp, bandwidth, order)
    [pix_num, aline_num] = size(cplx_data);
    bandwidth = aline_num * bandwidth
    if order ==1
        H0 = GetFreqComp(data',[1, bandwidth],3)';
        H0 = bsxfun(@times, H0, tan(phase_amp));
        H1 = data.*repmat([1,-1],[pix_num,aline_num/2]);%figure;plot(abs(fft(NewData,[],2))')
        H1 = GetFreqComp(H1',[1, bandwidth],3)';
        fr_data = H0 - 1i*H1 ;
    elseif order ==2
        cplx_img = fft(data);
        phase_diff = angle(cplx_img);
        phase_diff(:, 2:end)= angle(cplx_img(:, 2:end) .* conj(cplx_img(:, 1:end-1)));
        cplx_img = abs(cplx_img).* exp(1i*phase_diff);% 2nd order signal generation
        H0 = GetFreqComp(cplx_img',[1, bandwidth],3)';
        H0 = bsxfun(@times, H0, tan(phase_amp/2));
        H1 = cplx_img.*repmat([1,-1],[pix_num,aline_num/2]);
        H1 = GetFreqComp(H1',[1, bandwidth],3)';
        cplx_img = H0 - 1i*H1;
        fr_data = ifft(cplx_img);
    end
end

function fr_data = linear_mod_deconj(data, phase_amp, order)
    [pix_num, aline_num] = size(cplx_data);
    if order ==1
        data = update_alt_phase_amp(data, phase_amp, pi/4);
        data(:, 3:4:end) = -linear_signal(:, 3:4:end);
        data(:, 4:4:end) = -linear_signal(:, 4:4:end);
        fr_data = GetFreqComp(data',[1, aline_num/2],1)';
    elseif order == 2
        cplx_img = fft(data);
        phase_diff = angle(cplx_img);
        phase_diff(:, 2:end)= angle(cplx_img(:, 2:end) .* conj(cplx_img(:, 1:end-1)));
        cplx_img = abs(cplx_img).* exp(1i*phase_diff);% 2nd order signal generation
        cplx_img(:, 3:4:end) = -cplx_img(:, 3:4:end);
        cplx_img(:, 4:4:end) = -cplx_img(:, 4:4:end);
        cplx_img = GetFreqComp(cplx_img',[1, aline_num/2],1)';
        fr_data = ifft(cplx_img);
    end
end
