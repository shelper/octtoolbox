# -*- coding: utf-8 -*-
from numpy import *

# ### read in the raw data from file
def raw_data_fromfile(raw_data_file, oct_model, dtype='uint16'):
    '''
    read oct raw data, if 2 out put is needed, the 2nd is calibroate_ coeff
    if oct_model is a list, it can be [aline_num, pix_num]
    the raw_data is transposed to pixel to be vertical & aline to be horizontal
    '''
    raw_data = fromfile(raw_data_file, dtype)
    if oct_model == 'oct-2000':
        pix_num = int(raw_data[0])
        aline_num = int(raw_data[2])
        ## calibration data
        calibrate_coeff = raw_data[-(pix_num*6):-(pix_num*4)]
        calibrate_coeff.dtype = float32
        ## dc of reference from the raw data file
        refrence = raw_data[-(pix_num*4):]
        reference.dtype = double
        ## reshape raw data
        raw_data = raw_data[6:6 + pix_num * aline_num].reshape([-1, pix_num])
        raw_data = raw_data.astype('float32')
        raw_data -= reference
        raw_data = raw_data.transpose()
        return raw_data, calibrate_coeff
    elif oct_model == 'dri-oct1':
        pix_num = int(raw_data[0])
        aline_num = int(raw_data[2])
        raw_data = raw_data[6:6 + pix_num * aline_num]
        raw_data = raw_data.reshape([-1, pix_num])
        raw_data = raw_data.astype('float32')
        raw_data -= 2048
        raw_data = raw_data.transpose()
        return raw_data
    # oct_model as the shape should be [aline_num, pix_num]
    elif type(oct_model) is list and size(oct_model) == 2:
        raw_data = raw_data.reshape(oct_model)
        raw_data = raw_data.astype('float32')
        raw_data = raw_data.transpose()
        return raw_data
    else:
        print("file read error, cannot read the oct data")


# ### calibrate the sdOCT data
def calibrate_raw_data(raw_data, calibrate_coeff, method='linear'):
    from scipy.interpolate import interp1d
    pix_num = raw_data.shape[1]
    linear_k = linspace(0,1,pix_num)
    calibrated_data = interp1d(linear_k, raw_data, method)(calibrate_coeff)
    calibrated_data -= mean(calibrated_data, axis = 0)
    return calibrated_data

# ### calculate the DC reference by different approaches
def get_reference_by_cplx_median(raw_data):
    cplx_img = fft.fft(raw_data, axis = 0)
    cplx_reference = 1j* median(cplx_img.imag, 1) + median(cplx_img.real, 1)
    reference = fft.ifft(cplx_reference)
    return reference.real

# ### dispersion compensation
def eval_dispersion_coeff(dispersion_coeff, raw_data, method='entropy'):
    pix_num = raw_data.shape[0]
    #import pdb; pdb.set_trace()
    dispersion_coeff = concatenate((dispersion_coeff,[0,0]))
    dispersion_phase = polyval(dispersion_coeff, linspace(0, 1, pix_num))
    dispersion_phase -= linspace(dispersion_phase[0], dispersion_phase[-1],
                                 pix_num)
    #from matplotlib.pyplot import *
    #figure();plot(dispersion_phase)
    raw_data = raw_data * exp(-1j * dispersion_phase).reshape([pix_num, -1])
    img = reconstruct_img(raw_data)
    if method == 'entropy':
        img = img[50:,:]
        img = img / sum(img);
        print(-sum(img * log(img)))
        return -sum(img * log(img))
    elif method == 'variance':
        return sum((img-mean(img))**2)/mean(img)
    elif method == 'sum':
        return sum(img)
    elif method == 'logsum':
        return sum(log(img))
    else:
        print('method %s not avaliable' %method)
        return -1


# ### oct images reconstruction
def reconstruct_img(raw_data, fft_num=None,
                    img_range='half', img_kind='real', img_scale='default'):
    if fft_num is None:
        img = fft.fft(raw_data, axis = 0)
        fft_num = img.shape[0]
    else:
        img = fft.fft(raw_data, fft_num, axis = 0)
    if img_range == 'half':
        img = img[0:fft_num/2+1,:]
    if img_kind == 'real':
        if img_scale == 'default' or img_scale == 'log':
            img = 20*log10(abs(img)+1)
        else:
            img = abs(img)
    return img


# ### Doppler images reconstruction
def get_kasai_doppler_img(cplx_img, filter_size):
    doppler_img = cplx_img[:, 1:] * conjugate(cplx_img[:, 0:-1])
    if filter_size > 1 :
        from scipy.ndimage.filters import uniform_filter
        doppler_img = uniform_filter(doppler_img.real, filter_size) + \
                      1j * uniform_filter(doppler_img.imag, filter_size)
    doppler_img = angle(doppler_img)
    return doppler_img


def get_limp_doppler_img(cplx_img, filter_size, limp_ratio=2):
    aline_num = cplx_img.shape[1]
    oct_img = abs(cplx_img)
    limp_filter = zeros(aline_num)
    limp_filter[0:aline_num/2+1] = linspace(0, 1, aline_num/2+1)
    limp_img_pos = abs(fft.ifft(fft.fft(cplx_img, axis=1) \
                                * limp_filter, axis=1))

    limp_filter = flipud(limp_filter)
    limp_img_neg = abs(fft.ifft(fft.fft(cplx_img, axis=1) \
                                * limp_filter, axis=1))
    if filter_size > 1:
        from scipy.ndimage.filters import uniform_filter
        doppler_img_pos = uniform_filter(limp_img_pos) \
                          / uniform_filter(oct_img)
        doppler_img_neg = uniform_filter(limp_img_neg) \
                          / uniform_filter(oct_img)
    else:
        doppler_img_pos = limp_img_pos / oct_img
        doppler_img_neg = limp_img_neg / oct_img

    ## set pos/neg ratio threshold
    mask_pos = doppler_img_pos > doppler_img_neg \
               * (limp_ratio - doppler_img_neg*(limp_ratio - 1))
    mask_neg = doppler_img_neg > doppler_img_pos \
               * (limp_ratio - doppler_img_pos*(limp_ratio - 1))
    doppler_img = doppler_img_pos * mask_pos - doppler_img_neg * mask_neg

    doppler_img *= pi * ((sqrt(2) - 1) * abs(doppler_img) + 1)
    doppler_img[doppler_img > pi] = pi
    doppler_img[doppler_img < -pi] = -pi
    return doppler_img


# ### OCT micro-angiography (OMAG)
# 1. phase variance
# 2. intentsity variance
# 3. aline amplitude differences
# 4. aline complex value differences
def get_phase_variance_omag(cplx_img, filter_size, intensity_thresh):
    conj_img = cplx_img[:, 1:] * cplx_img[:, 0:-1]
    from scipy.ndimage.filters import uniform_filter
    conj_img = uniform_filter(conj_img.real, filter_size) + \
               1j * uniform_filter(conj_img.imag, filter_size)
    cc_img = abs(conj_img) * 2
    amplitude_img = abs(cplx_img)**2
    ac_img = amplitude_img[:, 0:-1] + amplitude_img[:, 1:]
    ac_img = uniform_filter(ac_img, filter_size)
    omag_img = 1-cc_img/ac_img
    omag_img *= (ac_img > intensity_thresh**2)
    return omag_img


def get_intensity_variance_omag(amplitude_img, filter_size, intensity_thresh):
    cc_img = 2*amplitude_img[:, 0:-1] * amplitude_img[:, 1:]
    from scipy.ndimage.filters import uniform_filter
    cc_img = uniform_filter(cc_img, filter_size)
    amplitude_img = amplitude_img**2
    ac_img = amplitude_img[:, 0:-1] + amplitude_img[:, 1:]
    ac_img = uniform_filter(ac_img, filter_size)
    omag_img = 1-cc_img/ac_img
    omag_img *= (ac_img > intensity_thresh**2)
    return omag_img


def get_diff_omag(amplitude_img, filter_size, intensity_thresh):
    diff_img = abs(diff(amplitude_img, axis=1))
    mean_img = (amplitude_img[:, 1:] + amplitude_img[:, 0:-1])/2

    from scipy.ndimage.filters import uniform_filter
    diff_img = uniform_filter(diff_img, filter_size)
    mean_img = uniform_filter(mean_img, filter_size)

    omag_img = diff_img / mean_img
    omag_img *= (mean_img > intensity_thresh)
    return omag_img


def get_cplx_diff_omag(cplx_img, filter_size, intensity_thresh):
    diff_img = abs(diff(cplx_img, axis=1))
    mean_img = abs(cplx_img[:, 1:] + cplx_img[:, 0:-1])/2

    from scipy.ndimage.filters import uniform_filter
    diff_img = uniform_filter(diff_img, filter_size)
    mean_img = uniform_filter(mean_img, filter_size)

    omag_img = diff_img / mean_img
    omag_img *= (mean_img > intensity_thresh)
    return omag_img


# ### correct bulky motion (crucial before applying LIMP method or diff method)
# 1. std of amplitude lines' diff to identify flows
# 2. remove flow pixels and calculate average phase differences as bulky motion
def correct_bulk_motion(cplx_img, bulk_mask):
    bulk_img = cplx_img[:, 1:] * conjugate(cplx_img[:, 0:-1]) * bulk_mask
    bulk_phase = angle(sum(bulk_img, axis = 0))
    bulk_phase = cumsum(bulk_phase)

    cplx_img[:,1:] = cplx_img[:,1:] * exp(-1j*bulk_phase)
    return cplx_img


# ### correct phase jittering by FBG signature spike
def correct_phase_jitter_by_FBG(raw_data, spike_range):
    aline_num = raw_data.shape[1]
    spike_data = raw_data[spike_range[0]: spike_range[1], :]
    spike_idx = abs(spike_data).argmax(axis=0)
    spike_idx -= spike_idx[0]
    # spike_range = arange(spike_range[0], spike_range[1])
    # spike_data = raw_data[spike_range, :]
    for i in xrange(1, aline_num):
        if spike_idx[i] != 0:
            raw_data[:, i] = roll(raw_data[:, i], -spike_idx[i])

    return raw_data


# ### correct phase jittering in swept source by minimize phase diff
# 1. substract each line with averaged phase (bulk_phase)
# 2. find out the difference of compensated complex alines
# 3. using binary search to reduce the loop number
def binary_correct_phase_jitter(raw_data, max_shift=2, mask=1):
    cplx_img = reconstruct_img(raw_data, img_range='half', img_kind='cplx')
    img_depth, aline_num = cplx_img.shape
    pixes_shift = zeros(aline_num-1, dtype=int64)
    jitter_phase = linspace(0, pi, img_depth).reshape([-1,1])

    for i in xrange(int(ceil(log2(aline_num)))):
        alines_left = arange(2**i - 1, aline_num-1, 2**(i+1))
        alines_right = alines_left + 1
        sub_img_left = cplx_img[:, alines_left]
        jitter_err = zeros([max_shift*2 + 1, alines_left.size])
        if mask is not 1:
            sub_mask = mask[:, alines_left]
        else:
            sub_mask = 1

        for n in xrange(max_shift*2+1):
            jitter_gain = exp(1j * (max_shift-n) * jitter_phase)
            sub_img_right = cplx_img[:, alines_right] * jitter_gain
            conj_img = sub_img_right * conjugate(sub_img_left) * sub_mask
            bulk_phase = angle(sum(conj_img, axis = 0))
            sub_img_diff = sub_img_right * exp(-1j * bulk_phase) - sub_img_left
            jitter_err[n, :] = sum(abs(sub_img_diff * sub_mask), axis = 0)

        pix_shift = max_shift - jitter_err.argmin(axis=0)
        for n in xrange(alines_left.size):
            ind1 = alines_left[n] - (2 ** i - 1)
            ind2 = alines_right[n]
            raw_data[:, ind1:ind2] = roll(raw_data[:, ind1:ind2],
                                          pix_shift[n], axis = 0)
        cplx_img = reconstruct_img(raw_data, img_range='half', img_kind='cplx')
        pixes_shift[alines_left] = pixes_shift[alines_left] + pix_shift

    return raw_data


# ### correct jittering line by line
def correct_phase_jitter(raw_data, max_shift=2, mask=1):
    aline_num = raw_data.shape[1]
    jitter_err = zeros(max_shift * 2 + 1)
    pixes_shift = zeros(aline_num, dtype=int)

    # correct the first alines
    aline1 = fft.rfft(raw_data[:, 0])
    if mask is not 1:
        line_mask = mask[:, 0]
    else:
        line_mask = 1
    for n in xrange(max_shift*2+1):
        aline2 = fft.rfft(roll(raw_data[:, 1], n-max_shift))
        conj_aline = aline1 * conjugate(aline2) * line_mask
        bulk_phase = angle(sum(conj_aline))
        aline_diff = aline2 * exp(1j * bulk_phase) - aline1
        jitter_err[n] = sum(abs(aline_diff * line_mask))
    pix_shift = jitter_err.argmin() - max_shift
    raw_data[:, 1] = roll(raw_data[:,1], pix_shift)

    for i in xrange(2, aline_num):
        aline1 = fft.rfft(raw_data[:, i-2])
        aline2 = fft.rfft(raw_data[:, i-1])
        if mask is not 1:
            line_mask = mask[:, i-1]
        else:
            line_mask = 1

        for n in xrange(max_shift*2+1):
            aline3 = fft.rfft(roll(raw_data[:, i], n-max_shift))
            phase_aline = angle(sum(aline1 * conjugate(aline2) * line_mask))
            #phase_aline = angle(aline1 * conjugate(aline2))
            line_diff = aline3 * exp(1j * phase_aline) - aline2
            jitter_err[n] = sum(abs(line_diff * line_mask), axis = 0)

        pix_shift = jitter_err.argmin() - max_shift
        raw_data[:, i] = roll(raw_data[:,i], pix_shift)

        pixes_shift[i] = pix_shift
    #figure();plot(pix_shift), plot(cumsum(pix_shift))
    return raw_data

# ### test the functions above
if __name__ == 'main':

    # ## test oct functions on oct-2000 model
    raw_data_file = 'D:\\Topcon\\Projects\\2010-10-28 New Doppler algorithm' \
                    '\\dataset\\pattern1-3460\\3489.RAW'
    raw_data, calibrate_coeff = raw_data_fromfile(raw_data_file, 'oct-2000')
    calibrated_data = calibrate_raw_data(raw_data, calibrate_coeff, 'linear')
    calibrated_data = calibrated_data.transpose()
    cplx_img = reconstruct_img(calibrated_data, img_kind='cplx')
    #amplitude_img = 20*log10(abs(cplx_img)+1)
    #figure();imshow(amplitude_img[100:600,:])
    amplitude_img = abs(cplx_img)
    bulk_mask = get_intensity_variance_omag(amplitude_img, 3, 900000)
    bulk_mask = (bulk_mask < 0.05) * (bulk_mask > 0)
    figure();imshow(bulk_mask[100:600,:])
    #cplx_img = correct_bulk_motion(cplx_img, bulk_mask)
    #bulk_mask = get_cplx_diff_omag(cplx_img, 2, 200)
    #bulk_mask = (bulk_mask < 0.5) * (bulk_mask > 0)
    #figure();imshow(bulk_mask[100:600,:])
    ## comparison of 4 different angiography methods
    #omag_img = get_intensity_variance_omag(amplitude_img, 5, 300)
    #figure();imshow(omag_img[100:600,:])
    #omag_img = get_diff_omag(amplitude_img, 5, 200)
    #figure();imshow(omag_img[100:600,:])
    #omag_img = get_phase_variance_omag(cplx_img, 5, 300)
    #figure();imshow(omag_img[100:600,:])
    #omag_img = get_cplx_diff_omag(cplx_img, 5, 200)
    #figure();imshow(omag_img[100:600,:])
    #doppler_img = get_kasai_doppler_img(cplx_img, 3)
    #figure();imshow(doppler_img[100:600,:])
    #doppler_img = get_kasai_doppler_img(diff(cplx_img, axis=1), 3)
    #figure();imshow(doppler_img[100:600,:])
    #doppler_img = get_limp_doppler_img(cplx_img, 3)
    #figure();imshow(doppler_img[100:600,:])
    #doppler_img = get_limp_doppler_img(diff(cplx_img, axis=1), 3)
    #figure();imshow(doppler_img[100:600,:])
    cplx_img = correct_bulk_motion(cplx_img, bulk_mask)
    doppler_img = get_kasai_doppler_img(cplx_img, 3)
    figure();imshow(doppler_img[100:600,:])
    doppler_img = get_kasai_doppler_img(diff(cplx_img, axis=1), 3)
    figure();imshow(doppler_img[100:600,:])
    doppler_img = get_limp_doppler_img(cplx_img, 3, 1.5)
    figure();imshow(doppler_img[100:600,:])
    doppler_img = get_limp_doppler_img(diff(cplx_img, axis=1), 3, 2)
    figure();imshow(doppler_img[100:600,:])

    # ## test oct functions on dri-oct1
    raw_data_file = 'D:\\Topcon\\Projects\\2010-10-28 New Doppler Algorithm' \
                    '\\Flow\\Octraw0002.Raw'
    raw_data = raw_data_fromfile(raw_data_file, 'dri-oct1')
    pix_num = raw_data.shape[1]
    raw_data -= mean(raw_data, axis=0)
    from scipy.signal import hann
    raw_data = raw_data * hann(pix_num)
    raw_data = raw_data.transpose()
    raw_data = raw_data[:, ::4]

    # correct dispersion
    from scipy.optimize import fmin as fminsearch
    dispersion_coeff = fminsearch(eval_dispersion_coeff,
                            [0, 0], [raw_data], xtol = 0.1)
    dispersion_coeff = concatenate((dispersion_coeff, [0,0]))
    dispersion_phase = polyval(dispersion_coeff, linspace(0, 1, pix_num))
    dispersion_phase -= linspace(dispersion_phase[0], dispersion_phase[-1],
                                 pix_num)
    raw_data *= exp(-1j * dispersion_phase).reshape([pix_num, -1])

    img = reconstruct_img(raw_data)
    figure();imshow(20*log10(img))
    #cplx_img = reconstruct_img(raw_data, img_range='half', img_kind='cplx')
    #doppler_img = get_kasai_doppler_img(cplx_img, 3)
    #figure();imshow(doppler_img)
    amplitude_img = reconstruct_img(raw_data, img_kind='real')
    #figure();imshow((amplitude_img))
    bulk_mask = get_intensity_variance_omag(amplitude_img, 5, 900000)
    #figure();imshow(bulk_mask)
    bulk_mask = (bulk_mask < 0.15) * (bulk_mask > 0)
    #figure();imshow(bulk_mask)
    bulk_mask[500:,:]=0
    raw_data = correct_phase_jitter(raw_data, 3, bulk_mask)
    cplx_img = reconstruct_img(raw_data, img_range='half', img_kind='cplx')
    doppler_img = get_kasai_doppler_img(cplx_img, 3)
    figure();imshow(doppler_img)
