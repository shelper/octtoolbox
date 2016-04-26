from sys import path as load_path
path_to_oct_tools = r'C:\Users\zyuan\Dropbox\research\OCTToolbox\PythonVersion'
load_path.append(path_to_oct_tools)
import python_oct as oct
raw_data_folder = r'D:\Topcon\Projects\2013-10-04 ssOCT doppler\RAW\0006'
raw_data_file = '\\'.join([raw_data_folder, 'octimg' + '%04d' %0 + '.raw'])
raw_data = oct.raw_data_fromfile(raw_data_file, 'dri-oct1')
pix_num = raw_data.shape[1]
from scipy.signal import hann
raw_data = raw_data * hann(pix_num)
raw_data = raw_data.transpose()
# raw_data = raw_data[:, ::4]

# #correct dispersion
# for i in xrange(-20, 20):
#     dispersion_coeff = [0, i * 20.0, 0, 0]
#     dispersion_phase = polyval(dispersion_coeff, linspace(0, 1, pix_num))
#     dispersion_phase -= linspace(dispersion_phase[0], dispersion_phase[-1],
#                                  pix_num)
#     oct_data = raw_data * exp(-1j * dispersion_phase).reshape([pix_num, -1])
#     img = oct.reconstruct_img(oct_data)  
#     img = img[50:, :]
#     print(dispersion_coeff, sum(img), sum(20*log10(img)))

## dispersion compensation
from scipy.optimize import fmin as fminsearch
dispersion_coeff = [0, 228, 0, 0]
dispersion_phase = polyval(dispersion_coeff, linspace(0, 1, pix_num))
dispersion_phase -= linspace(dispersion_phase[0], dispersion_phase[-1],
                             pix_num) 
raw_data = real(raw_data * exp(-1j * dispersion_phase).reshape([pix_num, -1]))  
raw_data -= mean(raw_data, axis=0) # dc removal
img = oct.reconstruct_img(raw_data)   
figure();imshow(20*log10(img));     
