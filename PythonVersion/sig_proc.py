'''
Created on Jun 8, 2012

@author: ZYUAN
'''
import numpy as np
from scipy import fftpack
#import matplotlib.pyplot as plt

#def plot_data_curve(data, draw_row = 1, row_slice = [], col_slice = []):
#    '''
#    function plot_data_curve(data, draw_row = 1, row_slice = [], col_slice = [])
#    plot the data with the following parameters:
#    data: input data
#    row_slice | col_slice: aray or list to retrieve part of the data
#    '''
#    if row_slice:
#        data = data[row_slice,:]
#    if col_slice:
#        data = data[:,col_slice]
#    if draw_row:
#        return plt.plot(data.T)
#    else:
#        return plt.plot(data)

def fft_filt_pass(data, filt_win, side = 0, real_output = True):
    '''
    function fft_filt_pass(data, filt_win, side = 0, real_output = True)
    filter the data with the following parameters:
    data: input data
    filt_win: pixel ind or normalized windnow within range [0, 1], where 1 is the half column size
    side: -1 for negative side, +1 for positive side, 0 for both side
    real_output: True if return real number, otherwise complex
    '''
    data_fft = fftpack.fft(data)
    filted_fft = np.zeros_like(data_fft,'complex64')
    if all(0<=foo<=1 for foo in filt_win):
        filt_win = np.int_(data_fft.shape[-1]  * np.array(filt_win)/2) #convert to integer
    else:
        filt_win = np.array(filt_win)

    if (len(data.shape) == 1):
        if side > 0:
            filted_fft[filt_win[0]:filt_win[1]]= data_fft[filt_win[0]:filt_win[1]]
        elif side < 0:
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            filted_fft[filt_win[0]:filt_win[1]]= data_fft[filt_win[0]:filt_win[1]]
        else:
            filted_fft[filt_win[0]:filt_win[1]]= data_fft[filt_win[0]:filt_win[1]]
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            filted_fft[filt_win[1]:filt_win[0]]= data_fft[filt_win[1]:filt_win[0]]
    else:
        if side > 0:
            filted_fft[:, filt_win[0]:filt_win[1]]= data_fft[:, filt_win[0]:filt_win[1]]
        elif side < 0:
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            filted_fft[:, filt_win[0]:filt_win[1]]= data_fft[:, filt_win[0]:filt_win[1]]
        else:
            filted_fft[:, filt_win[0]:filt_win[1]]= data_fft[:, filt_win[0]:filt_win[1]]
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            filted_fft[:, filt_win[1]:filt_win[0]]= data_fft[:, filt_win[1]:filt_win[0]]
    return (fftpack.ifft(filted_fft)).real if real_output else fftpack.ifft(filted_fft)


def fft_filt_block(data, filt_win, side = 0, real_output = True):
    '''
    function fft_filt_block(data, filt_win, side = 0, real_output = True)
    filter the data with the following parameters:
    data: input data
    filt_win: pixel ind or normalized windnow within range [0, 1], where 1 is the half column size
    side: -1 for negative side, +1 for positive side, 0 for both side
    real_output: True if return real number, otherwise complex
    '''
    data_fft = fftpack.fft(data)
    if all(0<=foo<=1 for foo in filt_win):
        filt_win = np.int_(data_fft.shape[-1]  * np.array(filt_win)/2) #convert to integer
    else:
        filt_win = np.array(filt_win)

    if (len(data.shape) == 1):
        if side > 0:
            data_fft[filt_win[0]:filt_win[1]]= 0
        elif side < 0:
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            data_fft[filt_win[0]:filt_win[1]]= 0
        else:
            data_fft[filt_win[0]:filt_win[1]]= 0
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            data_fft[filt_win[1]:filt_win[0]]= 0
    else:
        if side > 0:
            data_fft[:, filt_win[0]:filt_win[1]]= 0
        elif side < 0:
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            data_fft[:, filt_win[0]:filt_win[1]]= 0
        else:
            data_fft[:, filt_win[0]:filt_win[1]]= 0
            filt_win = (data_fft.shape[-1] + 1)- filt_win
            data_fft[:, filt_win[1]:filt_win[0]]= 0
    return (fftpack.ifft(data_fft)).real if real_output else fftpack.ifft(data_fft)


def hilbert(data):
    ''' function hilbert(data) returns the analytic signal of the input data
    '''
    data_fft = fftpack.fft(data)
    if data.shape[-1] % 2:
        heavi_win = [2]*(data.shape[-1]/2) + [1] + [0] * (data.shape[-1]/2)
    else:
        heavi_win = [2]*(data.shape[-1]/2) + [0] * (data.shape[-1]/2)
    heavi_win[0]=1;

    return fftpack.ifft(data_fft * heavi_win)
