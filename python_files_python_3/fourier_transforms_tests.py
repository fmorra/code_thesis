import numpy as np
import scipy.signal as signal
import pdb
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
# Script to test different signal functions and the Fourier transform function from Python

def custom_sawtooth_1(x_1, period_1, offset_1):
    return ((x_1 - offset_1) % period_1) + offset_1


def custom_sawtooth_2(x_2, period_2, offset_2):
    interval_1 = (((-1/8)*period_2 + np.abs(offset_2) > ((x_2 - offset_2) % period_2) + offset_2) & (((x_2 - offset_2) % period_2) +
                                                                                             offset_2 >= (-5/8)*period_2 - np.abs(offset_2)))
    interval_2 = (((3/8)*period_2 + np.abs(offset_2) > ((x_2 - offset_2) % period_2) + offset_2) & (((x_2 - offset_2) % period_2) +
                                                                                 offset_2 >= (-1/8)*period_2-np.abs(offset_2)))
    function_values = []
    for i in range(len(x_2)):
        if interval_1[i]:
            function_values.append(1/4 + ((x_2[i] - offset_2) % period_2) + offset_2)
        elif interval_2[i]:
            function_values.append(5/4 + ((x_2[i] - offset_2) % period_2) + offset_2)
        else:
            function_values.append(0)
    return function_values

    # return np.piecewise(x_2, [interval_1, interval_2], [lambda x: 1/4 + ((x_2 - offset_2) % period_2) + offset_2,
    #                                                     lambda x: 5/4 + ((x_2 - offset_2) % period_2) + offset_2])


T_1 = 0.2
freq_1 = 1/T_1
t_1 = np.linspace(-3, 3, 300)
function_1 = signal.square((2 * np.pi * freq_1 * t_1) + 1)
transform_1 = fft(function_1)

# plt.plot(t_1, function_1)
# plt.plot(t_1, np.abs(transform_1))
# plt.xlabel('t')
# plt.ylabel('Amplitude')
# plt.grid()
# plt.show()


T_2 = 1
freq_2 = 1/T_2
t_2 = np.linspace(-3, 3, 300)
function_2 = signal.sawtooth((2 * np.pi * freq_2 * t_2), 0.5)
transform_2 = fft(function_2)

# plt.plot(t_2, function_2)
# plt.plot(t_2, np.abs(transform_2))
# plt.xlabel('t')
# plt.ylabel('Amplitude')
# plt.grid()
# plt.show()

T_coeff_1 = 2
T_3 = T_coeff_1*np.pi
offset_coefficient_1 = -(T_coeff_1/2)
normalizer_1 = 1  # 1.5*T_coeff_1
offset_1 = offset_coefficient_1*np.pi  # Centered on the X-axis too
t_3 = np.linspace(-3*T_3, 3*T_3, 300)
function_3 = custom_sawtooth_1(t_3, T_3, offset_1)
transform_3 = fft(function_3)

# plt.plot(t_3, function_3/normalizer_1)
# # plt.plot(t_3, np.abs(transform_3))
# plt.xlabel('t')
# plt.ylabel('Amplitude')
# plt.grid()
# plt.show()

T_coeff_2 = 2
T_4 = T_coeff_2*np.pi
offset_coefficient_2 = -(T_coeff_2/3)
normalizer_2 = 1.5*T_coeff_2
offset_2 = offset_coefficient_2*np.pi  # Centered on the X-axis too
t_4 = np.linspace(-3*T_4, 3*T_4, 300)
function_4 = custom_sawtooth_2(t_4, T_4, offset_2)
transform_4 = fft(function_4)
transform_4 = transform_4[30:-30]
plt.plot(t_4, function_4)
# plt.plot(t_4[30:-30], np.abs(transform_4))
plt.xlabel('t')
plt.ylabel('Amplitude')
plt.grid()
plt.show()

