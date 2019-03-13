# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#   K. Yablunovskiy <kirill-yablunovskii@mail.ru>       #
#                                                       #
# ----------------------------------------------------- #
"""
Contributions to UV/vis extinction spectra other
then obtained from MSTM.
"""
#TODO: add backgtound contributions and adjust corresponding GUI panels. Remove background from mstm_spectrum.py
#TODO: add Mie contribution
#TDOD: add Kirill's neural network contribution
from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# use input in both python2 and python3
try:
   input = raw_input
except NameError:
   pass
# use xrange in both python2 and python3
try:
    xrange
except NameError:
    xrange = range

try:
    from film_exctinction import gold_film_ex  # for gold film background
except:
    pass


class Contribution(object):
    """
    Abstract class to account for contributions other then
    calculated by MSTM. Here should come all lightweight calculated
    contribtions: constant background, lorentz and guass peaks, Mie, etc.
    Each should be implemented as a child class.
    """
    number_of_params = 0  # Should be another value in child class

    def __init__(self, wavelengths=[], name='ExtraContrib'):
        """
        Parameter object

        name : str
        wavelengths : list or np.array
            wavelengths
        values : list
            adjustable parameters, like constant bakcground, or peak parameters
        """
        self.name = name
        self.wavelengths = np.array(wavelengths)

    def calculate(self, values):
        """
        return np.array of calculated contribution at wavelength.

        This method should be overriden in child classes.
        """
        self._check(values)
        return np.zeros(len(self.wavelengths))

    def _check(self, values):
        if len(values) < self.number_of_params:
            raise Exception('Too few values! '+str(values))

    def plot(self, values, fig=None, axs=None):
        """
        plot contribution
        values : list of parameters
        fig, axs : matplotlib objects
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.wavelengths
        y = self.calculate(values)
        # if isinstance(y, float):  # if not an array
            # y = y * np.ones(len(x))
        axs.plot(x, y, 'g--', label=self.name)
        axs.set_ylabel('Intensity')
        axs.set_xlabel('Wavelength, nm')
        axs.legend()
        if flag:
            plt.show()


class ConstantBackground(Contribution):
    """
    Simple background contribution ruled by single parameter.
    """
    number_of_params = 1

    def calculate(self, values):
        self._check(values)
        return values[0] * np.ones(len(self.wavelengths))


class LinearBackground(Contribution):
    """
    Two-parameter background $ a * {\lambda} + b $.
    """
    number_of_params = 2

    def calculate(self, values):
        self._check(values)
        return values[0] + values[1] * self.wavelengths


class LorentzPeak(Contribution):
    """
    Lorentz function
    """
    number_of_params = 3

    def calculate(self, values):
        self._check(values)
        return values[0] / ((self.wavelengths-values[1])**2 + (values[2])**2)


class LorentzBackground(Contribution):
    """
    Lorentz peak in background. Peak center is fixed.
    """
    number_of_params = 3

    def calculate(self, values, center=250):
        self._check(values)
        return values[0] + values[1] / ((self.wavelengths-center)**2 + (values[2])**2)


class FilmBackground(Contribution):
    """
        Background interpolated from experimental spectra of gold foil
    """
    number_of_params = 3

    def calculate(self, values):
        self._check(values)
        return values[0] + values[1] * gold_film_ex(values[2], self.wavelengths)


if __name__=='__main__':
    # tests come here
    cb = ConstantBackground(name='const', wavelengths=[300,400,500,600,700,800])
    print(cb.calculate([3]))
    cb.plot([3])
    print('See you!')
