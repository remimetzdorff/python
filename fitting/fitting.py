# add self.params_key for each function

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas
from inspect import signature
from scipy.fftpack import fft as fft

def make_fft(x,y):
    N = len(x)
    T = x[1]-x[0]
    yf = fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    return xf, 2.0/N * np.abs(yf[0:N//2])

def find_exponent(value):
    if value == 0:
        exponent = 0
    else:
        exponent = np.floor(np.log10(abs(value)))
    return exponent


def display_readable(param, uparam):
    exponent   = find_exponent(param)
    uexponent  = find_exponent(uparam)
    poweroften = int(exponent - uexponent)
    if np.abs(param) > np.abs(uparam):
        new_uparam = np.round(uparam * 10. ** (-uexponent)) * 10. ** -poweroften
        new_param = np.round(param * 10. ** (-uexponent)) * 10. ** -poweroften
    else:
        poweroften = np.abs(poweroften)
        new_uparam = np.round(uparam * 10. ** (-exponent)) * 10. ** -poweroften
        new_param = np.round(param * 10. ** (-exponent)) * 10. ** -poweroften
    if not new_param < 10:
        # Correct value to avoid 10 * 1e-1
        new_param /= 10.
        new_uparam /= 10.
        exponent += 1
        poweroften += 1
    row = (str("(%." + str(poweroften) + "f +/- %." + str(poweroften) + "f) * 1e%.0f") % (new_param, new_uparam, exponent))
    return row


class Fit():

    def __init__(self, func, data=None, x=None, y=None, uy=None, ux=None, verbosemode=True):
        if data is not None:
            self.data = data
        else:
            self.data = pandas.Series(y, index=x)
        self.x = self.data.index.values
        self.y = self.data.values
        self.uy = uy
        self.ux = ux

        self.func_name = func
        self.func      = (getattr(self, self.func_name))

        guess_function = "_guess" + self.func_name
        self._guessfunction = (getattr(Fit, guess_function))

        self.verbosemode = verbosemode
        self.fit_params = None
        self.number_of_fitparams = len(signature(self.func).parameters)-1

    def peak_detection(self):
        center = self.data.idxmax()
        offset = self.data.min()
        scale = self.data.max() - offset
        threshold = offset + scale / 2
        for xx, yy in zip(self.x, self.y):
            if yy > threshold:
                x1 = xx
                break
            else:
                pass
        for xx, yy in zip(self.x[::-1], self.y[::-1]):
            if yy > threshold:
                x2 = xx
                break
            else:
                pass
        width = x2-x1
        return center, width, scale, offset

    def linear(self, x, slope, y0):
        """
        slope   : slope
        y0      : y0 at x=0
        """
        self.params_key = ["slope", "y0"]
        return slope*x + y0
    def _difflinear(self):
        if self.fit_params is not None:
            slope, y0 = self.fit_params
        else:
            slope, y0 = self._guesslinear()
        return slope
    def _guesslinear(self):
        ymax = self.data.max()
        ymin = self.data.min()
        xmax = self.data.idxmax()
        xmin = self.data.idxmin()
        slope = (ymax-ymin)/(xmax-xmin)
        y0 = self.data.mean() - self.x.mean()*slope
        fit_params = [slope, y0]
        return fit_params

    def pendulum(self, x, T0, a, b):
        self.params_key = ["T0", "a", "b"]
        X = x-0
        return T0 * (1 + a * X**2 + b * X**4)
    def _guesspendulum(self):
        guess_params = np.ones(self.number_of_fitparams)
        return guess_params

    def exponential(self, x, scale, alpha):
        """
        scale   : scale
        alpha   : alpha
        offset  : offset (background)
        """
        self.params_key = ["scale", "alpha"]
        return scale*np.exp(alpha*x)
    def _guessexponential(self):
        guess_params = self._guessexponential_offset()[:2]
        return guess_params

    def exponential_offset(self, x, scale, alpha, y0):
        """
        scale   : scale
        alpha   : alpha
        offset  : offset (background)
        """
        self.params_key = ["scale", "alpha", "offset"]
        return scale*np.exp(alpha*x) + y0
    def _guessexponential_offset(self):
        # https://fr.scribd.com/doc/14674814/Regressions-et-equations-integrales p16-17
        for i in range(len(self.data)):
            if i == 0:
                s = [0]
            else:
                val = s[i - 1] + .5 * (self.y[i] + self.y[i - 1]) * (self.x[i] - self.x[i - 1])
                s.append(val)
        s = np.array(s)
        matrix = [[((self.x - self.x[0]) ** 2).sum(), ((self.x - self.x[0]) * s).sum()],
                  [((self.x - self.x[0]) * s).sum(), (s ** 2).sum()]]
        vector = [((self.y - self.y[0]) * (self.x - self.x[0])).sum(), ((self.y - self.y[0]) * s).sum()]
        a, b = np.linalg.inv(matrix).dot(vector)
        alpha = b
        y0 = -a / b
        # other set of matrix vector required to get scale value
        theta = np.exp(alpha*self.x)
        matrix = [[len(self.x), theta.sum()],
                  [theta.sum(), (theta**2).sum()]]
        vector = [self.y.sum(), (self.y*theta).sum()]
        y0, scale = np.linalg.inv(matrix).dot(vector)
        guess_params = [alpha, scale, y0]
        return guess_params

    def sin(self, x, freq, scale, phi, offset):
        """
        freq   : frequency in 1/x unit
        scale  : amplitude (half peak to peak)
        phi    : phase at x=0
        offset : offset (mean value)
        """
        self.params_key = ["freq", "scale", "phi", "offset"]
        return scale*np.sin(2*np.pi*freq*x + phi) + offset
    def _guesssin(self):
        xf, yf = make_fft(self.x, self.y)
        offset = yf[0]/2*np.sign(self.y.mean())
        maxpos = yf[1:].argmax()
        freq = xf[maxpos+1]
        scale = yf[1:].max()
        phases = np.linspace(0, 2 * np.pi,7)
        diff = []
        for phase in phases:
            diff.append(((self.y-self.sin(self.x, freq, scale, phase, offset))**2).sum())
        phase = phases[np.array(diff).argmin()]
        fit_params = [freq, scale, phase, offset]
        return fit_params

    def sinc(self, x, center, freq, scale, offset):
        """
        center : center
        freq   : frequency in 1/x unit
        scale  : amplitude (half peak to peak)
        offset : offset (mean value)
        """
        self.params_key = ["center", "freq", "scale", "offset"]
        X = 2*np.pi*freq*(x-center)
        return offset + scale*np.sinc(X/np.pi)
    def _guesssinc(self):
        center, width, scale, offset = self.peak_detection()
        scale /= 1.2
        offset += 0.2*scale
        freq = 2 / (np.pi * (width))
        guess_params = [center, freq, scale, offset]
        return guess_params

    def sinc_squarred(self, x, center, freq, scale, offset):
        """
        center : center
        freq   : frequency in 1/x unit
        scale  : amplitude (half peak to peak)
        offset : offset (mean value)
        """
        self.params_key = ["center", "freq", "scale", "offset"]
        return offset + scale * (self.sinc(x, center, freq, 1, 0)) ** 2
    def _guesssinc_squarred(self):
        center, width, scale, offset = self.peak_detection()
        freq = 1/(np.pi*(width))
        guess_params = [center, freq, scale, offset]
        return guess_params

    def gaussian(self, x, center, bandwidth, scale, offset):
        self.params_key = ["center", "bandwidth", "scale", "offset"]
        return scale * np.exp(-1 * ((x - center) / abs(bandwidth)) ** 2) + offset
    def _guessgaussian(self):
        center, width, scale, offset = self.peak_detection()
        guess_params = [center, width/np.sqrt(2), scale, offset]
        return guess_params

    def guess_params(self):
        guess_params=self._guessfunction(self)
        return guess_params

    def fit(self, manualguess_params=None, verbosemode=None):
        if manualguess_params is not None:
            guess_params = manualguess_params
        else:
            try:
                autoguess_params = self.guess_params()
                guess_params = autoguess_params
            except:
                guess_params = None
        if self.ux is not None:
            diff_function = "_diff" + self.func_name
            self._difffunction = (getattr(Fit, diff_function))
            uy = np.sqrt(self.uy**2 + (self._difffunction(self)*self.ux)**2)
        else:
            uy = self.uy
        fit_params, pcov = curve_fit(self.func, self.x, self.y, sigma=uy, p0=guess_params, maxfev=50000, absolute_sigma=True)
        self.fit_params = fit_params
        self.fit_uparams = np.sqrt(np.abs(np.diagonal(pcov)))
        if verbosemode is None:
            if self.verbosemode:
                self.report()
        elif verbosemode:
            self.report()
        return fit_params, self.fit_uparams

    def chi2(self):
        """computed after fit"""
        if self.uy is not None:
            uy = self.uy
        else:
            uy = np.ones(len(self.x))
        if self.ux is not None:
            diff_function = "_diff" + self.func_name
            self._difffunction = (getattr(Fit, diff_function))
            uy = np.sqrt(uy ** 2 + (self._difffunction(self)*self.ux) ** 2)
        return np.sum(((self.y - self.func(self.x, *self.fit_params))/uy)**2)

    def chi2r(self):
        return self.chi2()/(len(self.x)-len(self.fit_params))

    def r2(self):
        num   = np.sum((self.y - self.func(self.x, *self.fit_params))**2)
        denom = np.sum((self.y - self.y.mean())**2)
        return 1 - num / denom

    def report(self):
        print("##### Fit results #####")
        print("FIT FUNCTION:", self.func_name)
        print("RAW")
        print("        params_key      :", self.params_key)
        print("        optimised params:", self.fit_params)
        print("        uncertainties   :", self.fit_uparams)
        print("        chi2r           :", self.chi2r())
        print("        r2              :", self.r2())
        print("READABLE")
        for param, uparam, key in zip(self.fit_params, self.fit_uparams, self.params_key):
            new_key = key
            while (len(new_key) < 10):
                new_key += " "
            print(str("        "+new_key+" = "+display_readable(param, uparam)))
        print("#######################")
