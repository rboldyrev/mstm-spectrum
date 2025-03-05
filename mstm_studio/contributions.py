# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  Этот код является частью проекта подгонки T-матрицы  #
#  Вкладчики:                                           #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#   K. Yablunovskiy <kirill-yablunovskii@mail.ru>       #
#                                                       #
# ----------------------------------------------------- #
"""
Вклады в спектры экстинкции UV/vis, отличные от
полученных с помощью MSTM.
"""
from __future__ import print_function
from __future__ import division
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass

# используем input как в python2, так и в python3
try:
   input = raw_input
except NameError:
   pass
# используем xrange как в python2, так и в python3
try:
    xrange
except NameError:
    xrange = range

try:
    from film_exctinction import gold_film_ex  # для фонового вклада золотой пленки
except:
    pass

try:
    from mstm_studio.mie_theory import calculate_mie_spectra
except:
    pass


class Contribution(object):
    """
    Абстрактный класс для включения вкладов, отличных от
    рассчитанных с помощью MSTM. Все легковесные расчетные
    вклады (постоянный фон, пики Лоренца и Гаусса, Ми и т.д.)
    должны наследоваться от него.
    """
    number_of_params = 0  # Должно быть переопределено в дочернем классе

    def __init__(self, wavelengths=[], name='ExtraContrib'):
        """
        Параметры:

            wavelengths: список или массив numpy
                длины волн в нм

            name: строка
                необязательная метка

        """
        self.name = name
        self.set_wavelengths(wavelengths)

    def set_wavelengths(self, wavelengths):
        """
        Изменить длины волн
        """
        self.wavelengths = np.array(wavelengths)

    def calculate(self, values):
        """
        Этот метод должен быть переопределен в дочерних классах.

        Параметры:

            values: список управляющих параметров

        Возвращает:

            массив numpy значений вклада на указанных длинах волн
        """
        self._check(values)
        return np.zeros(len(self.wavelengths))

    def _check(self, values):
        if len(values) < self.number_of_params:
            raise Exception('Слишком мало значений! '+str(values))

    def plot(self, values, fig=None, axs=None):
        """
        Построить график вклада

        Параметры:

            values: список параметров

            fig: фигура matplotlib

            axs: оси matplotlib

        Возвращает:

            созданные/заполненные объекты fig и axs
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.wavelengths
        y = self.calculate(values)
        axs.plot(x, y, 'g--', label=self.name)
        axs.set_ylabel('Интенсивность')
        axs.set_xlabel('Длина волны, нм')
        axs.legend()
        if flag:
            plt.show()
        return fig, axs


class ConstantBackground(Contribution):
    """
    Постоянный фоновый вклад :math:`f(\lambda) = bkg`.
    """
    number_of_params = 1

    def calculate(self, values):
        """
        Параметры:

            values: [bkg]

        Возвращает:

            массив numpy
        """
        self._check(values)
        return values[0] * np.ones(len(self.wavelengths))


class LinearBackground(Contribution):
    """
    Двухпараметрический фон :math:`f(\lambda) = a \cdot \lambda + b`.
    """
    number_of_params = 2

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `mu` и `Gamma`

        Возвращает:

            массив numpy
        """
        self._check(values)
        return values[0] + values[1] * self.wavelengths


class LorentzPeak(Contribution):
    """
    Функция Лоренца

    .. math::

        L(\lambda) = \\frac {scale} {(\lambda-\mu)^2 + \Gamma^2}

    """
    number_of_params = 3

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `mu` и `Gamma`

        Возвращает:

            массив numpy
        """
        self._check(values)
        return values[0] / ((self.wavelengths - values[1])**2 + (values[2])**2)


class GaussPeak(Contribution):
    """
    Функция Гаусса

    .. math::

        G(\lambda) = scale \cdot \exp\left( - \\frac{(\lambda-\mu)^2}{2\sigma^2} \\right)

    """
    number_of_params = 3

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `mu` и `sigma`

        Возвращает:

            массив numpy
        """
        self._check(values)
        return values[0] * np.exp(-(self.wavelengths - values[1])**2 / (2 * values[2]**2))


class LorentzBackground(Contribution):
    """
    Пик Лоренца на фоне. Центр пика фиксирован.

    .. math::

        L(\lambda) = \\frac {scale} {(\lambda-center)^2 + \Gamma^2}

    """
    number_of_params = 2
    center = 250

    def calculate(self, values):
        self._check(values)
        return values[0] / ((self.wavelengths-self.center)**2 + (values[1])**2)


class FilmBackground(Contribution):
    """
    Фон, интерполированный из экспериментальных спектров золотой пленки
    """
    number_of_params = 3

    def calculate(self, values):
        """ TODO """
        self._check(values)
        return values[0] + values[1] * gold_film_ex(values[2], self.wavelengths)


class MieSingleSphere(Contribution):
    """
    Вклад Ми от одной сферы.

    Подробности широко обсуждаются, см., например, [Kreibig_book1995]_
    """
    number_of_params = 2
    material = None  # экземпляр mstm_spectrum.Material
    matrix = 1.0     # показатель преломления среды

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `diameter`

        Возвращает:

            массив эффективности экстинкции Ми для сферы
        """
        self._check(values)
        if self.material is None:
            raise Exception('Для расчета Ми требуются данные о материале. Остановка.')
        D = np.abs(values[1])
        self.material.D = D
        _, _, mie_extinction, _ = calculate_mie_spectra(
            self.wavelengths, D/2.0, self.material, self.matrix)
        return values[0] * mie_extinction

    def set_material(self, material, matrix=1.0):
        """
        Определить материал сферы и окружающей среды

        Параметры:

            material: объект Material
                материал сферы

            matrix: float, строка или объект Material
                материал окружающей среды

        Возвращает:

            True, если свойства были изменены, False - в противном случае.

        """
        changed = False
        try:
            matr = float(matrix)
        except:
            matr = matrix.get_n(550)  # предполагаем, что это экземпляр Material
        if matr != self.matrix:
            self.matrix = matr
            changed = True
        try:
            material.get_n(550)
            material.get_k(550)
        except:
            raise Exception('Некорректный объект материала')
        if material is not self.material:
            self.material = material
            changed = True
        return changed


class MieLognormSpheres(MieSingleSphere):
    """
    Вклад Ми от ансамбля сфер
    с размерами, распределенными по логнормальному закону
    """
    number_of_params = 3
    diameters = np.logspace(0, 3, 301)
    MAX_DIAMETER_TO_PLOT = 100

    def lognorm(self, x, mu, sigma):
        """
        Форма логнормального распределения:

        .. math::

            LN(D) = \\frac {1}{D \sigma \sqrt{2\pi}} \exp\left( - \\frac{(\log(D)-\mu)^2}{2\sigma^2} \\right)
        """
        return (1.0/(x*sigma*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `mu` и `sigma`

        Возвращает:

            эффективность экстинкции Ми для логнормально распределенных сфер
        """
        self._check(values)
        dD = np.ediff1d(self.diameters, to_begin=1e-3)
        distrib = self.lognorm(self.diameters, np.abs(values[1]), np.abs(values[2]))
        result = np.zeros_like(self.wavelengths)
        for diameter, count in zip(self.diameters, distrib*dD):
            self.number_of_params = 2  # иначе получим ошибку при проверке
            mie_ext = super(MieLognormSpheres, self).calculate(values=[1.0, diameter])
            self.number_of_params = 3  # некрасиво, но за все надо платить
            result += count * mie_ext * D**2  # эффективность -> сечение
        av_diameter = np.sum(self.diameters * distrib * dD) / np.sum(distrib * dD)
        return values[0] * result / av_diameter**2

    def plot_distrib(self, values, fig=None, axs=None):
        """
        Построить график распределения размеров

        Параметры:

            values: список управляющих параметров

            fig: фигура matplotlib

            axs: оси matplotlib

        Возвращает:

            созданные/заполненные объекты fig и axs
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.diameters[self.diameters < self.MAX_DIAMETER_TO_PLOT]
        y = self.lognorm(x, np.abs(values[1]), np.abs(values[2]))
        axs.plot(x, y, 'b', label=self.name)
        axs.set_ylabel('Количество')
        axs.set_xlabel('Диаметр, нм')
        axs.legend()
        if flag:
            plt.show()
        return fig, axs

    def set_material(self, material, matrix=1.0):
        print('matrix: %s' % matrix)
        if super().set_material(material, matrix):
            self._M = None  # очистить кэш при изменении материалов


class MieLognormSpheresCached(MieLognormSpheres):
    """
    Вклад Ми от ансамбля сфер
    с размерами, распределенными по логнормальному закону.

    Кэшированная версия - используйте для ускорения подгонки.
    """
    number_of_params = 3
    diameters = np.logspace(0, 3, 301)
    MAX_DIAMETER_TO_PLOT = 100
    _M = None  # кэшированная матрица, None после инициализации

    def calculate(self, values):
        """
        Параметры:

            values: список управляющих параметров `scale`, `mu` и `sigma`

        Возвращает:

            эффективность экстинкции Ми для логнормально распределенных сфер
        """
        self._check(values)
        if self._M is None:  # инициализация кэшированной матрицы
            print('Building cache...')
            self._M = np.zeros(shape=(len(self.wavelengths), len(self.diameters)))
            self.number_of_params = 2  # иначе получим ошибку при проверке
            for i, diameter in enumerate(self.diameters):
                self._M[:, i] = super(MieLognormSpheres, self).calculate(values=[1.0, diameter])  # D или D/2 ?
            self.number_of_params = 3  # некрасиво, но за все надо платить
            print('Building cache... done')

        dD = np.ediff1d(self.diameters, to_begin=1e-3)
        distrib = self.lognorm(self.diameters, np.abs(values[1]), np.abs(values[2]))
        result =  np.dot(self._M, distrib * self.diameters**2 * dD) / np.sum(distrib * dD)
        return values[0] * result


if __name__=='__main__':
    # тесты здесь
    # ~ cb = ConstantBackground(name='const', wavelengths=[300,400,500,600,700,800])
    # ~ print(cb.calculate([3]))
    # ~ cb.plot([3])
    # ~ mie = MieSingleSphere(name='mie', wavelengths=np.linspace(300,800,50))
    # ~ mie = MieLognormSpheres(name='mie', wavelengths=np.linspace(300,800,50))
    # ~ from mstm_studio.contributions import MieLognormSpheresCached
    from mstm_studio.alloy_AuAg import AlloyAuAg
    mie = MieLognormSpheresCached(name='mie', wavelengths=np.linspace(300, 800, 50))
    mie.set_material(AlloyAuAg(x_Au=1), 1.66)
    mie.plot([1,1.5,0.5])  # scale mu sigma
    print('See you!')