# -*- coding: utf-8 -*-
#-----------------------------------------------------#
#                                                     #
# Этот код является частью проекта подбора T-матрицы  #
# Вкладчики:                                          #
#  L. Avakyan <laavakyan@sfedu.ru>                    #
#  A. Skidanenko <ann.skidanenko@ya.ru>               #
#                                                     #
#-----------------------------------------------------#
"""
  Подбор спектра T-матрицы агрегата частиц
  к экспериментальному спектру SPR.
"""
from __future__ import print_function
import os
from mstm_studio.mstm_spectrum import SPR, ExplicitSpheres, SpheresOverlapError
from mstm_studio.contributions import ConstantBackground
import numpy as np
from scipy import interpolate
import scipy.optimize as so
from datetime import datetime
import threading

try:
    import matplotlib.pyplot as plt
except:
    pass

# использование input в Python2 и Python3
try:
   input = raw_input
except NameError:
   pass
# использование xrange в Python2 и Python3
try:
    xrange
except NameError:
    xrange = range

class Parameter(object):
    """
    Класс для объекта параметра, используемого для хранения
    имени параметра, его значения и пределов вариации.

    Соглашения по именованию параметров:

        `scale` - внешний общий множитель

        `ext%i` - дополнительный параметр, например, фон, пики или вклады Ми

        `a%i` - радиус сферы

        `x%i`, `y%i`, `z%i` - координаты центра сферы

    где `%i` - это число (0, 1, 2, ...)
    """
    def __init__(self, name, value=1, min=None, max=None, internal_loop=False):
        """
        Параметры:

            name: строка
                имя параметра, используемое для ограничений и т.д.

            value: float
                начальное значение параметра

            min, max: float
                границы вариации параметра (опционально)

            internal_loop : bool
                если `True`, параметр будет разрешен для вариации во внутреннем
                (быстром) цикле, который не требует пересчета MSTM.
                Примечание: этот флаг будет удален в будущем.

            varied: bool
                если `True` - будет изменяться в процессе подбора
        """
        self.name = name
        self.value = self.ini_value = value
        self.min = min
        self.max = max
        self.internal_loop = internal_loop
        self.varied = True
        # ...

    def __str__(self):
        return '%s : %f' % (self.name, self.value)


class Constraint(object):
    """
    Абстрактный класс ограничения. Все остальные должны наследоваться от него.
    """
    def apply(self, params):
        """
        Изменяет словарь params
        в соответствии с заданным алгоритмом ограничения.

        Примечание: Абстрактный метод!
        """
        pass


class FixConstraint(Constraint):
    def __init__(self, prm, value=None):
        """
        Фиксирует значение параметра с именем `prm` на значении `value`.

        Параметры:

            prm: строка
                имя параметра

            value: float
                если `None`, то будет использовано начальное значение.
        """
        self.prm = prm.lower()
        self.value = value

    def apply(self, params):
        """ Применяет фиксирующее ограничение """
        assert self.prm in params
        if self.value is not None:
            params[self.prm].value = self.value
        params[self.prm].varied = False


class EqualityConstraint(Constraint):
    def __init__(self, prm1, prm2):
        """
        Фиксирует два параметра с именами `prm1` и `prm2` равными друг другу
        """
        self.prm1 = prm1.lower()
        self.prm2 = prm2.lower()

    def apply(self, params):
        """ Применяет ограничение равенства """
        assert self.prm1 in params
        assert self.prm2 in params
        params[self.prm2].value = params[self.prm1].value
        params[self.prm2].varied = False


class ConcentricConstraint(Constraint):
    def __init__(self, i1, i2):
        """
        Две сферы с общими центрами.

        `i1` и `i2` -- индексы сфер
        """
        self.constraints = [EqualityConstraint('x%02i'%i1, 'x%02i'%i2),
                            EqualityConstraint('y%02i'%i1, 'y%02i'%i2),
                            EqualityConstraint('z%02i'%i1, 'z%02i'%i2)]

    def apply(self, params):
        """ Применяет ограничение концентричности """
        for c in self.constraints:
            c.apply(params)


class RatioConstraint(Constraint):
    def __init__(self, prm1, prm2, ratio=1):
        """
        Поддерживает соотношение двух переменных, `prm1`/`prm2` = `ratio`
        """
        self.prm1 = prm1.lower()
        self.prm2 = prm2.lower()
        self.set_ratio(ratio)

    def apply(self, params):
        """ Применяет ограничение соотношения """
        assert self.prm1 in params
        assert self.prm2 in params
        params[self.prm2].value = params[self.prm1].value / self.ratio
        params[self.prm2].varied = False

    def set_ratio(self, ratio):
        """
        Устанавливает соотношение :math:`prm1/prm2 = ratio`.
        """
        assert np.abs(ratio) > 1e-10
        self.ratio = ratio


class Fitter(threading.Thread):
    """
    Класс для выполнения подбора экспериментального спектра экстинкции

    Поле:

        tolerance: float
            критерий остановки, по умолчанию 1e-4
    """

    tolerance = 1e-4  # критерий остановки

    def __init__(self, exp_filename, wl_min=300, wl_max=800, wl_npoints=51,
                 extra_contributions=None, plot_progress=False):
        """
        Параметры:

            exp_filename: строка
                имя файла с экспериментальными данными

            wl_min, wl_max: float
                границы длин волн для подбора (в нм).

            wl_npoints: int
                количество длин волн, на которых будут рассчитываться и сравниваться спектры.

            extra_contributions: список объектов Contribution
                Если `None`, то будет использован ConstantBackground.
                Предполагается, что первый элемент - это фон.
                Если вы не хотите использовать дополнительные вклады, установите пустой список `[]`.

            plot_progress: bool
                Показывать прогресс подбора с использованием matplotlib.
                Должно быть отключено при запуске на параллельном кластере без графического интерфейса.
        """
        super(Fitter, self).__init__()
        self._stop_event = threading.Event()  # для возможности остановки извне

        self.exp_filename = exp_filename
        data = np.loadtxt(self.exp_filename)    # загрузка данных
        data = data[np.argsort(data[:,0]),:]    # сортировка по 0-му столбцу
        if np.max(data[:,0]) < 10:      # если значения действительно низкие
            print('WARNING: Data X column is probably in mum, automatilcally rescaling to nm.')
            data[:,0] = data[:,0] * 1000

        self.wl_min = max(wl_min, data[ 0, 0])
        self.wl_max = min(wl_max, data[-1, 0])
        self.wl_npoints = wl_npoints
        print('Wavelength limits are setted to: %f < wl < %f'% (self.wl_min, self.wl_max))
        self.wls, self.exp = self._rebin(self.wl_min, self.wl_max, self.wl_npoints,
                                         data[:,0], data[:,1])
        self.params = {}             # словарь объектов параметров
        self.spheres = None          # объект сфер
        self.constraints = []        # список объектов ограничений
        self.calc = np.zeros_like(self.wls)  # рассчитанный спектр
        self.chisq = -1              # хи-квадрат (квадрат невязки)
        # установка масштаба по умолчанию
        self.set_scale()
        # добавление дополнительных вкладов
        self.extra_contributions = []
        self.set_extra_contributions(extra_contributions)
        # установка материала матрицы по умолчанию
        self.set_matrix()
        # отображение прогресса, если указано
        self.plot_progress = plot_progress
        if self.plot_progress:
            plt.ion()
            self.fig = plt.figure()
            ax = self.fig.add_subplot(111)
            ax.plot(self.wls, self.exp, 'ro')
            self.line1, = ax.plot(self.wls, self.calc, 'b-')
            self.fig.canvas.draw()
            #~ self.lock = threading.Lock()  # используется для синхронизации с основным потоком, где происходит отрисовка
        # функция обратного вызова, задаваемая извне
        self._cbuser = None

    def _rebin(self, xmin, xmax, N, x, y):
        """
        Скрытый метод, используемый для пересчета данных на равномерную шкалу
        """
        f = interpolate.interp1d(x, y)
        xnew = np.linspace(xmin, xmax, N)
        ynew = f(xnew)
        return xnew, ynew

    def _print_params(self):
        for key in sorted(self.params):
            print('params[ %s ] \t %s ' % (key, self.params[key]))

    def set_matrix(self, material='AIR'):
        """
        Устанавливает показатель преломления материала матрицы

        material : {'AIR'|'WATER'|'GLASS'} или float
            название материала или
            значение показателя преломления.
        """
        self.MATRIX_MATERIAL = material

    def set_scale(self, value=1):
        if 'scale' in self.params:
            self.params['scale'].value = value
            self.params['scale'].ini_value = value
        else:
            self.params['scale'] = Parameter('scale', value=value, internal_loop=True)

    def set_extra_contributions(self, contributions, initial_values=None):
        """
        Добавляет дополнительные вклады и инициализирует соответствующие параметры.

        Параметры:

            contributions: список объектов Contribution

            initial_values: массив float
        """
        # удаление старых параметров
        i_tot = 0
        for contribution in self.extra_contributions:
            if contribution is not None:
                n = contribution.number_of_params
                for _ in range(n):
                    self.params.pop('ext%02i' % i_tot)
                    i_tot += 1

        if contributions is None:
            contributions = [ConstantBackground(self.wls, 'ConstBkg')]
        self.extra_contributions = contributions[:]

        # создание новых объектов параметров
        n_tot = 0
        i_tot = 0
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            for _ in range(n):
                self.params['ext%02i' % i_tot] = Parameter('ext%02i' % i_tot, value=0.1, internal_loop=True)
                i_tot += 1
            n_tot += n
        self.extra_contrib_params_count = n_tot
        if initial_values is not None:
            assert len(initial_values) == n_tot
            for i in range(n_tot):
                if initial_values[i] is not None:
                    self.params['ext%02i' % i].value = initial_values[i]
                    self.params['ext%02i' % i].ini_value = initial_values[i]
        # print(self.params)

    def set_spheres(self, spheres):
        """
        Задает сферы для подбора.

        Параметр:

            spheres: список объектов mstm_spectrum.Sphere
                Если `None`, то MSTM не будет запущен.
        """
        if self.spheres is not None:  # удаление параметров старых сфер
            for i in xrange(self.spheres.N):
                self.params.pop('a%02i' % i)
                self.params.pop('x%02i' % i)
                self.params.pop('y%02i' % i)
                self.params.pop('z%02i' % i)
        if spheres is not None:
            self.spheres = spheres
            for i in xrange(self.spheres.N):
                self.params['a%02i' % i] = Parameter('a%02i' % i, self.spheres.a[i])
                self.params['x%02i' % i] = Parameter('x%02i' % i, self.spheres.x[i])
                self.params['y%02i' % i] = Parameter('y%02i' % i, self.spheres.y[i])
                self.params['z%02i' % i] = Parameter('z%02i' % i, self.spheres.z[i])
        else:
            self.set_spheres(ExplicitSpheres())  # пустой объект сфер

    def _update_spheres(self):
        """
        Устанавливает радиусы и положения сфер в соответствии со значениями из словаря params
        """
        assert self.spheres is not None
        for i in xrange(len(self.spheres)):
            self.spheres.a[i] = self.params['a%02i' % i].value
            self.spheres.x[i] = self.params['x%02i' % i].value
            self.spheres.y[i] = self.params['y%02i' % i].value
            self.spheres.z[i] = self.params['z%02i' % i].value

    def _update_params(self, values, internal=False):
        """
        Устанавливает значения из оптимизированных параметров в params

        internal : bool
            если True, то будут обновлены внутренние переменные (масштаб, фон, ..)
        """
        try:  # если не итерируемый (одно значение в values)
            len(values)
        except:
            print('WARNING: values is not a list')
            values = [values]
        if internal:  # внутренние (быстрые) параметры цикла
            self.params['scale'].value = values[0]
            for i in range(self.extra_contrib_params_count):
                self.params['ext%02i' % i].value = values[i+1]  # 0-й - это масштаб
            print('inner: ', (self.extra_contrib_params_count+1), values)
        else:
            # применение ограничений, -- пока работает только для MSTM
            for c in self.constraints:  # применение до
                c.apply(self.params)
            # обновление параметров
            i_tot = 0
            for i in range(len(self.spheres)):
                for key in ('a%02i'%i,'x%02i'%i,'y%02i'%i,'z%02i'%i):
                    if self.params[key].varied:
                        self.params[key].value = values[i_tot]
                        i_tot += 1
            assert i_tot == len(values)
            for c in self.constraints:  # и применение после
                c.apply(self.params)
            self.report_result(msg='[%s] Scale: %.3f Bkg: %.2f\n' % (str(datetime.now()),
                self.params['scale'].value, self.params['ext00'].value))  # может быть многословным!

    def add_constraint(self, cs):
        """
        Добавляет ограничения на параметры.
        Полезно для случая структур с оболочкой и слоистых структур.

        Параметр:

            cs: объект Contraint или список объектов Contraint
        """
        try:
            _ = iter(cs)
        except TypeError:
            cs = [cs]
        for c in cs:
            self.constraints.append(c)

    def _get_spectrum(self):
        """
        Вычисляет спектр агрегатов с использованием модуля mstm_spectrum.
        """
        if self.stopped():
            raise Exception('Fitting interrupted')

        #~ self._apply_constraints
        if len(self.spheres) > 0:
            spr = SPR(self.wls)
            spr.environment_material = self.MATRIX_MATERIAL

            self._update_spheres()
            spr.set_spheres(self.spheres)
            #~ self.lock.acquire()
            try:
                _, extinction = spr.simulate()
                self.result = np.array(extinction)
            except SpheresOverlapError as e:
                self.chisq = 666  # большое злое значение
                return np.zeros_like(self.wls)
            except Exception as e:
                print(e)  # пусть пользователь решает
                raise e
            #~ finally:
                #~ self.lock.release()
        else:  # пустой список сфер
            self.result = np.zeros_like(self.wls)

        # выполнение быстрого подбора по внутренним переменным (масштаб, фон, ..)
        values_internal = []
        values_internal.append(self.params['scale'].value)
        for i in range(self.extra_contrib_params_count):
            values_internal.append(self.params['ext%02i' % i].value)

        def _target_func_int(values):
            """ целевая функция для внутреннего подбора (быстрый цикл) """
            self._update_params(values, internal=True)

            y_dat = self.exp
            assert self.params['scale'].value == values[0]
            y_fit = values[0] * self.result
            n_tot = 1  # масштаб - это values[0]
            for contribution in self.extra_contributions:
                n = contribution.number_of_params
                y_fit += contribution.calculate(values[n_tot:n_tot+n])
                n_tot += n
            self.chisq = np.sum((y_fit - y_dat)**2)
            #~ self.chisq = np.sum((y_fit - y_dat)**2 * (y_dat/np.max(y_dat)+0.001)) / np.sum((y_dat/np.max(y_dat)+0.001))
            #~ self.chisq = np.sum( (y_fit - y_dat)**2 * y_dat**3 ) * 1E3
            print(self.chisq)
            return self.chisq

        #~ print('/ Internal fit loop /')
        result_int = so.minimize(fun=_target_func_int, x0=values_internal, method='BFGS', tol=self.tolerance,
                                 options={'maxiter':100, 'disp':False})
        values_internal = result_int.x

        self._update_params(values_internal, internal=True)

        self.calc = self.params['scale'].value * self.result
        n_tot = 1  # 0-й - это масштаб
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            self.calc += contribution.calculate(values_internal[n_tot:n_tot+n])
            n_tot += n
        return self.calc

    def get_extra_contributions(self):
        '''
        Возвращает список текущих дополнительных вкладов в спектр
        '''
        result = []
        values_internal = []
        values_internal.append(self.params['scale'].value)
        for i in range(self.extra_contrib_params_count):
            values_internal.append(self.params['ext%02i' % i].value)
        n_tot = 1  # масштаб - это values[0]
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            result.append(contribution.calculate(values_internal[n_tot:n_tot+n]))
            n_tot += n
        return result

    def _target_func(self, values):
        """ основная целевая функция """
        self._update_params(values)

        y_dat = self.exp
        y_fit = self._get_spectrum()
        self.chisq = np.sum((y_fit - y_dat)**2)
        #~ self.chisq = np.sum((y_fit - y_dat)**2 * (y_dat/np.max(y_dat)+0.001)) / np.sum((y_dat/np.max(y_dat)+0.001))
        #~ self.chisq = np.sum( (y_fit - y_dat)**2 * (y_dat + np.max(y_dat*0.001))**3)
        #~ self.chisq = np.sum((y_fit - y_dat)**2 * y_dat**3) * 1E3
        #print(chisq)
        return self.chisq

    def set_callback(self, func):
        """
        Устанавливает функцию обратного вызова, которая будет вызываться на
        каждом шаге внешнего цикла оптимизации.

        Параметр:

            func: функция(values)
                где values -- список значений, передаваемых из оптимизационной процедуры
        """
        self._cbuser = func

    def _cbplot(self, values):
        """
        функция обратного вызова
        """
        #~ self.lock.acquire()  # будет ждать здесь
        #~ try:
        #print('Scale: %0.3f Bkg: %0.3f ChiSq: %.8f'% (self.params['scale'].value,
        #      self.params['bkg0'].value, self.chisq) )
        if self._cbuser is not None:  # вызов пользовательской функции
            self._cbuser(self, values)
        if self.plot_progress:
            self.line1.set_ydata(self.calc)
            self.fig.canvas.draw()
            #~ plt.pause(0.05)  # это приводит к захвату фокуса окном графика
            self.fig.canvas.start_event_loop(0.05)
            #from:
            #https://stackoverflow.com/questions/45729092/make-interactive-matplotlib-window-not-pop-to-front-on-each-update-windows-7
            #~ backend = plt.rcParams['backend']
            #~ import matplotlib
            #~ if backend in matplotlib.rcsetup.interactive_bk:
                #~ figManager = matplotlib._pylab_helpers.Gcf.get_active()
                #~ if figManager is not None:
                    #~ canvas = figManager.canvas
                    #~ if canvas.figure.stale:
                        #~ canvas.draw()
                    #~ canvas.start_event_loop(0.05)
        #~ finally:
            #~ self.lock.release()
        #~ input('pe')

    def _apply_constraints(self):
        for c in self.constraints:
            c.apply(self.params)

    def run(self, maxsteps=400):
        """
        Запускает подбор.

        Параметры:
            maxsteps: int
                ограничивает количество выполняемых шагов
        """
        self._apply_constraints()
        # упаковка параметров в values
        values = []
        for i in range(len(self.spheres)):
            for key in ('a%02i'%i,'x%02i'%i,'y%02i'%i,'z%02i'%i):
                if self.params[key].varied:
                    values.append(self.params[key].value)
        # запуск оптимизатора
        result = so.minimize(fun=self._target_func, x0=values, method='Powell', tol=self.tolerance,
                             options={'maxiter':maxsteps, 'disp':True}, callback=self._cbplot)
        self._update_params(result.x)

    def stop(self):
        # https://stackoverflow.com/questions/323972/is-there-any-way-to-kill-a-thread-in-python
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def report_freedom(self):
        """
        Возвращает строку с кратким описанием перед подбором
        """
        self._apply_constraints()
        N = len(self.spheres)
        s = 'Количество сфер:\t%i\n' % N
        n_tot = 0
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            s += 'Дополнительный вклад %s с %i параметрами\n' % (contribution.name, n)
            n_tot += n
        s += 'Общее количество дополнительных параметров:\t%i\n' % n_tot
        n_fast = n_slow = n_fix = 0
        for key in self.params:
            if self.params[key].varied:
                if self.params[key].internal_loop:
                    n_fast += 1
                else:
                    n_slow += 1
            else: # не изменяется
                n_fix += 1
        assert n_fast + n_slow + n_fix == 1 + 4*N + n_tot
        s += 'Степени свободы\n'
        s += '\tбыстрый цикл:\t%i\n' % n_fast
        s += '\tмедленный цикл:\t%i\n' % n_slow
        print(s)
        return s

    def report_result(self, msg=None):
        """
        Возвращает строку с кратким описанием результатов подбора
        """
        s = 'Хи-квадрат:\t%f\n' % self.chisq
        if msg is None:
            s += 'Оптимальные параметры'
        else:
            s += msg
        for key in sorted(self.params):
            s += '\n\t%s:\t%f\t(Изменяемый:%s)' % (key, self.params[key].value, str(self.params[key].varied))
        print(s)
        return s

if __name__ == '__main__':
    fitter = Fitter('../example/experiment.dat')
    # тест подбора Ми
    from mstm_studio.contributions import LinearBackground, MieSingleSphere, MieLognormSpheresCached
    from mstm_studio.alloy_AuAg import AlloyAuAg
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg'),
                                    MieLognormSpheresCached(fitter.wls, 'LN Mie')],
                                    [0.02, -0.001,
                                    0.1, 1.5, 0.5])
                                    #~ MieSingleSphere(fitter.wls, 'Mie')],
                                    #~ [0.02, -0.001,
                                    #~ 0.1, 10])
    fitter.extra_contributions[1].set_material(AlloyAuAg(1.), 1.66)
    fitter.extra_contributions[1].plot([0.1, 1.5, 0.5])
    fitter.extra_contributions[1].plot_distrib([0.1, 1.5, 0.5])
    #~ fitter.extra_contributions[1].plot([0.1, 10])
    fitter.set_spheres(None)  # нет сфер, нет запусков MSTM
    fitter.report_freedom()
    input('Press enter to run peak fitting')
    fitter.run()
    fitter.report_result()
    contribs = fitter.get_extra_contributions()
    print(contribs)
    input('Press enter to continue')
    # тест подбора пика
    from contributions import LinearBackground, LorentzPeak
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg'),
                                    LorentzPeak(fitter.wls, 'lorentz peak')],
                                    [0.02, -0.001,
                                     100, 550, 50])
    # fitter.extra_contributions[1].plot([100, 550, 50])
    fitter.set_spheres(None)  # нет сфер, нет запусков MSTM
    fitter.report_freedom()
    input('Press enter to run peak fitting')
    fitter.run()
    fitter.report_result()
    input('Press enter to continue')
    # тест подбора MSTM
    fitter.set_matrix('glass')
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg')], [0.02, -0.001])
    #                         N    X      Y      Z    radius    materials
    spheres = ExplicitSpheres(2, [-1,1], [-2,2], [-3,3], [14,20], [AlloyAuAg(1.), AlloyAuAg(0.)])
    fitter.set_spheres(spheres)
    fitter.add_constraint(ConcentricConstraint(0, 1))
    fitter.add_constraint(FixConstraint('x00', 0))
    fitter.add_constraint(FixConstraint('y00', 0))
    fitter.add_constraint(FixConstraint('z00', 0))
    fitter.add_constraint(RatioConstraint('a00', 'a01', spheres.a[0]/spheres.a[1]))
    fitter.report_freedom()
    input('Press enter to run MSTM fitting')

    fitter.run()
    #~ fitter.start()  # метод потока
    #~ fitter.join()   # ожидание завершения
    fitter.report_result()
    #fitter.plot_result()
    #~ y_fit = _get_spectrum( wavelengths, values )
    #~ plt.plot( wavelengths, exp, wavelengths, y_fit )
    #~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
    #~ plt.xlabel('Wavelength, nm')
    #~ plt.ylabel('Exctinction, a.u.')
    #~ plt.show()
    input('Press enter to finish')
    print('It is over.')