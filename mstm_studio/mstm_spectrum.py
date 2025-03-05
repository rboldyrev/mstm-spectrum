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
Основан на сильно переписанном коде MSTM-GUI
<URL:https://github.com/dmayerich/mstm-gui>
<https://git.stim.ee.uh.edu/optics/mstm-gui.git>
автор Dr. David Mayerich

Оптимизирован для спектральных расчетов (для многих длин волн)
с целью использования для подгонки к эксперименту
"""
from __future__ import print_function
from __future__ import division
import numpy as np
from numpy.random import lognormal
from scipy import interpolate
import subprocess
import os   # для удаления файлов после расчета
import sys  # для проверки, работает ли на Linux или Windows
import datetime
import time
import tempfile  # для запуска mstm во временной директории
try:
    import matplotlib.pyplot as plt
except ImportError:
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


class Profiler(object):
    '''
    Этот класс для бенчмаркинга взят из
    http://onesteptospace.blogspot.pt/2013/01/python.html
    Использование:
    >>> with Profiler() as p:
    >>>     // ваш код для профилирования здесь
    '''
    def __enter__(self):
        self._startTime = time.time()

    def __exit__(self, type, value, traceback):
        print('Elapsed time: {:.3f} sec'.format(time.time() - self._startTime))


class SpheresOverlapError(Exception):
    pass


class SPR(object):
    '''
    Класс для расчета поверхностного плазмонного резонанса (SPR),
    запуска внешнего кода MSTM.
    Исполняемый файл MSTM должен быть установлен в переменной окружения MSTM_BIN.
    По умолчанию это ~/bin/mstm.x
    '''

    environment_material = 'Air'

    paramDict = {
      'number_spheres': 0,
      'sphere_position_file': '',          # радиус, X,Y,Z [nm], n ,k
      'length_scale_factor': 1.0,          # 2π/λ[nm]
      'real_ref_index_scale_factor': 1.0,  # множитель для сфер
      'imag_ref_index_scale_factor': 1.0,
      'real_chiral_factor': 0.0,        # хиральные пассивные сферы
      'imag_chiral_factor': 0.0,
      'medium_real_ref_index': 1.0,     # показатель преломления среды
      'medium_imag_ref_index': 0.0,
      'medium_real_chiral_factor': 0.0,
      'medium_imag_chiral_factor': 0.0,
      'target_euler_angles_deg': [0.0, 0.0, 0.0],  # игнорируется для расчетов со случайной ориентацией

      'mie_epsilon': 1.0E-12,           # Критерий сходимости для определения количества порядков
                                        # в разложениях Ми. Отрицательное значение - количество порядков.
      'translation_epsilon': 1.0E-8,    # Критерий сходимости для оценки максимального порядка кластерной T-матрицы
      'solution_epsilon': 1.0E-8,       # Точность решения системы линейных уравнений
      't_matrix_convergence_epsilon': 1.0E-6,
      'plane_wave_epsilon': 1E-3,       # Точность разложения падающего поля (как для плоских, так и для гауссовых волн)
      'iterations_per_correction': 20,  # игнорируется для больших 'near_field_translation_distance'
      'max_number_iterations': 2000,    # с учетом всех итераций
      'near_field_translation_distance': 1.0E6,  # может быть большим действительным, малым действительным или отрицательным. НАСТРОЙКА ДЛЯ ПРОИЗВОДИТЕЛЬНОСТИ
      'store_translation_matrix': 0,
      'fixed_or_random_orientation': 1,  # 0 - фиксированная, 1 - случайная
      'gaussian_beam_constant': 0,       # CB = 1/(k ω0). CB = 0 - плоская волна
      'gaussian_beam_focal_point': [0.0, 0.0, 0.0],  # не влияет на результаты для плоской волны и случайных ориентаций
      'run_print_file': '',              # если пусто, будет использован stdout
      'write_sphere_data': 0,            # 1 - подробно, 0 - кратко

      'output_file': 'test.dat',         # должно изменяться для каждого запуска

      'incident_or_target_frame': 0,     # используется для вывода матрицы рассеяния
      'min_scattering_angle_deg': 0.0,
      'max_scattering_angle_deg': 180.0,
      'min_scattering_plane_angle_deg': 0.0,   # выбирает плоскость для фиксированной ориентации
      'max_scattering_plane_angle_deg': 0.0,   # выбирает плоскость для фиксированной ориентации
      'delta_scattering_angle_deg': 1.0,
      'calculate_near_field': 0,       # без расчетов ближнего поля
      'calculate_t_matrix': 1,         # 1 - новый расчет, 0 - использовать старый, 2 - продолжить расчет
      't_matrix_file': 'tmatrix-temp.dat',
      'sm_number_processors': 10,      # фактическое количество процессоров
                                       # минимум до этого значения и предоставляется mpi
    }

    local_keys = ['output_file', 'length_scale_factor',
                  'medium_real_ref_index', 'medium_imag_ref_index',
                  't_matrix_file']

    def __init__(self, wavelengths):
        '''
        Параметр:
            wavelengths: numpy array
                Длины волн в нм
        '''
        self.wavelengths = wavelengths
        self.command = os.environ.get('MSTM_BIN', '~/bin/mstm.x')

    def set_spheres(self, spheres):
        self.spheres = spheres
        # считаем сферы с положительным радиусом:
        self.paramDict['number_spheres'] = np.sum(self.spheres.a > 0)

    def simulate(self, outfn=None):
        '''
        Начать симуляцию.

        Входные параметры читаются из словаря объекта `paramDict`.
        Программа подготовит входной файл `scriptParams.inp` во временной папке,
        который будет удален после расчета.

        После расчета результат зависит от настройки поляризации.
        Для поляризованного света поля объекта будут заполнены:

            extinction_par, extinction_ort,
            absorbtion_par, absorbtion_ort,
            scattering_par, scattering_ort.

        В то время как для усредненной по ориентациям расчетов просто:

            extinction, absorbtion и scattering.
        '''
        if self.paramDict['number_spheres'] == 0:  # нет сфер
            return self.wavelengths, np.zeros_like(self.wavelengths)
        if self.spheres.check_overlap():
            raise SpheresOverlapError('Сферы пересекаются!')
        if isinstance(self.environment_material, Material):
            material = self.environment_material
        else:
            print(self.environment_material)
            material = Material(self.environment_material)
        with tempfile.TemporaryDirectory() as tmpdir:
            print('Using temporary directory: %s' % tmpdir)
            outFID = open(os.path.join(tmpdir, 'scriptParams.inp'), 'w')
            outFID.write('begin_comment\n')
            outFID.write('**********************************\n')
            outFID.write('  MSTM input for SPR calculation\n')
            outFID.write('  Generated by python script\n')
            outFID.write('  %s\n' %
                         datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
            outFID.write('**********************************\n')
            outFID.write('end_comment\n')
            for key in self.paramDict.keys():
                if key not in self.local_keys:
                    outFID.write(key + '\n')
                    if isinstance(self.paramDict[key], str):
                        svalue = self.paramDict[key]
                    else:
                        if isinstance(self.paramDict[key], list):
                            svalue = '  '.join(map(str, self.paramDict[key]))
                        else:
                            svalue = str(self.paramDict[key])
                        # заменяем символ экспоненты
                        svalue = svalue.replace('e', 'd', 1)
                    outFID.write('%s \n' % svalue)

            for l in self.wavelengths:
                outFID.write('begin_comment\n')
                outFID.write('**********************************\n')
                outFID.write('  Wavelength  %.3f \n' % l)
                outFID.write('**********************************\n')
                outFID.write('end_comment\n')
                outFID.write('output_file\n')
                outFID.write('mstm_l%.0f.out\n' % (l * 1000))
                outFID.write('length_scale_factor\n')
                outFID.write('  %.6f\n' % (2.0 * 3.14159 / l))
                outFID.write('medium_real_ref_index\n')
                outFID.write('  %f\n' % material.get_n(l))
                outFID.write('medium_imag_ref_index\n')
                outFID.write('  %f\n' % material.get_k(l))

                outFID.write('sphere_sizes_and_positions\n')

                for i in xrange(len(self.spheres)):
                    a = self.spheres.a[i]
                    if a > 0:  # учитываем только положительные радиусы
                        x = self.spheres.x[i]
                        y = self.spheres.y[i]
                        z = self.spheres.z[i]
                        self.spheres.materials[i].D = 2 * a
                        n = self.spheres.materials[i].get_n(l)
                        k = self.spheres.materials[i].get_k(l)
                        outFID.write('  %.4f  %.4f  %.4f  %.4f  %.3f  %.3f \n' %
                                     (a, x, y, z, n, k))
                outFID.write('new_run\n')

            outFID.write('end_of_options\n')
            outFID.close()

            # запуск бинарного файла
            if sys.platform == 'win32':
                si = subprocess.STARTUPINFO()
                si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                subprocess.call('%s scriptParams.inp > NUL' % self.command,
                                shell=True, startupinfo=si, cwd=tmpdir)
            else:
                subprocess.call('%s scriptParams.inp > /dev/null' % self.command,
                                shell=True, cwd=tmpdir)

            # разбор результатов симуляции
            if self.paramDict['fixed_or_random_orientation'] == 0:  # фиксированная ориентация
                self.extinction_par = []  # параллельная поляризация (\hat \alpha)
                self.absorbtion_par = []
                self.scattering_par = []
                self.extinction_ort = []  # перпендикулярная поляризация (\hat \beta)
                self.absorbtion_ort = []
                self.scattering_ort = []
                for l in self.wavelengths:
                    inFID = open(os.path.join(tmpdir,
                                              'mstm_l%.0f.out' % (l * 1000)),
                                 'r')
                    while True:
                        line = inFID.readline()
                        if 'scattering matrix elements' in line:
                            break
                        elif 'parallel total ext, abs, scat efficiencies' in line:
                            values = map(float,
                                         inFID.readline().strip().split())
                            values = list(values)
                            self.extinction_par.append(float(values[0]))
                            self.absorbtion_par.append(float(values[1]))
                            self.scattering_par.append(float(values[2]))
                        elif 'perpendicular total ext' in line:
                            values = map(float,
                                         inFID.readline().strip().split())
                            values = list(values)
                            self.extinction_ort.append(float(values[0]))
                            self.absorbtion_ort.append(float(values[1]))
                            self.scattering_ort.append(float(values[2]))
                    inFID.close()
                    os.remove(os.path.join(tmpdir,
                                           'mstm_l%.0f.out' % (l * 1000)))
                self.extinction_par = np.array(self.extinction_par)
                self.absorbtion_par = np.array(self.absorbtion_par)
                self.scattering_par = np.array(self.scattering_par)
                self.extinction_ort = np.array(self.extinction_ort)
                self.absorbtion_ort = np.array(self.absorbtion_ort)
                self.scattering_ort = np.array(self.scattering_ort)
                return (self.wavelengths,
                        (self.extinction_par + self.extinction_ort))
            else:    # случайная ориентация
                self.extinction = []
                self.absorbtion = []
                self.scattering = []
                for lam in self.wavelengths:
                    fnl = os.path.join(tmpdir, 'mstm_l%.0f.out' % (lam * 1000))
                    with open(fnl, 'r') as fout:
                        while True:
                            line = fout.readline()
                            if 'scattering matrix elements' in line:
                                break
                            elif 'total ext, abs, scat efficiencies' in line:
                                values = map(float,
                                             fout.readline().strip().split())
                                values = list(values)  # python3 is evil
                                self.extinction.append(float(values[0]))
                                self.absorbtion.append(float(values[1]))
                                self.scattering.append(float(values[2]))
                    os.remove(fnl)
                self.extinction = np.array(self.extinction)
                self.absorbtion = np.array(self.absorbtion)
                self.scattering = np.array(self.scattering)
            if outfn is not None:
                self.write(outfn)
        return self.wavelengths, self.extinction


    def plot(self):
        '''
        Построить результаты с использованием matplotlib.pyplot
        '''
        if self.paramDict["fixed_or_random_orientation"] == 1:  # случайная ориентация
            plt.plot(self.wavelengths, self.extinction, 'r-', label='extinction')
        else:
            plt.plot(self.wavelengths, self.extinction_par, 'r-', label='extinction par.')
            plt.plot(self.wavelengths, self.extinction_ort, 'b-', label='extinction ort.')
        plt.legend()
        plt.show()
        return plt

    def write(self, filename):
        '''
        Сохранить результаты в файл
        '''
        if self.paramDict["fixed_or_random_orientation"] == 1:  # случайная ориентация
            fout = open(filename, 'w')
            fout.write('#Wavel.\tExtinct.\n')
            for i in range(len(self.wavelengths)):
                fout.write('%.4f\t%.8f\r\n' % (self.wavelengths[i],
                                               self.extinction[i]))
            fout.close()
        else:   # фиксированная ориентация
            fout = open(filename, 'w')
            fout.write('#Wavel.\tExt_par\tExt_ort\n')
            for i in range(len(self.wavelengths)):
                fout.write('%.4f\t%.8f\t%.8f\r\n' % (self.wavelengths[i],
                           self.extinction_par[i], self.extinction_ort[i]))
            fout.close()

    def set_incident_field(self, fixed=False, azimuth_angle=0.0,
                           polar_angle=0.0, polarization_angle=0.0):
        '''
            Установить ориентацию и поляризацию падающей волны

            Параметры:

                fixed: bool
                    True  - фиксированная ориентация и поляризованный свет
                    False - усреднение по всем ориентациям и поляризациям

                azimuth_angle, polar_angle: float (градусы)

                polarization_angle: float (градусы)
                    !имеет смысл только для расчета ближнего поля!
                    угол поляризации относительно плоскости `k-z`.
                    0 - X-поляризация, 90 - Y-поляризация (если `azimuth` и
                    `polar` углы равны нулю).
        '''
        if not fixed:
            self.paramDict['fixed_or_random_orientation'] = 1  # случайная ориентация
        else:
            self.paramDict['fixed_or_random_orientation'] = 0  # фиксированная ориентация
            self.paramDict['incident_azimuth_angle_deg'] = azimuth_angle
            self.paramDict['incident_polar_angle_deg'] = polar_angle
            self.paramDict['polarization_angle_deg'] = polarization_angle


class Material(object):
    r"""
    Класс материала.

    Используйте методы `get_n()` и `get_k()` для получения значений показателя
    преломления на произвольной длине волны (в нм).
    """
    def __init__(self, file_name, wls=None, nk=None, eps=None):
        r"""
        Параметры:

        file_name:
            1. комплексное значение, записанное в формате numpy или как строка;
            2. одна из предопределенных строк (air, water, glass);
            3. имя файла с оптическими константами.

            Заголовок файла должен указывать столбцы `lambda`, `n` и `k`
            Если указаны либо `nk= n + 1j*k`, либо `eps = re + 1j*im` массивы,
            то данные из одного из них будут использованы,
            а содержимое файла будет проигнорировано.

        wls: массив float
            массив длин волн (в нм), используемых для интерполяции данных.
            Если None, то будет использован ``np.linspace(300, 800, 500)``.

        """
        if isinstance(file_name, str):
            self.__name__ = 'Mat_%s' % os.path.basename(file_name)
        else:
            self.__name__ = 'Mat_%.3f' % file_name

        if wls is None:
            wl_min = 200   # 149.9
            wl_max = 1200  # 950.1
            wls = np.array([wl_min, wl_max])
        k = np.array([0.0, 0.0])
        if nk is not None:
            n = np.real(nk)
            k = np.imag(nk)
        elif eps is not None:
            mod = np.absolute(eps)
            n = np.sqrt((mod + np.real(eps)) / 2)
            k = np.sqrt((mod - np.real(eps)) / 2)
        else:
            try:
                np.cdouble(file_name)
                is_complex = True
            except ValueError:
                is_complex = False
            if is_complex:
                nk = np.cdouble(file_name)
                n = np.array([np.real(nk), np.real(nk)])
                k = np.array([np.imag(nk), np.imag(nk)])
            else:
                if file_name.lower() == 'air':
                    n = np.array([1.0, 1.0])
                elif file_name.lower() == 'water':
                    n = np.array([1.33, 1.33])
                elif file_name.lower() == 'glass':
                    n = np.array([1.66, 1.66])
                else:
                    optical_constants = np.genfromtxt(file_name, names=True)
                    wls = optical_constants['lambda']
                    if np.max(wls) < 100:  # длины волн в микрометрах
                        wls = wls * 1000   # конвертируем в нм
                    n = optical_constants['n']
                    k = optical_constants['k']
                    if wls[0] > wls[1]:  # от большего к меньшему
                        wls = np.flipud(wls)  # обратный порядок
                        n = np.flipud(n)
                        k = np.flipud(k)
                    n = n[wls > wl_min]
                    k = k[wls > wl_min]
                    wls = wls[wls > wl_min]
                    n = n[wls < wl_max]
                    k = k[wls < wl_max]
                    wls = wls[wls < wl_max]
        wl_step = np.abs(wls[1] - wls[0])
        if (wl_step > 1.1) and (wl_step < 500):
            interp_kind = 'cubic'                # кубическая интерполяция
        else:  # слишком плотная или слишком редкая сетка, нужна линейная интерполяция
            interp_kind = 'linear'
        # print('Interpolation kind : %s'%interp_kind)
        self._get_n_interp = interpolate.interp1d(wls, n, kind=interp_kind)
        self._get_k_interp = interpolate.interp1d(wls, k, kind=interp_kind)

    def get_n(self, wl):
        return self._get_n_interp(wl)

    def get_k(self, wl):
        return self._get_k_interp(wl)

    def __str__(self):
        return self.__name__

    def plot(self, wls=None, fig=None, axs=None):
        r"""
        построить зависимость ``n`` и ``k`` от длины волны

        Параметры:

            wls: массив float
                массив длин волн (в нм). Если None, то
                будет использован ``np.linspace(300, 800, 500)``.

            fig: фигура matplotlib

            axs: оси matplotlib

        Возвращает:

            созданные/заполненные объекты fig и axs
        """
        if wls is None:
            wls = np.linspace(300, 800, 500)
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        axs.plot(wls, self.get_n(wls), label='Real')
        axs.plot(wls, self.get_k(wls), label='Imag')
        axs.set_ylabel('Показатель преломления')
        axs.set_xlabel('Длина волны, нм')
        axs.legend()
        if flag:
            plt.show()
        return fig, axs


# class MaterialManager():
    # """
    # Кэш для материалов, чтобы уменьшить количество операций ввода-вывода
    # """
    # def __init__(self, wavelengths):
        # self.materials = {}


class Spheres(object):
    """
    Абстрактная коллекция сфер

    Поля объекта:
        N: int
            количество сфер
        x, y, z: массивы numpy
            координаты центров сфер
        a: список или массив
            радиусы сфер
        materials: массив numpy
            объекты Material или строки
    """
    def __init__(self):
        """
        Создает пустую коллекцию сфер. Используйте дочерние классы для непустых!
        """
        self.N = 0
        self.x = []
        self.y = []
        self.z = []
        self.a = []  # радиус
        self.materials = []

    def __len__(self):
        return self.N

    def check_overlap(self, eps=0.001):
        """
        Проверить, пересекаются ли сферы
        """
        result = False
        n = len(self.x)
        for i in xrange(n):
            for j in xrange(i + 1, n):
                dx = abs(self.x[j] - self.x[i])
                dy = abs(self.y[j] - self.y[i])
                dz = abs(self.z[j] - self.z[i])
                Ri = self.a[i]
                Rj = self.a[j]
                dist = np.sqrt(dx * dx + dy * dy + dz * dz)
                if dist < Ri + Rj + eps:
                    # расстояние между сферами меньше суммы их радиусов
                    # но все еще могут быть вложенные сферы, проверяем это
                    if Ri > Rj:
                        result = Ri < dist + Rj + eps
                    else:  # Rj < Ri
                        result = Rj < dist + Ri + eps
                if result:  # избегаем лишних шагов
                    return True
        return result

    def append(self, sphere):
        """
        Добавить данные из объекта SingleSphere

        Параметр:

            sphere: SingleSphere
        """
        self.a = np.append(self.a, sphere.a[0])
        self.x = np.append(self.x, sphere.x[0])
        self.y = np.append(self.y, sphere.y[0])
        self.z = np.append(self.z, sphere.z[0])
        self.materials.append(sphere.materials[0])
        self.N += 1

    def delete(self, i):
        """
        Удалить элемент с индексом `i`
        """
        self.a = np.delete(self.a, i)
        self.x = np.delete(self.x, i)
        self.y = np.delete(self.y, i)
        self.z = np.delete(self.z, i)
        self.materials.pop(i)
        self.N -= 1

    def extend(self, spheres):
        """
        Добавить все элементы из объекта `spheres`
        """
        for i in xrange(len(spheres)):
            self.append(SingleSphere(spheres.x[i], spheres.y[i],
                        spheres.z[i], spheres.a[i], spheres.materials[i]))

    def get_center(self, method=''):
        """
        рассчитать центр масс в предположении равномерной плотности

        Параметр:

            method: строка {''|'mass'}
                Если method == 'mass', то рассчитывается центр масс
                (строго говоря, объемов).
                В противном случае все сферы усредняются равномерно.
        """
        weights = np.ones(self.N)
        if method.lower() == 'mass':
            weights = self.a**3
        Xc = np.sum(np.dot(self.x, weights)) / np.sum(weights)
        Yc = np.sum(np.dot(self.y, weights)) / np.sum(weights)
        Zc = np.sum(np.dot(self.z, weights)) / np.sum(weights)
        return np.array((Xc, Yc, Zc))

    def load(self, filename, mat_filename='etaGold.txt', units='nm'):
        """
            Читает координаты и радиусы сфер из файла.

            Параметры:

                filename: строка
                    файл для чтения

                mat_filename: строка
                    все сферы будут иметь этот материал (хранение
                    материала сферы еще не реализовано)

                units: строка {'mum'|'nm'}
                    единицы измерения расстояния.
                    Если 'mum', то координаты будут масштабированы (x1000)
        """
        x = []
        y = []
        z = []
        a = []
        try:
            f = open(filename, 'r')
            text = f.readlines()
            for line in text:
                if line[0] != '#':  # пропустить комментарии и заголовок
                    words = [w.strip() for w in line.replace(',', '.').split()]
                    data = [float(w) for w in words]
                    a.append(data[0])
                    x.append(data[1])
                    y.append(data[2])
                    z.append(data[3])
            f.close()
        except Exception as err:
            print('Load failed \n %s' % err)
        self.N = len(a)
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.a = np.array(a)
        if units == 'mum':
            self.x = self.x * 1000.0
            self.y = self.y * 1000.0
            self.z = self.z * 1000.0
            self.a = self.a * 1000.0
        self._set_material(mat_filename)

    def save(self, filename):
        """
        Сохраняет координаты и радиусы сфер в файл.

        Параметр:

            filename: строка
        """
        try:
            f = open(filename, 'w')
            f.write('#radius\tx\ty\tz\tn\tk\r\n')
            for i in xrange(self.N):
                wl = 555
                a = self.a[i]
                x = self.x[i]
                y = self.y[i]
                z = self.z[i]
                n = self.materials[i].get_n(wl)
                k = self.materials[i].get_k(wl)
                f.write('%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\r\n' %
                        (a, x, y, z, n, k))
        except Exception as err:
            print('Save failed \n %s' % err)
        finally:
            f.close()


class SingleSphere(Spheres):
    """
    Коллекция сфер с одной сферой
    """
    def __init__(self, x, y, z, a, mat_filename='etaGold.txt'):
        """
        Параметры:

            x, y, z: float
                координаты центров сфер

            a: float
                радиусы сфер

            mat_filename: строка, float, комплексное значение или объект Material
                спецификация материала
        """
        self.N = 1
        self.x = np.array([x])
        self.y = np.array([y])
        self.z = np.array([z])
        self.a = np.array([a])
        if isinstance(mat_filename, Material):
            self.materials = [mat_filename]
        else:
            self.materials = [Material(mat_filename)]


class LogNormalSpheres(Spheres):
    """
    Набор сфер, расположенных на регулярной сетке
    с размерами, распределенными по логнормальному закону.
    В случае пересечения сфер размеры
    должны(?) быть перегенерированы.
    """
    def __init__(self, N, mu, sigma, d, mat_filename='etaGold.txt'):
        """
        Параметры:

            N: int
                количество сфер
            mu, sigma: float
                параметры логнормального распределения
            d: float
                среднее пустое пространство между центрами сфер
            mat_filename: строка или объект Material
                спецификация материала сфер
        """
        # оцениваем размер коробки:
        a = mu  # средний радиус сферы
        A = (N**(1. / 3) + 1) * (d + 2 * a)
        print('Box size estimated as: %.1f nm' % A)
        # A = A*1.5
        Xc = []
        Yc = []
        Zc = []
        x = -A / 2.0
        while x < A / 2.0:
            y = -A / 2.0
            while y < A / 2.0:
                z = -A / 2.0
                while z < A / 2.0:
                    if (x * x + y * y + z * z < A * A / 4.0):
                        Xc.append(x)
                        Yc.append(y)
                        Zc.append(z)
                    z = z + (2 * a + d)
                y = y + (2 * a + d)
            x = x + (2 * a + d)
        print('Desired number of particles: %i' % N)
        print('Number of particles in a box: %i' % len(Xc))
        self.N = min([N, len(Xc)])
        print('Resulted number of particles: %i' % self.N)
        self.x = np.array(Xc)
        self.y = np.array(Yc)
        self.z = np.array(Zc)
        random_a = lognormal(np.log(mu), sigma, self.N)  # nm
        random_a = random_a
        self.a = np.array(random_a)
        if isinstance(mat_filename, Material):
            mat = mat_filename
        else:
            mat = Material(mat_filename)
        self.materials = [mat for i in xrange(self.N)]


class ExplicitSpheres (Spheres):
    def __init__(self, N=0, Xc=[], Yc=[], Zc=[], a=[],
                 mat_filename='etaGold.txt'):
        """
        Создать явно определенные сферы

        Параметры:
            N: int
                количество сфер
            Xc, Yc, Zc: списки или массивы numpy
                координаты центров сфер
            a: список или массив numpy
                радиусы сфер
            mat_filename: строка, список строк, Material или список
                Materials спецификация материала сфер

            Примечание: Если указан только первый массив Xc, то все данные
            предполагаются упакованными в нем,
            т.е.: `Xc = [X1, Y1, Z1, a1, ..., XN, YN, ZN, aN]`
        """
        super(ExplicitSpheres, self).__init__()
        self.N = N
        if N == 0:  # специальный случай пустого объекта
            self.x = []
            self.y = []
            self.z = []
            self.a = []
            return
        if N < len(Xc):  # данные упакованы в Xc
            assert(4 * N == len(Xc))
            self.x = np.zeros(N)
            self.y = np.zeros(N)
            self.z = np.zeros(N)
            self.a = np.zeros(N)
            i = 0
            while i < len(Xc):
                self.x[i // 4] = Xc[i + 0]
                self.y[i // 4] = Xc[i + 1]
                self.z[i // 4] = Xc[i + 2]
                self.a[i // 4] = abs(Xc[i + 3])
                i = i + 4
        else:
            self.x = np.array(Xc)
            self.y = np.array(Yc)
            self.z = np.array(Zc)
            self.a = np.abs(np.array(a))

        if isinstance(mat_filename, (Material, str)):
            # одно имя файла материала для всех сфер
            self._set_material(mat_filename)
        elif isinstance(mat_filename, list):
            # список имен файлов материалов для всех сфер
            if len(mat_filename) == 1:
                self._set_material(mat_filename[0])
            else:
                assert(len(mat_filename) == self.N)
                for mat_fn in mat_filename:
                    # TODO: использовать менеджер материалов, чтобы избежать повторного создания
                    # и лишних чтений файлов
                    if isinstance(mat_fn, Material):
                        self.materials.append(mat_fn)
                    else:
                        self.materials.append(Material(mat_fn))
        else:
            raise Exception('Bad material variable: %s' % str(mat_filename))

        # if self.check_overlap():
            # print('Warning: Spheres are overlapping!')

    def _set_material(self, mat_filename):
        if isinstance(mat_filename, Material):
            mat = mat_filename
        else:
            mat = Material(mat_filename)
        self.materials = [mat for i in xrange(self.N)]


if __name__ == '__main__':
    print('Overlap tests')
    spheres = Spheres()
    print('  Test not overlapped... ')
    spheres.x = [-5, 5]
    spheres.y = [0, 0]
    spheres.z = [0, 0]
    spheres.a = [4, 4]
    assert(not spheres.check_overlap())
    print('  Test overlapped... ')
    spheres.a = [5, 5]
    assert(spheres.check_overlap())
    print('  Test nested... ')
    spheres.x = [0, 0]
    spheres.a = [2, 5]
    assert(not spheres.check_overlap())
    spheres.a = [5, 3]
    assert(not spheres.check_overlap())
    # input('Press enter')

    print('Materials test')
    mat = Material(os.path.join('nk', 'etaGold.txt'))
    # mat.plot()
    mat1 = Material(os.path.join('nk', 'etaSilver.txt'))
    mat3 = Material('glass')
    mat5 = Material(1.5)
    mat6 = Material('2.0+0.5j')
    mat7 = Material('mat7', wls=np.linspace(300, 800, 100),
                    nk=np.linspace(-10, 5, 100) + 1j * np.linspace(0, 10, 100))
    mat8 = Material('mat7', wls=np.linspace(300, 800, 100),
                    eps=np.linspace(-10, 5, 100) + 1j * np.linspace(0, 10, 100))
    print('etaGold ', mat.get_n(800))
    print('etaSilver ', mat1.get_n(800))
    print('Glass (constant) ', mat3.get_n(800), mat3.get_k(800))
    print('n=1.5 material ', mat5.get_n(550))
    print('n=2.0+0.5j material ', mat6.get_n(550), mat6.get_k(550))
    print('nk material ', mat7.get_n(550), mat7.get_k(550))
    print('eps material ', mat8.get_n(550), mat8.get_k(550))
    # input('Press enter')
    with Profiler() as p:
        wls = np.linspace(300, 800, 100)
        # create SPR object
        spr = SPR(wls)
        spr.environment_material = 'glass'
        # spr.set_spheres(SingleSphere(0.0, 0.0, 0.0, 25.0, 'etaGold.txt'))
        spheres = ExplicitSpheres(2, [-20, 0, 0, 10, 10, 0, 0, 12],
                                  mat_filename=['nk/etaGold.txt',
                                                'nk/etaSilver.txt'])
        # spheres = ExplicitSpheres(2, [0,0,0,20,0,0,0,21],
        #                           mat_filename='etaGold.txt')
        spr.set_spheres(spheres)
        spr.set_incident_field(fixed=True, azimuth_angle=90, polar_angle=90,
                               polarization_angle=45)
        # spr.set_spheres(LogNormalSpheres(27, 0.020, 0.9, 0.050 ))
        # calculate!
        # spr.command = ''
        spr.simulate()
    spr.write('test.dat')
    spr.plot()
    input('Press enter')