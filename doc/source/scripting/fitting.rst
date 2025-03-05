...

Подгонка
--------

Подгонка экспериментальных спектров - это мощный инструмент для изучения плазмонных наночастиц.
Подгонка с помощью теории Ми обычно используется для определения размеров частиц, но подгонка с помощью MSTM может
анализировать даже агломераты (скопления) наночастиц, где теория Ми неприменима, см. [Avakyan2017]_ для примера.
Другим применением является подгонка ядерно-оболочечных или многослойных частиц.

MSTM-studio использует жёстко заданную целевую (штрафную) функцию, которая минимизируется в процессе подгонки (ChiSq):

.. math::

    \chi^2 = \sum_i \left( y_i^\text{(fit)} - y_i^\text{(dat)} \right)^2,    
    
где индекс `i` обозначает длины волн.


Пример: подгонка по теории Ми
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Пример экспериментального файла, включённого в дистрибутив, - это спектр экстинкции золотых частиц, внедрённых в стекло лазером,
синтезированных и изученных Максимилианом Хайнцем [Avakyan2017]_.

.. literalinclude:: fit_by_Mie.py
   :lines: 3-27

Выходные данные (финальная часть)::

    ChiSq:	0.000219
    Optimal parameters
        ext00:	0.035177	(Varied:True)
        ext01:	-0.000049	(Varied:True)
        ext02:	0.007908	(Varied:True)
        ext03:	4.207724	(Varied:True)
        ext04:	0.284066	(Varied:True)
        scale:	7030.322097	(Varied:True)


.. image:: fit_by_Mie.png

Низкое значение `ChiSq` и визуальная проверка соответствия теоретических и экспериментальных кривых
указывает на *приемлемую* подгонку.
Названия параметров подгонки объясняются в разделе Ограничения (см. :class:`.Parameter`).
В этом примере `ext00` и `ext01` представляют параметры `a` и `b` линейного вклада,
`ext02` - коэффициент масштабирования для вклада Ми, `ext03` и `ext04` соответствуют параметрам `mu` и `sigma`
логнормального распределения (см. :class:`mstm_spectrum.MieLognormSpheres`).
Последний параметр, общий коэффициент `scale`, на 100% коррелирует с `ext02`, что приводит к ложным абсолютным значениям.
При необходимости концентрация частиц может быть оценена по их произведению :math:`scale \times ext02` или путём наложения ограничений на один из них.


Класс Fitter
^^^^^^^^^^^^

.. autoclass:: mstm_studio.fit_spheres_optic.Fitter
    :members:


Ограничения
-----------

Ограничения позволяют ускорить или направить подгонку.
Их настройка требует указания имен переменных, которые описаны в документации к классу Parameter:

.. autoclass:: mstm_studio.fit_spheres_optic.Parameter
    :members:


Пример: подгонка по модели ядро-оболочка
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Подгонка того же эксперимента, что и выше, но с использованием модели частицы ядро-оболочка, просто для демонстрации методики.

.. literalinclude:: fit_by_core-shell.py
   :lines: 3-35

Выходные данные (финальная часть)::

    ChiSq:	0.002354
    Оптимальные параметры
        a00:	1.284882	(Изменяемый:True)
        a01:	1.958142	(Изменяемый:True)
        ext00:	0.186312	(Изменяемый:True)
        ext01:	-0.000247	(Изменяемый:True)
        scale:	-0.063814	(Изменяемый:True)
        x00:	0.000000	(Изменяемый:False)
        x01:	0.000000	(Изменяемый:False)
        y00:	0.000000	(Изменяемый:False)
        y01:	0.000000	(Изменяемый:False)
        z00:	0.000000	(Изменяемый:False)
        z01:	0.000000	(Изменяемый:False)

.. image:: fit_by_core-shell.png

Качество подгонки, демонстрируемое параметром ChiSq, в ~10 раз хуже по сравнению с использованием ансамбля невзаимодействующих золотых частиц.
График также показывает неприемлемое качество подгонки.

Классы ограничений
^^^^^^^^^^^^^^^^^^

.. autoclass:: mstm_studio.fit_spheres_optic.Constraint
    :members:

.. autoclass:: mstm_studio.fit_spheres_optic.FixConstraint
    :members:
    
.. autoclass:: mstm_studio.fit_spheres_optic.EqualityConstraint
    :members:
    
.. autoclass:: mstm_studio.fit_spheres_optic.ConcentricConstraint
    :members:
    
.. autoclass:: mstm_studio.fit_spheres_optic.RatioConstraint
    :members:

.. [Avakyan2017] L. Avakyan, M. Heinz, A. Skidanenko, K. Yablunovskiy, J. Ihlemann, J. Meinertz, C. Patzig, M. Dubiel, L. Bugaev "Insight on agglomerates of gold nanoparticles in glass based on surface plasmon resonance spectrum: Study by multi-spheres T-matrix method" J. Phys.: Condens. Matter (2018) *30*, 045901-045909 <https://doi.org/10.1088/1361-648X/aa9fcc>

