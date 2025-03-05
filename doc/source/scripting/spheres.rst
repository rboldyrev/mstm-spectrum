.. _spheres:  

Настройка сфер  
-------------  

Код MSTM требует явного указания позиций и размеров сфер. В MSTM-studio разработано несколько классов для упрощения этой настройки.  

Пример: агрегат частиц  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^  

Скрипт для создания сфер случайного размера, размещённых на регулярной сетке:  

.. literalinclude:: spheres_aggr.py  
   :lines: 1-8  

Выходные данные::  

    Box size estimated as: 77.0 nm  
    Desired number of particles: 9  
    Number of particles in a box: 8  
    Resulted number of particles: 8  
    spheres are overlapping, regenerating...  
    Box size estimated as: 77.0 nm  
    Desired number of particles: 9  
    Number of particles in a box: 8  
    Resulted number of particles: 8  
    [ 9.307773    8.61185299  9.92867988  8.84140858  9.87175352  8.71090184  
      9.71505038 12.40459688]  

Классы  
^^^^^^^  

.. autoclass:: mstm_studio.mstm_spectrum.Spheres  
    :members:  

.. autoclass:: mstm_studio.mstm_spectrum.SingleSphere  
    :members:  

.. autoclass:: mstm_studio.mstm_spectrum.ExplicitSpheres  
    :members:  

.. autoclass:: mstm_studio.mstm_spectrum.LogNormalSpheres  
    :members:  

Запуск MSTM  
--------  

Формализм T-матрицы, предложенный Уотерманом [Khlebtsov2013]_, является одной из обобщений теории Ми для множества сферических объектов.  
Fortran-код Multi Sphere T-matrix (MSTM) разработан Мищенко и Маковски [Mackowski2011]_.  
Класс `SPR` реализует функциональность, необходимую для расчёта спектра экстинкции в видимом диапазоне.  
Обратите внимание, что Fortran-код обладает более широким функционалом, включая расчёты ближнего поля, угловые зависимости и т. д., которые пока не реализованы.  
Для подробностей обратитесь к сайту MSTM <http://eng.auburn.edu/users/dmckwski/scatcodes/>.  

Пример: частица с ядром и оболочкой  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^  

Коэффициент поглощения (нормированное сечение) частицы с золотым ядром и серебряной оболочкой.  

.. literalinclude:: core-shell_mstm.py  
   :lines: 3-21  

.. image:: core-shell_mstm.png  

Класс  
^^^^^  

.. autoclass:: mstm_studio.mstm_spectrum.SPR  
    :members:  

.. [Khlebtsov2013] N. Khlebtsov, "T-matrix method in plasmonics: An overview" J. Quant. Spectrosc. Radiat. Transfer (2013) *123*, 184-217, Peter C. Waterman and his scientific legacy  

.. [Mackowski2011] D. Mackowski, M. Mishchenko, "A Multiple Sphere T-matrix Fortran Code for Use on Parallel Computer Clusters" J. Quant. Spectrosc. Radiat. Transfer (2011) *112*, 2182–2192  
