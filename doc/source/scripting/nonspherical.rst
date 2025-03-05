.. _nonspherical:  

Несферические частицы  
-----------------------  

Формулировка T-матрицы позволяет эффективно выполнять вычисления для одиночных аксиально-симметричных объектов.  
В настоящее время поддерживается только сфероидальная форма (вращающийся эллипсоид) с использованием внешней библиотеки ScatterPy <https://github.com/TCvanLeth/ScatterPy>.  
Подробнее об этом можно узнать в [Mishchenko1998]_.  

Пример  
^^^^^^^  

Расчет экстинкции для уплощенного сфероида с соотношением осей :math:`\alpha=a/c=1.5`.  
Размер задается диаметром :math:`a_{eff}` эквивалентной по объему сферы.  

Размер сфероида можно вычислить следующим образом:  

.. math::  

   a = a_{eff} \alpha^{1/3} \\  
   c = a / \alpha  

.. literalinclude:: spheroid_contrib.py  
   :lines: 2-12  

.. image:: spheroid_shape.png  

.. literalinclude:: spheroid_contrib.py  
   :lines: 14-19  

.. image:: spheroid_ext.png  

Классы  
^^^^^^^  

.. automodule:: mstm_studio.contrib_spheroid  
    :members: SpheroidSP  

.. [Mishchenko1998] M. Mishchenko, L. Travis, "Capabilities and limitations of a current FORTRAN implementation of the T-matrix method for randomly oriented, rotationally symmetric scatterers " (1998) 309  
