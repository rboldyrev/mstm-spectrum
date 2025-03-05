   .. _contribs:


   Простые функции и теория Ми
   -------------------------------

   Пример
   ^^^^^^^

   Пример того, как получить вклад в экстинкцию от сфер с логнормальным распределением. Другие вклады рассчитываются аналогичным образом.

   .. literalinclude:: mie_contrib.py
      :lines: 3-12


   .. image:: mie_contrib.png


   .. literalinclude:: mie_contrib.py
      :lines: 14-16


   .. image:: mie_distrib.png


   Классы
   ^^^^^^^

   .. automodule:: mstm_studio.contributions
      :members: Contribution, ConstantBackground, LinearBackground, LorentzBackground, LorentzPeak, GaussPeak, MieSingleSphere, MieLognormSpheres, MieLognormSpheresCached 



   .. [Kreibig_book1995] U. Kreibig, M. Vollmer, "Optical Properties of Metal Clusters" (1995) 553

