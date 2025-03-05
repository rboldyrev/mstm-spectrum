.. _materials:

Материалы
---------

Материалы характеризуются показателем преломления, который является квадратным корнем (комплекснозначным) от макроскопической диэлектрической функции. 
В общем случае это спектральная функция, то есть:

.. math::

    n_c(\lambda) = n(\lambda) + i \cdot k(\lambda) = \sqrt{ \epsilon_1(\omega) + i \cdot \epsilon_2(\omega) }

таким образом,

.. math::

    n = \sqrt{(|\epsilon|+\epsilon_1)/2} \\
    k = \sqrt{(|\epsilon|-\epsilon_1)/2}

при :math:`\hbar \omega = 2 \pi \hbar c / \lambda`.

Постоянный материал
^^^^^^^^^^^^^^^^^

Материал с постоянным показателем преломления можно задать в первом аргументе конструктора (`file_name`)::

    >>> from mstm_studio.mstm_spectrum import Material
    >>> mat_glass = Material('1.5')
    >>> mat_glass.get_n(500)
    array(1.5)

Можно также задать комплексное значение::

    >>> mat_lossy = Material('3+1j')
    >>> mat_lossy.get_n(500)
    array(3.)
    >>> mat_lossy.get_k(500)
    array(1.)

Также можно использовать предопределенные названия: `air`, `water`, `glass`.

Загрузка из файла
^^^^^^^^^^^^^^^^^

Табличные данные о показателе преломления удобно хранить в файле. Заголовок файла должен содержать специальные метки: ``lambda n k``. Пример файла "etaGold.txt" можно найти в каталоге "nk" исходного кода.

Если файл "etaGold.txt" находится в той же директории, что и скрипт, его можно загрузить следующим образом:

.. literalinclude:: load_gold.py
   :lines: 3-5

Полученный график:

.. image:: loaded_gold.png

.. Note:: Расширенная база данных показателей преломления материалов <https://refractiveindex.info/>.

Материал из numpy массива
^^^^^^^^^^^^^^^^^^^^^^^^^

Данные о материале можно задать напрямую с помощью массива numpy (complex), передав `nk` или `eps`. 
Следующий пример показывает загрузку диэлектрической функции Drude:

.. literalinclude:: drude_gold.py
   :lines: 1-12

Члены класса Material
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mstm_studio.mstm_spectrum.Material
   :members:

Аналитическая формула для AuAg
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Серебро, золото и их сплавы можно задать с помощью аналитического выражения, предложенного в исследовании [Rioux2014]_.
Пример для сплава Au:Ag = 1:2:

.. literalinclude:: mat_au1ag2.py
   :lines: 1-5

Полученный график:

.. image:: mat_au1ag2.png

.. autoclass:: mstm_studio.alloy_AuAg.AlloyAuAg
   :members:

Материалы из RefractionIndex.Info 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Онлайн-база данных материалов RII [RII]_ <https:\refractiveindex.info> содержит 
диэлектрические функции сотен материалов. По состоянию на декабрь 2024 года 
она включает 445 различных материалов ("books" в терминологии RII), из которых 
163 подходят для расчетов в MSTM-Studio (то есть содержат табличные данные для 
длины волны от 300 до 800 нм).

Материалы можно загружать из локального дампа онлайн-базы, который можно 
скачать с официального сайта RII (About->Resources, прямая ссылка: 
<https://refractiveindex.info/download/database/rii-database-2024-12-31.zip>).

По умолчанию `mstm_studio` ищет архив `rii-database-*.zip` в домашнем каталоге 
и каталоге данных приложения. Точное местоположение можно задать в аргументе конструктора класса.

Пример: сравнение показателей преломления серебра у разных авторов::

.. literalinclude:: mat_rii_ag.py
   :lines: 1-29

Полученный график:

.. image:: mat_rii_ag.png

.. autoclass:: mstm_studio.rii_materials.RiiMaterial
   :members:

Коррекция размеров для диэлектрических функций
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Макроскопическая диэлектрическая функция, полученная для объемных образцов, 
может быть применена к наночастицам с осторожностью. 
Считается, что только частицы с радиусом более 10 нм можно считать корректными. 
Однако рассмотрение можно расширить до размеров ~2 нм, 
учитывая наиболее значимый эффект — уменьшение 
средней длины свободного пробега электронов из-за конечного размера наночастиц. 
Коррекция применяется к параметру :math:`\gamma` функции Drude, добавляя вклад:

.. math::

    \Delta \epsilon(\omega, D) = 
    {\epsilon}_{Drude, corr.}(\omega, D) - 
    {\epsilon}_{Drude}(\omega, D=\infty)

к экспериментальной диэлектрической функции, заданной в таблице.

Пример для золотой наночастицы радиусом 3 нм:

.. literalinclude:: size_correction.py
   :lines: 1-19

Полученный график:

.. image:: size_correction.png

В настоящее время реализованы коррекции для золота и серебра:

.. autoclass:: mstm_studio.diel_size_correction.SizeCorrectedGold
   :members:

.. autoclass:: mstm_studio.diel_size_correction.SizeCorrectedSilver
   :members:

Также доступен общий класс коррекции:

.. autoclass:: mstm_studio.diel_size_correction.SizeCorrectedMaterial
   :members:

.. [Rioux2014] D. Rioux, et al., "An Analytic Model for the Dielectric Function of Au, Ag, and their Alloys" Adv. Opt. Mater. (2014) *2* 176-182 <http://dx.doi.org/10.1002/adom.201300457>

.. [RII] M. N. Polyanskiy, "Refractiveindex.info database of optical constants" Sci. Data (2024) *11*, 94 <https://doi.org/10.1038/s41597-023-02898-2>
