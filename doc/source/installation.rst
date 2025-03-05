Установка
============

Исходный код
------------

Исходный код Python-обертки доступен на GitHub <https://github.com/lavakyan/mstm-spectrum>. 
Стабильная версия опубликована на PyPi <https://pypi.org/project/mstm-studio/>.

Исходный код MSTM не включен и должен быть получен с <https://scattport.org/index.php/light-scattering-software/multiple-particle-scattering/468-mstm>. 
MSTM Studio может быть запущен без бинарного файла MSTM, но с ограниченной функциональностью.

Для не сферических частиц (в настоящее время доступны только сфероиды) используется библиотека ScatterPy (см. :ref:`binding-scatterpy`).


Установка на Linux
------------------

Установка из PyPi:

``pip install mstm_studio``

На системах без прав root:

``pip install mstm_studio --user``

Запуск графического интерфейса:

``python -m mstm_studio``

Может потребоваться явно указать версию Python, т.е. использовать ``pip3`` и ``python3`` в вышеуказанных командах.


Связывание с MSTM
^^^^^^^^^^^^^^^^^

MSTM-studio будет искать бинарный файл ``mstm.x`` в директории ``~/bin``.

Это можно изменить, установив переменную окружения `MSTM_BIN`, например, в bash:

``export MSTM_BIN=~/my_compiled_mstm/mstm_v3.bin``


.. Примечание::   MSTM может быть скомпилирован с помощью gfortran следующим образом::
      
    gfortran  mpidefs-serial.f90 mstm-intrinsics-v3.0.f90 mstm-modules-v3.0.f90 mstm-main-v3.0.f90 -O2  -o mstm.x
   
   Это последовательная компиляция, для параллельной файл ``mpidefs-serial.f90`` должен быть заменен. Подробности см. на сайте MSTM.


Установка на Windows
--------------------

Проверенный способ — использование дистрибутива Python Anaconda <https://www.anaconda.com/>. 

1. Откройте "Anaconda Prompt". Появится новое окно терминала. 
2. Введите ``pip install mstm_studio``. Это может занять некоторое время, так как зависимый код будет загружен и установлен.
3. Проверьте графический интерфейс, введя ``python -m mstm_studio`` в Anaconda Prompt, 
   или проверьте работу сценариев Python, введя ``import mstm_studio`` в консоли Python.


Связывание с MSTM
^^^^^^^^^^^^^^^^^

4. Получите бинарный файл MSTM. Поместите его в какую-либо папку. 
5. Установите переменную окружения ``MSTM_BIN``, чтобы она указывала на бинарный файл. 
   Команда оболочки ``SETX MSTM_BIN="path_to_your_mstm_bin"`` 
   выполнит временную настройку, что полезно для создания ``*.cmd`` скриптов. 
   Постоянная настройка переменной окружения должна быть выполнена с помощью графического интерфейса, см., например, 
   <https://docs.oracle.com/en/database/oracle/r-enterprise/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html>.

.. Примечание:: Если вы пишете \*.cmd скрипт для запуска графического интерфейса, не забудьте обновить переменную ``PATH``, чтобы она указывала на дистрибутив Python. 
    Самый простой способ — ввести ``echo %PATH%`` в Anaconda Prompt и использовать вывод в вашем скрипте.
    Пример скрипта для запуска графического интерфейса:
    
    .. code-block:: cmd
    
        @ECHO OFF
        PATH=C:\ProgramData\Anaconda3;C:\ProgramData\Anaconda3\Library\mingw-w64\bin;C:\ProgramData\Anaconda3\Library\usr\bin;C:\ProgramData\Anaconda3\Library\bin;C:\ProgramData\Anaconda3\Scripts;C:\ProgramData\Anaconda3\bin;C:\ProgramData\Anaconda3\condabin;%PATH%
        set MSTM_BIN="C:\Users\L\Desktop\mstm_studio old\mstm-spectrum\mstm.exe"
        python.exe -m mstm_studio
        PAUSE
        
    Последняя команда (``PAUSE``) добавлена, чтобы предотвратить закрытие окна консоли после завершения программы.


.. _binding-scatterpy:

Связывание с ScatterPy
----------------------

Для расчета спектров экстинкции изолированных не сферических частиц может быть использована библиотека ScatterPy. Эта библиотека доступна на GitHub <https://github.com/TCvanLeth/ScatterPy> и в репозитории PyPi.

Установка из PyPi: ``pip install scatterpy`` или ``pip install scatterpy --user``


ScatterPy без Numba
^^^^^^^^^^^^^^^^^^^

ScatterPy требует библиотеку Numba для ускорения расчетов. Однако, возможно установить ScatterPy без Numba:

1. Скачайте исходный код ScatterPy.
2. Отредактируйте файл ``scatterpy/special.py``.
   Удалите строку:
   
   .. code-block:: python
   
        import numba as nb
   
   и добавьте строки:
   
   .. code-block:: python
       
       try:
           import numba as nb
       except ImportError:
           print('WARNING: Numba support is disabled in ScatterPy')


3. Соберите и установите: ``python setup.py install`` (Требуется setuptools и, возможно, другие dev-пакеты)