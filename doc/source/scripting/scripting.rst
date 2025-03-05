.. _scripting:  

Руководство: Скрипты  
=================  

Используйте Python-скрипты для полного контроля над вычислениями.  

.. note:: Простой способ указать путь к бинарному файлу MSTM в скрипте — использовать модуль `os`:  

    .. code-block:: python  

        import os  
        os.environ['MSTM_BIN'] = 'укажите путь к бинарному файлу mstm'  

    Путь по умолчанию: `'~/bin/mstm.x'`  

.. toctree::  
   materials  
   contribs  
   spheres  
   nearfield  
   nonspherical  
   fitting  
