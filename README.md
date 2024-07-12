# PatternMatchingDTW
Repositorio con codigo para la calibración de lamparas de comparación en longitudes de onda utilizando DTW para identificar los puntos de referencia de las lámparas con las señales de referencia


## Estructura


El código principal en donde se pueden realizar pruebas de lámparas contra datos de referencia se encuentra en `run.py`.


El archivo `utils.py` contiene un conjunto de funciones variadas útiles para que el código sea más sencillo al momento de ejecutar.


El archivo `requirements.txt` contiene un listado de las dependencias funcionales que tiene la aplicación.


En la carpeta `Lamps` se puede encontrar una colección de lámparas de distintas fuentes. Están agrupadas en carpetas acorde a su origen.


En la carpeta `Reference_Data` se almacenan los archivos de datos correspondientes a los espectros de ciertos materiales. Están agrupados en carpetas según su origen.


## Instalación de dependencias


Para preparar el entorno para la ejecución deben estar instaladas todas las librerías de las que depende la aplicación, para ello:
```
pip install requirements.txt
```


## Ejecutar


Para ejecutar el algoritmo de calibrado:
```
python run.py
```