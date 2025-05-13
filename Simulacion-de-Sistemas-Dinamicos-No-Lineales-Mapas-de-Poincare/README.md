# Simulación de Sistemas Dinámicos No Lineales: Mapas de Poincaré (TFG)

**Autor:** Pedro Escudero Murcia
**Institución:** Trabajo Fin de Grado (Código FS23-38-FSC), Grado en Física, Universidad de Córdoba.

## Introducción

Este repositorio contiene el código fuente, los resultados de las simulaciones y el documento final asociados al Trabajo Fin de Grado (TFG) titulado "Simulación de un sistema dinámico no lineal: mapas de Poincaré".

El objetivo principal de este trabajo es estudiar los fundamentos de los sistemas dinámicos no lineales y el fenómeno del caos determinista. Para ello, se simula y analiza el comportamiento de dos sistemas físicos paradigmáticos:

1.  **Péndulo Forzado Amortiguado:** Un sistema **disipativo** cuyo comportamiento (periódico, caótico) se estudia variando la amplitud del forzamiento externo ($A$).
2.  **Péndulo Doble:** Un sistema **conservativo Hamiltoniano** cuya dinámica (regular, mixta, caótica) se explora en función de su energía total ($E$).

El análisis se centra en la visualización de trayectorias en el espacio de fases, la construcción e interpretación de **mapas de Poincaré** (estroboscópicos y de sección transversal) y la cuantificación de la sensibilidad a condiciones iniciales mediante la estimación del **exponente máximo de Lyapunov**.

## Sistemas Estudiados y Código

* **Péndulo Forzado Amortiguado:**
    * Ecuación de movimiento: $\frac{d^{2}\theta}{dt^{2}}+\gamma\frac{d\theta}{dt}+(\omega_{0}^{2}+2A \cos(\omega t))\sin \theta=0$
    * Script principal: `PenduloForzado.py`
    * Parámetros clave: Coeficiente de amortiguamiento ($\gamma$), frecuencia natural ($\omega_0$), amplitud de forzamiento ($A$), frecuencia de forzamiento ($\omega$).
* **Péndulo Doble:**
    * Descrito mediante formalismo Hamiltoniano. Ver TFG y `EcuacionesMovimiento.py` para detalles.
    * Script principal: `PenduloDoble.py`
    * Parámetros clave: Masa ($m$), longitud ($L$), gravedad ($g$), energía total ($E$).

## Estructura del Repositorio

* `PenduloForzado.py`: Script principal para simular el péndulo forzado.
* `PenduloDoble.py`: Script principal para simular el péndulo doble.
* `Funciones.py`: Contiene funciones comunes como el integrador Runge-Kutta (RK4) y la normalización de ángulos.
* `EcuacionesMovimiento.py`: Script auxiliar (usando `sympy`) para derivar simbólicamente las ecuaciones de Hamilton del péndulo doble.
* `Simulación_de_un_sistema_dinámico_no_lineal__mapas_de_Poincaré.pdf`: Documento completo del TFG.
* `requirements.txt`: Lista de dependencias de Python.
* `README.md`: Este archivo.
* `data/`: Carpeta que contiene los resultados de las simulaciones.
    * `penduloforzado/`: Resultados del péndulo forzado.
        * `A_baja/`, `A_media/`, `A_alta/`: Subcarpetas con resultados pre-calculados para los casos discutidos en el TFG (¡**Nota:** los valores exactos de A para estos casos deben consultarse/ajustarse en el TFG o el código!).
        * Contiene archivos `.csv` (datos de trayectoria, Poincaré, error, discrepancia, bifurcaciones) y `.png` (gráficas).
        * `Resultados.txt`: Archivo resumen de texto con métricas clave de la simulación para esa configuración específica.
    * `pendulodoble/`: Resultados del péndulo doble.
        * `E_baja/`, `E_media/`, `E_alta/`: Subcarpetas con resultados pre-calculados para las energías discutidas en el TFG (E=1.0 J, 15.0 J, 40.0 J).
        * Contiene archivos `.csv` (datos de trayectoria, Poincaré, error, discrepancia, energía) y `.png` (gráficas).
        * `Resultados.txt`: Archivo resumen de texto con métricas clave de la simulación para esa configuración específica.

## Instalación y Dependencias

1.  **Clonar el repositorio:**
    ```bash
    git clone [https://github.com/tu_usuario/Mapas-de-Poincare.git](https://github.com/tu_usuario/Mapas-de-Poincare.git)
    cd Mapas-de-Poincare
    ```
2.  **Entorno Virtual (Recomendado):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # En Windows: venv\Scripts\activate
    ```
3.  **Instalar dependencias:**
    ```bash
    pip install -r requirements.txt
    ```
    Las dependencias principales son `numpy` y `matplotlib`. `sympy` es necesario si se quiere ejecutar `EcuacionesMovimiento.py`.

## Ejecución de las Simulaciones

Puedes ejecutar las simulaciones principales directamente desde la terminal:

* **Péndulo Forzado:**
    ```bash
    python PenduloForzado.py
    ```
* **Péndulo Doble:**
    ```bash
    python PenduloDoble.py
    ```

**Para modificar los parámetros y explorar diferentes regímenes:**

1.  **Abre el script** correspondiente (`PenduloForzado.py` o `PenduloDoble.py`) en un editor de texto o IDE.
2.  **Localiza la "SECCIÓN 2: PARÁMETROS DEL SISTEMA"** dentro del script.
3.  **Modifica los valores** de los parámetros físicos y de simulación según desees:
    * En `PenduloForzado.py`: Cambia `CoeficienteAmortiguamiento`, `FrecuenciaNatural`, `Amplitud`, `FrecuenciaFuerzaExterna`, `CondicionesIniciales`, `TFinal`, `TamañoPaso`, etc. Para reproducir los casos $A_{baja}, A_{media}, A_{alta}$ del TFG, deberás ajustar el valor de `Amplitud` a los valores correspondientes (que deben especificarse en el TFG).
    * En `PenduloDoble.py`: Cambia `m`, `L`, `Gravedad`, `EnergiaDeseada`, `CondicionesIniciales` (se recomienda modificar solo $\theta_1(0)$, $\theta_2(0)$, $p_1(0)$ y dejar que $p_2(0)$ se calcule a partir de `EnergiaDeseada`), `TFinal`, `TamañoPaso`, etc.
4.  **Controla la ejecución:** Puedes activar/desactivar partes del análisis modificando las variables booleanas cerca del final de la Sección 2:
    * `calcular_errores`: Activa/desactiva el cálculo detallado del error RK4.
    * `generar_graficas`: Activa/desactiva la generación de archivos `.png`.
    * `realizar_analisis_discrepancia`: Activa/desactiva la simulación de la trayectoria perturbada y el cálculo del exponente de Lyapunov.
    * `realizar_diagrama_bifurcaciones` (solo en `PenduloForzado.py`): Activa/desactiva la generación del diagrama de bifurcaciones (¡puede ser computacionalmente intensivo!).
5.  **Guarda los cambios** en el script.
6.  **Ejecuta el script** modificado desde la terminal como se indicó anteriormente. Los nuevos resultados se guardarán en la carpeta `data/` (sobrescribiendo los anteriores si no se cambian los nombres de archivo o la estructura de carpetas).

## Interpretación de los Resultados

Al ejecutar los scripts, se generarán (si están activadas las opciones correspondientes):

* **Archivos `.csv`:** Contienen datos numéricos detallados. Pueden ser importados en otras herramientas para análisis posteriores.
    * `*_Datos.csv`: Coordenadas de la trayectoria a lo largo del tiempo.
    * `*_Poincare.csv`: Coordenadas $(\theta, \dot{\theta})$ o $(\theta_1, p_1)$ etc., en los cruces de la sección de Poincaré.
    * `*_ErrorRK4.csv`: Estimación del error local del método RK4 en cada paso.
    * `*_Discrepancia.csv`: Distancia entre la trayectoria original y la perturbada, y su logaritmo (usado para $\lambda_{max}$).
    * `*_Energia.csv` (Péndulo Doble): Valor de la energía Hamiltoniana y su desviación a lo largo del tiempo.
    * `*_Bifurcaciones.csv` (Péndulo Forzado): Pares (Amplitud, $\theta_{Poincare}$) para construir el diagrama.
* **Archivos `.png`:** Visualizaciones generadas con Matplotlib.
    * Diagramas de fase.
    * Mapas de Poincaré (a veces superpuestos o con densidad).
    * Evolución temporal de las variables.
    * Gráfica para la estimación del exponente de Lyapunov.
    * Diagrama de bifurcaciones (péndulo forzado).
* **Salida en Terminal y `Resultados.txt`:** Al finalizar la ejecución, el script imprime un resumen de los resultados clave en la terminal. Este mismo resumen se guarda en el archivo `Resultados.txt` dentro de la carpeta `data` correspondiente (o subcarpetas si se ejecutan casos específicos), proporcionando una visión rápida de:
    * Número de puntos de Poincaré encontrados.
    * Errores de conservación de energía (péndulo doble).
    * Errores estimados del método RK4.
    * Métricas de sensibilidad a condiciones iniciales (discrepancia, tiempo de duplicación, $\lambda_{max}$).
    * Correlación entre trayectorias (si se calculó discrepancia).

## Documento del TFG

Para una comprensión profunda del marco teórico, la derivación detallada de las ecuaciones, la justificación de los métodos numéricos y un análisis e interpretación completos de los resultados presentados, por favor consulta el documento PDF incluido en este repositorio:

`Simulación_de_un_sistema_dinámico_no_lineal__mapas_de_Poincaré.pdf`

---