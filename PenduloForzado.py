import numpy as np
import matplotlib.pyplot as plt
import os
import Funciones as Funciones
import subprocess
subprocess.run('clear', shell=True) #Para limpiar la terminal cada vez que se ejecute el programa.

"""
PÉNDULO FORZADO AMORTIGUADO
===========================
Simulación completa del péndulo forzado amortiguado, que incluye:
- Integración numérica mediante método Runge-Kutta de 4º orden
- Cálculo de puntos de Poincaré
- Análisis de errores numéricos
- Modo de discrepancia para comparación de trayectorias
"""

# =====================================================================
#           SECCIÓN 1: DEFINICIÓN DE FUNCIONES
# =====================================================================

def PenduloAmortiguadoForzado(Tiempo, Variables, CoeficienteAmortiguamiento=0.2, FrecuenciaNatural=1.0, Amplitud=0.15, FrecuenciaFuerzaExterna=2.0):
    """
    Define las ecuaciones de movimiento del péndulo forzado amortiguado.
    
    Entradas:
    Tiempo --> es el instante en que calculamos la función (s)
    Variables --> Es un vector fila con [Theta (rad), dTheta/dt (rad/s)]
    CoeficienteAmortiguamiento --> Coeficiente de amortiguamiento del sistema (1/s)
    FrecuenciaNatural --> Frecuencia de oscilación del péndulo sin fuerzas externas (rad/s)
    Amplitud --> Amplitud de la fuerza externa (1/s²)
    FrecuenciaFuerzaExterna --> Frecuencia de oscilación de la fuerza externa (rad/s)
    
    Salidas:
    Du1 -> Resultado de la primera ecuacion de primer grado (rad/s)
    Du2 -> Resultado de la segunda ecuación de segundo grado (rad/s²)
    
    -------------------------------------------------------------
    Como la ecuación del péndulo forzado amortiguado es de segundo orden,
    para aplicarle el método de Runge-Kutta tenemos que reducirla a un sistema 
    de dos ecuaciones de primer orden. Siendo
    Du1 = dTheta/dt  (rad/s)
    Du2 = d²Theta/dt²  (rad/s²)
    """
    
    Du1 = Variables[1]

    Du2 = -CoeficienteAmortiguamiento * Variables[1] - (FrecuenciaNatural**2 + 2* Amplitud * np.cos((FrecuenciaFuerzaExterna * Tiempo))) * np.sin(Variables[0])
    
    return np.array([Du1, Du2])

def PoincarePenduloForzado(Tiempo, Variables, PeriodoFuerzaExterna, NumeroPeriodosTransitorio, TamañoPaso, nombreArchivoDatos='data/penduloforzado/PenduloForzado_Datos.csv', nombreArchivoPoincare='data/penduloforzado/PenduloForzado_Poincare.csv'):
    """
    Genera la sección de Poincaré estroboscópica para el péndulo forzado.
    Selecciona el estado del sistema (theta, dtheta/dt) en tiempos t = nT, 
    donde T es el PeriodoFuerzaExterna, después de descartar un número 
    especificado de periodos transitorios.

    Encuentra el punto de la simulación más cercano a cada t = nT.
    
    Parámetros:
    -----------
    Tiempo: array
        Array con los tiempos de la simulación (s).
    Variables: array
        Array con las variables del sistema [theta (rad), dtheta/dt (rad/s)].
    PeriodoFuerzaExterna: float
        El periodo (T) de la fuerza externa (s).
    NumeroPeriodosTransitorio: int
        Número de periodos iniciales a descartar como transitorio.
    TamañoPaso: float
        El paso de tiempo usado en la simulación RungeKutta (s).
    nombreArchivoDatos: str
        Nombre del archivo CSV donde se guardarán los datos completos (tiempo, theta, dtheta/dt).
    nombreArchivoPoincare: str
        Nombre del archivo CSV donde se guardarán los puntos de Poincaré (theta, dtheta/dt).
        
    Devuelve:
    --------
    datos_completos: array
        Array con tiempo (s), theta (rad), dtheta/dt (rad/s).
    puntos_poincare: array
        Array con los puntos de Poincaré [theta (rad), dtheta/dt (rad/s)].
    """
    
    PuntosPoincare = []
    indices_poincare_usados = set() # Para evitar duplicados si varios nT caen cerca del mismo índice

    TiempoInicioMuestreo = NumeroPeriodosTransitorio * PeriodoFuerzaExterna
    
    # Determinar el primer múltiplo de T después del transitorio
    n_inicial = int(np.ceil(TiempoInicioMuestreo / PeriodoFuerzaExterna))
    if n_inicial <= NumeroPeriodosTransitorio:
        n_inicial = NumeroPeriodosTransitorio + 1

    n = n_inicial
    while True:
        t_target = n * PeriodoFuerzaExterna
        
        # Si el tiempo objetivo excede el tiempo de simulación, paramos
        if t_target > Tiempo[-1]:
            break
            
        # Encontrar el índice del tiempo de simulación más cercano a t_target
        idx = np.argmin(np.abs(Tiempo - t_target))
        
        # Verificar que este índice no se haya usado ya y que esté dentro de una tolerancia razonable (medio paso de tiempo)
        if idx not in indices_poincare_usados and np.abs(Tiempo[idx] - t_target) < TamañoPaso / 2:
            # Asegurarse de que el punto encontrado está realmente después del tiempo de inicio del muestreo
            if Tiempo[idx] >= TiempoInicioMuestreo:
                 PuntosPoincare.append(Variables[idx, :])
                 indices_poincare_usados.add(idx)
        
        n += 1 # Pasar al siguiente múltiplo del periodo

    # Guardar los datos completos en un fichero csv:
    datos_completos = np.column_stack((Tiempo, Variables))
    if nombreArchivoDatos is not None:
        np.savetxt(nombreArchivoDatos, datos_completos, delimiter='  ', header='Tiempo Theta dTheta_dt') # Cabecera actualizada
    
    puntos_poincare = np.array(PuntosPoincare)
    if puntos_poincare.size > 0:
        if nombreArchivoPoincare is not None:
            np.savetxt(nombreArchivoPoincare, puntos_poincare, delimiter='  ', header='Theta_Poincare dTheta_dt_Poincare') # Cabecera actualizada
    else:
        if nombreArchivoPoincare is not None:
            # Crear archivo vacío si no hay puntos
            open(nombreArchivoPoincare, 'w').close()
            print(f"Advertencia: No se encontraron puntos de Poincaré para guardar en {nombreArchivoPoincare}")

    return datos_completos, puntos_poincare

def calcular_error_rk4(func, y, t, h, args=()):
    """
    Calcula una estimación rigurosa del error local del método RK4.
    
    Parámetros:
    -----------
    func: callable
        La función que define el sistema de ecuaciones diferenciales.
    y: array_like
        El estado actual del sistema.
    t: float
        El tiempo actual.
    h: float
        El tamaño del paso.
    args: tuple, opcional
        Argumentos adicionales para la función func.
        
    Devuelve:
    --------
    error_estimado: array
        Estimación del error local en cada componente.
    """
    # Método de Richardson extrapolation
    # Calculamos la solución con paso h y con paso h/2
    
    # Solución con paso h
    k1 = np.array(func(t, y, *args))
    k2 = np.array(func(t + h/2, y + h*k1/2, *args))
    k3 = np.array(func(t + h/2, y + h*k2/2, *args))
    k4 = np.array(func(t + h, y + h*k3, *args))
    y_h = y + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    # Solución con paso h/2 (dos pasos)
    h_half = h / 2
    
    # Primer paso h/2
    k1 = np.array(func(t, y, *args))
    k2 = np.array(func(t + h_half/2, y + h_half*k1/2, *args))
    k3 = np.array(func(t + h_half/2, y + h_half*k2/2, *args))
    k4 = np.array(func(t + h_half, y + h_half*k3, *args))
    y_temp = y + h_half * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    # Segundo paso h/2
    k1 = np.array(func(t + h_half, y_temp, *args))
    k2 = np.array(func(t + h_half + h_half/2, y_temp + h_half*k1/2, *args))
    k3 = np.array(func(t + h_half + h_half/2, y_temp + h_half*k2/2, *args))
    k4 = np.array(func(t + h, y_temp + h_half*k3, *args))
    y_h_half = y_temp + h_half * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    # Estimación del error usando extrapolación de Richardson
    # Para RK4, el error es O(h^5), así que el factor es 2^4 = 16
    error_estimado = abs(y_h_half - y_h) / 15  # (2^4 - 1)
    
    return error_estimado

def calcular_discrepancia(Variables1, Variables2):
    """
    Calcula la discrepancia entre dos trayectorias.
    
    Parámetros:
    -----------
    Variables1: ndarray
        Primera trayectoria.
    Variables2: ndarray
        Segunda trayectoria con una pequeña perturbación inicial.
        
    Devuelve:
    --------
    discrepancia: ndarray
        Array con la discrepancia (distancia euclidiana) en cada paso de tiempo.
    """
    # Asegurar que ambas trayectorias tengan la misma longitud
    min_len = min(len(Variables1), len(Variables2))
    Variables1 = Variables1[:min_len]
    Variables2 = Variables2[:min_len]
    
    # Calcular la discrepancia como la distancia euclidiana en el espacio de fases
    discrepancia = np.sqrt(np.sum((Variables1 - Variables2)**2, axis=1))
    
    return discrepancia

def generar_diagrama_bifurcaciones(rango_amplitudes, parametros_sistema, parametros_simulacion, nombre_archivo=None):
    """
    Genera un diagrama de bifurcaciones para el péndulo forzado variando la amplitud
    de la fuerza externa.
    
    Parámetros:
    -----------
    rango_amplitudes: tuple
        Tupla (A_min, A_max, A_step) con el rango de amplitudes a explorar.
    parametros_sistema: tuple
        Tupla (CoeficienteAmortiguamiento, FrecuenciaNatural, FrecuenciaFuerzaExterna) 
        con los parámetros del sistema.
    parametros_simulacion: tuple
        Tupla (TInicial, TFinal, NumeroPeriodosTransitorio, TamañoPaso, CondicionesIniciales)
        con los parámetros de la simulación.
    nombre_archivo: str, opcional
        Nombre del archivo donde guardar el diagrama de bifurcaciones.
    
    Devuelve:
    --------
    A_bif: ndarray
        Array con los valores de amplitud para cada punto del diagrama.
    theta_bif: ndarray
        Array con los valores de theta para cada punto del diagrama.
    """
    # Desempaquetar parámetros
    A_min, A_max, A_step = rango_amplitudes
    CoeficienteAmortiguamiento, FrecuenciaNatural, FrecuenciaFuerzaExterna = parametros_sistema
    TInicial, TFinal, NumeroPeriodosTransitorio, TamañoPaso, CondicionesIniciales = parametros_simulacion
    
    # Periodo de la fuerza externa
    PeriodoFuerzaExterna = 2.0 * np.pi / FrecuenciaFuerzaExterna
    
    print("Generando diagrama de bifurcaciones para el péndulo forzado...")
    print(f"Rango de amplitudes: {A_min} a {A_max} con paso {A_step}")
    
    # Crear la secuencia de valores de amplitud a explorar
    valores_A = np.arange(A_min, A_max + A_step, A_step)
    num_valores = len(valores_A)
    print(f"Número total de amplitudes a calcular: {num_valores}")
    
    # Para guardar los resultados
    A_bif = []
    theta_bif = []
    
    for i, Amplitud in enumerate(valores_A):
        # Mostrar progreso cada cierto número de iteraciones
        if i % max(1, num_valores // 20) == 0 or i == num_valores - 1:
            print(f"Calculando para A = {Amplitud:.5f} ({i+1}/{num_valores}, {(i+1)/num_valores*100:.1f}%)")
        
        # Integrar el sistema para este valor de A
        Tiempo, Variables, _ = Funciones.RungeKutta(
            lambda Tiempo, Variables: PenduloAmortiguadoForzado(
                Tiempo,
                Variables,
                CoeficienteAmortiguamiento,
                FrecuenciaNatural,
                Amplitud,
                FrecuenciaFuerzaExterna
            ),
            CondicionesIniciales,
            TInicial,
            TFinal,
            TamañoPaso
        )
        
        # Normalizar ángulos
        Variables_norm = Funciones.NormalizarAngulos(Variables)
        
        # Obtener puntos de Poincaré
        _, puntos_poincare = PoincarePenduloForzado(
            Tiempo,
            Variables_norm,
            PeriodoFuerzaExterna,
            NumeroPeriodosTransitorio,
            TamañoPaso,
            None,  # No guardar CSV para cada iteración
            None   # No guardar CSV para cada iteración
        )
        
        # Guardar los valores de theta de los puntos de Poincaré
        if puntos_poincare is not None and len(puntos_poincare) > 0:
            for punto in puntos_poincare:
                A_bif.append(Amplitud)
                theta_bif.append(punto[0])  # theta
    
    # Convertir a arrays de NumPy
    A_bif = np.array(A_bif)
    theta_bif = np.array(theta_bif)
    
    print(f"Diagrama de bifurcaciones completado con {len(A_bif)} puntos.")
    
    # Guardar los resultados en un archivo CSV si se proporciona un nombre de archivo
    if nombre_archivo is not None:
        datos_bifurcacion = np.column_stack((A_bif, theta_bif))
        np.savetxt(nombre_archivo, datos_bifurcacion, delimiter='  ', 
                   header='Amplitud Theta_Poincare')
        print(f"Datos del diagrama de bifurcaciones guardados en: {nombre_archivo}")
    
    return A_bif, theta_bif

# =====================================================================
#           SECCIÓN 2: PARÁMETROS DEL SISTEMA
# =====================================================================

# Parámetros físicos del sistema
CoeficienteAmortiguamiento = 0.2    # Coeficiente de amortiguamiento del sistema (1/s)
FrecuenciaNatural = 1.0             # Frecuencia de oscilación del péndulo cuando no está sometido a ninguna fuerza (rad/s)
Amplitud = 1.3                     # Amplitud de la fuerza externa (s^-2)
FrecuenciaFuerzaExterna = 2         # Frecuencia de oscilación de la fuerza externa (rad/s)

# Tipos de amortiguamiento:
# - Subamortiguado: CoeficienteAmortiguamiento < FrecuenciaNatural 
# - Sobreamortiguado: CoeficienteAmortiguamiento > FrecuenciaNatural 
# - Amortiguamiento crítico: CoeficienteAmortiguamiento = FrecuenciaNatural 

# Periodo de la fuerza externa
PeriodoFuerzaExterna = 2.0 * np.pi / FrecuenciaFuerzaExterna  # (s)

# Declaro las variables para la simulación
# Intervalo de Tiempo:
TInicial = 0                        # Tiempo inicial de la simulación (s)
TFinal = 120                        # Tiempo total de simulación (s)
NumeroPeriodosTransitorio = 10      # Número de periodos a descartar

# Nombres de los archivos CSV de salida
carpetaDatos = "data/penduloforzado"
nombreArchivoDatos = f'{carpetaDatos}/PenduloForzado_Datos.csv'
nombreArchivoPoincare = f'{carpetaDatos}/PenduloForzado_Poincare.csv'
nombreArchivoError = f'{carpetaDatos}/PenduloForzado_ErrorRK4.csv'
nombreArchivoDiscrepancia = f'{carpetaDatos}/PenduloForzado_Discrepancia.csv'
nombreArchivoBifurcaciones = f'{carpetaDatos}/PenduloForzado_Bifurcaciones.csv'

# Asegurarse de que la carpeta existe
os.makedirs(carpetaDatos, exist_ok=True)

# Tamaño de paso para la integración numérica
TamañoPaso = PeriodoFuerzaExterna / 150  # Incremento de tiempo (s)

# Establezco las condiciones iniciales de la posición y la velocidad en un vector fila.
CondicionesIniciales = np.array([0.8, 0.8])  # [theta (rad), dtheta/dt (rad/s)]
print(f"Condiciones Iniciales: theta = {CondicionesIniciales[0]} rad, dtheta/dt = {CondicionesIniciales[1]} rad/s")

# Parámetros para análisis de discrepancia
realizar_analisis_discrepancia = True  # Bandera para activar/desactivar el análisis
perturbacion = 1e-10  # Tamaño de la perturbación para el análisis de discrepancia
CondicionesInicialesPerturbadas = CondicionesIniciales + np.array([perturbacion, 0])  # Perturbación en theta

# Flags para controlar la ejecución
calcular_errores = True        # Si es False, omite el cálculo de errores para una ejecución más rápida
generar_graficas = True        # Si es False, no generará gráficas (solo datos en archivos CSV)

# Parámetros para generar diagrama de bifurcaciones
realizar_diagrama_bifurcaciones = False  # Bandera para activar/desactivar el diagrama de bifurcaciones
rango_amplitudes = (0.1, 3.0, 0.01)  # Rango de amplitudes (A_min, A_max, A_step)

# =====================================================================
#           SECCIÓN 3: SIMULACIÓN Y ANÁLISIS
# =====================================================================

print("Iniciando integración numérica con el método Runge-Kutta 4...")

# Aplico el método de Runge-Kutta con las constantes definidas:
Tiempo, Variables, Derivadas = Funciones.RungeKutta(
    lambda Tiempo, Variables: PenduloAmortiguadoForzado(
        Tiempo, 
        Variables, 
        CoeficienteAmortiguamiento, 
        FrecuenciaNatural, 
        Amplitud, 
        FrecuenciaFuerzaExterna
    ), 
    CondicionesIniciales, 
    TInicial, 
    TFinal, 
    TamañoPaso
)
print(f"Integración completada: {len(Tiempo)} pasos de tiempo calculados")

# Calcular errores del método Runge-Kutta en puntos seleccionados
if calcular_errores:
    print("Calculando estimación rigurosa de errores del método RK4...")
    # Usar todos los puntos para el cálculo de error
    num_puntos_error = len(Tiempo)  # Usar todos los puntos disponibles
    indices_error = np.linspace(0, len(Tiempo)-1, num_puntos_error, dtype=int)
    errores_locales = np.zeros((num_puntos_error, 2))  # 2 componentes del sistema [theta, dtheta/dt]
    tiempos_error = np.zeros(num_puntos_error)

    for i, idx in enumerate(indices_error):
        
        t = Tiempo[idx]
        y = Variables[idx]
        
        # Calcular error local del método RK4
        error_local = calcular_error_rk4(
            lambda t, y: PenduloAmortiguadoForzado(
                t, y, 
                CoeficienteAmortiguamiento, 
                FrecuenciaNatural, 
                Amplitud, 
                FrecuenciaFuerzaExterna
            ),
            y, t, TamañoPaso
        )
        
        errores_locales[i] = error_local
        tiempos_error[i] = t

    # Calcular normas de los errores para cada componente
    norma_error_theta = np.linalg.norm(errores_locales[:, 0])
    norma_error_dtheta_dt = np.linalg.norm(errores_locales[:, 1])

    print(f"Norma del error en theta: {norma_error_theta:.8e} rad")
    print(f"Norma del error en dtheta/dt: {norma_error_dtheta_dt:.8e} rad/s")

# Normalizar los ángulos después de la integración y antes de la visualización
Variables_norm = Funciones.NormalizarAngulos(Variables)

# Análisis de discrepancia (si está habilitado)
if realizar_analisis_discrepancia:
    print("\nRealizando análisis de discrepancia (sensibilidad a condiciones iniciales)...")
    print(f"Perturbación aplicada: {perturbacion} rad en theta")
    
    # Integrar trayectoria con condición inicial perturbada
    Tiempo_pert, Variables_pert, _ = Funciones.RungeKutta(
        lambda Tiempo, Variables: PenduloAmortiguadoForzado(
            Tiempo, 
            Variables, 
            CoeficienteAmortiguamiento, 
            FrecuenciaNatural, 
            Amplitud, 
            FrecuenciaFuerzaExterna
        ), 
        CondicionesInicialesPerturbadas, 
        TInicial, 
        TFinal, 
        TamañoPaso
    )
    
    # Normalizar ángulos de la trayectoria perturbada
    Variables_pert_norm = Funciones.NormalizarAngulos(Variables_pert)
    
    # Calcular discrepancia entre trayectorias
    discrepancia = calcular_discrepancia(Variables_norm, Variables_pert_norm)
    
    # Calcular logaritmo natural de la discrepancia para análisis de exponente de Lyapunov
    log_discrepancia = np.log(discrepancia)
    log_discrepancia[np.isneginf(log_discrepancia)] = -20  # Reemplazar valores -inf
    
    # Guardar datos de discrepancia
    datos_discrepancia = np.column_stack((Tiempo, discrepancia, log_discrepancia))
    np.savetxt(nombreArchivoDiscrepancia, datos_discrepancia, delimiter='  ', 
               header='Tiempo(s) Discrepancia log(Discrepancia)')
    print(f"Datos de discrepancia guardados en: {nombreArchivoDiscrepancia}")

# Generar diagrama de bifurcaciones (si está habilitado)
if realizar_diagrama_bifurcaciones:
    print("\nGenerando diagrama de bifurcaciones...")
    parametros_sistema = (CoeficienteAmortiguamiento, FrecuenciaNatural, FrecuenciaFuerzaExterna)
    parametros_simulacion = (TInicial, TFinal, NumeroPeriodosTransitorio, TamañoPaso, CondicionesIniciales)
    A_bif, theta_bif = generar_diagrama_bifurcaciones(
        rango_amplitudes, 
        parametros_sistema, 
        parametros_simulacion, 
        nombreArchivoBifurcaciones
    )
    print(f"Diagrama de bifurcaciones guardado en: {nombreArchivoBifurcaciones}")

# =====================================================================
#           SECCIÓN 4: GUARDAR DATOS
# =====================================================================

print("\nGuardando datos en archivos CSV...")

# Guardar datos de error RK4
if calcular_errores:
    error_data = np.column_stack((
        tiempos_error, 
        errores_locales[:, 0], errores_locales[:, 1]
    ))
    np.savetxt(nombreArchivoError, error_data, delimiter='  ', 
               header='Tiempo(s) Error_theta(rad) Error_dtheta_dt(rad/s)')
    print(f"Estimaciones de error guardadas en: {nombreArchivoError}")

# Guardo la lista de resultados en archivos CSV y obtengo los datos para visualización directa:
print("Generando puntos de Poincaré...")
datos_completos, puntos_poincare = PoincarePenduloForzado(
    Tiempo, 
    Variables_norm, 
    PeriodoFuerzaExterna, # Pasar el periodo
    NumeroPeriodosTransitorio, # Pasar el número de periodos a descartar
    TamañoPaso, # Pasar el tamaño del paso
    nombreArchivoDatos,
    nombreArchivoPoincare
)
print(f"Datos completos guardados en: {nombreArchivoDatos}")
print(f"Puntos de Poincaré guardados en: {nombreArchivoPoincare}")

# =====================================================================
#           SECCIÓN 5: VISUALIZACIÓN
# =====================================================================

print("\nResumen de la simulación:")
print(f"Número de puntos de Poincaré encontrados: {len(puntos_poincare)}")
if calcular_errores:
    print(f"Error RK4 estimado (norma):")
    print(f"  - theta: {norma_error_theta:.8e} rad")
    print(f"  - dtheta/dt: {norma_error_dtheta_dt:.8e} rad/s")

if generar_graficas:
    # Configuración para todas las gráficas
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 14,
        'figure.titlesize': 20,
        'figure.figsize': (8, 6),
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })

    print("Generando y guardando gráficas...")

    # Diagrama de fases completo (sin filtrar por periodo transitorio)
    plt.figure(figsize=(8, 6))
    plt.scatter(Variables_norm[:, 0], Variables_norm[:, 1], s=1, color='blue', label='Trayectoria Completa')
    plt.xlabel('$\\theta$ (rad)') # Usar LaTeX para theta
    plt.ylabel('$\\dot{\\theta}$ (rad/s)') # Usar LaTeX para dtheta/dt
    plt.title('Diagrama de Fases del Péndulo Forzado')
    plt.legend()
    plt.grid(True)
    plt.axis('tight')
    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloForzado_DiagramaFases.png')
    plt.close()

    # Mapa de Poincaré (afectado por NumeroPeriodosTransitorio)
    if len(puntos_poincare) > 0:
        # Crear figura para el mapa de Poincaré superpuesto con diagrama de densidades
        plt.figure(figsize=(10, 8))
        
        # Graficar primero la trayectoria completa (diagrama de fases)
        plt.scatter(Variables_norm[:, 0], Variables_norm[:, 1], s=0.5, color='blue', alpha=0.6, label='Trayectoria')
        
        # Superponer el mapa de densidad de puntos de Poincaré
        if len(puntos_poincare) >= 10:
            # Calcular los límites para el eje Y basado en los datos de Poincaré
            y_min = min(puntos_poincare[:, 1])
            y_max = max(puntos_poincare[:, 1])
            # Añadir margen para mejor visualización
            y_margin = 0.1 * (y_max - y_min)
            y_min -= y_margin
            y_max += y_margin
            
            # Crear hexbin con colores que destaquen sobre el diagrama de fases
            hb = plt.hexbin(puntos_poincare[:, 0], puntos_poincare[:, 1], 
                          gridsize=50, cmap='plasma_r', bins='log',  # Escala logarítmica e invertida
                          mincnt=1, alpha=0.7)  # Alpha para ver el diagrama de fases debajo
            cb = plt.colorbar(hb, label='log10(N) puntos')
        else:
            # Si hay pocos puntos, mostrar solo puntos de Poincaré destacados
            plt.scatter(puntos_poincare[:, 0], puntos_poincare[:, 1], s=20, 
                       color='red', alpha=0.8, marker='o', label='Puntos de Poincaré')
            
        plt.xlabel('$\\theta$ (rad)')
        plt.ylabel('$\\dot{\\theta}$ (rad/s)')
        plt.title(f'Mapa de Poincaré sobre Diagrama de Fases\n(Descartando {NumeroPeriodosTransitorio} periodos iniciales)')
        plt.grid(True)
        plt.axis('tight')
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(f'{carpetaDatos}/PenduloForzado_MapaPoincare_Combinado.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Mantener la visualización de densidad en alta resolución por separado
        if len(puntos_poincare) >= 20:
            plt.figure(figsize=(10, 8))
            h = plt.hist2d(puntos_poincare[:, 0], puntos_poincare[:, 1], bins=60, 
                         cmap='YlOrRd', norm=plt.matplotlib.colors.LogNorm())
            plt.colorbar(h[3], label='Número de puntos')
            plt.xlabel('$\\theta$ (rad)')
            plt.ylabel('$\\dot{\\theta}$ (rad/s)')
            plt.title('Concentración de Puntos del Mapa de Poincaré (hist2d)')
            plt.grid(True)
            plt.axis('tight')
            plt.close()
            
            # Mostrar estadísticas de concentración de puntos en la terminal
            print("\n----- Análisis de Concentración de Puntos Poincaré -----")
            print(f"Total de puntos: {len(puntos_poincare)}")
            print(f"Rango θ: [{min(puntos_poincare[:,0]):.4f}, {max(puntos_poincare[:,0]):.4f}] rad")
            print(f"Rango θ̇: [{min(puntos_poincare[:,1]):.4f}, {max(puntos_poincare[:,1]):.4f}] rad/s")
        
        print(f"Mapa de Poincaré combinado guardado con {len(puntos_poincare)} puntos, mostrando densidad sobre diagrama de fases")

    # Theta y Velocidad vs tiempo (en la misma figura pero en subgráficas separadas)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Gráfica de Theta vs tiempo
    ax1.plot(Tiempo, Variables_norm[:, 0], color='blue')
    ax1.set_ylabel('$\\theta$ (rad)')
    ax1.set_title('Evolución de $\\theta$ y $\\dot{\\theta}$ en el tiempo')
    ax1.grid(True)
    ax1.axis('tight')

    # Gráfica de Velocidad vs tiempo
    ax2.plot(Tiempo, Variables_norm[:, 1], color='red')
    ax2.set_xlabel('Tiempo (s)')
    ax2.set_ylabel('$\\dot{\\theta}$ (rad/s)')
    ax2.grid(True)
    ax2.axis('tight')

    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloForzado_ThetaYVelocidadTiempo.png')
    plt.close()

    # Error del método RK4 para cada componente en subplots separados
    if calcular_errores:
        fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # Gráfica de error en theta
        axs[0].semilogy(tiempos_error, errores_locales[:, 0], 'b-')
        axs[0].set_ylabel('Error $\\theta$ (rad)')
        axs[0].set_title('Análisis Riguroso de Error del Método RK4')
        axs[0].grid(True)
        axs[0].axis('tight')
        
        # Gráfica de error en dtheta/dt
        axs[1].semilogy(tiempos_error, errores_locales[:, 1], 'r-')
        axs[1].set_xlabel('Tiempo (s)')
        axs[1].set_ylabel('Error $\\dot{\\theta}$ (rad/s)')
        axs[1].grid(True)
        axs[1].axis('tight')
        
        # En lugar de guardar la gráfica, mostrar los errores en terminal
        plt.close()

        # Mostrar estadísticas de error RK4 en la terminal
        max_error_theta = np.max(errores_locales[:, 0])
        max_error_dtheta_dt = np.max(errores_locales[:, 1])
        avg_error_theta = np.mean(errores_locales[:, 0])
        avg_error_dtheta_dt = np.mean(errores_locales[:, 1])

        print("\n----- Análisis de Error del Método RK4 -----")
        print(f"Error máximo en theta: {max_error_theta:.8e} rad")
        print(f"Error máximo en dtheta/dt: {max_error_dtheta_dt:.8e} rad/s")
        print(f"Error promedio en theta: {avg_error_theta:.8e} rad")
        print(f"Error promedio en dtheta/dt: {avg_error_dtheta_dt:.8e} rad/s")
        print(f"Norma RMS en theta: {norma_error_theta:.8e} rad")
        print(f"Norma RMS en dtheta/dt: {norma_error_dtheta_dt:.8e} rad/s")

    # Gráficas de discrepancia (si el análisis está habilitado)
    if realizar_analisis_discrepancia:
        # Gráfica de discrepancia vs tiempo (escala log)
        plt.figure(figsize=(10, 8))
        plt.semilogy(Tiempo, discrepancia, 'b-')
        plt.xlabel('Tiempo (s)')
        plt.ylabel('Discrepancia (escala log)')
        plt.title('Análisis de Sensibilidad a Condiciones Iniciales')
        plt.grid(True)
        plt.axis('tight')
        plt.tight_layout()
        plt.close()

        # Mostrar estadísticas de discrepancia en la terminal
        discrepancia_max = np.max(discrepancia)
        discrepancia_final = discrepancia[-1]
        tiempo_duplicacion = np.nan
        
        # Estimar tiempo de duplicación (tiempo para que la discrepancia aumente por un factor de 2)
        if discrepancia[0] > 0:
            idx_duplicacion = np.where(discrepancia > 2*discrepancia[0])[0]
            if len(idx_duplicacion) > 0:
                tiempo_duplicacion = Tiempo[idx_duplicacion[0]] - Tiempo[0]
        
        print("\n----- Análisis de Sensibilidad a Condiciones Iniciales -----")
        print(f"Perturbación inicial: {perturbacion:.2e} rad en θ")
        print(f"Discrepancia máxima: {discrepancia_max:.4e}")
        print(f"Discrepancia final (t={Tiempo[-1]}s): {discrepancia_final:.4e}")
        if not np.isnan(tiempo_duplicacion):
            print(f"Tiempo estimado de duplicación: {tiempo_duplicacion:.4f} s")
        
        # Gráfica para estimar el exponente de Lyapunov
        plt.figure(figsize=(12, 9))
        plt.plot(Tiempo, log_discrepancia, 'r-', label='Datos')
        
        # Estimar el exponente de Lyapunov con regresión lineal
        # Usar solo una parte específica de la curva, entre t=0s y t=25s
        tiempo_min_regresion = 0.0  # Omitir los primeros 5 segundos (transitorio)
        tiempo_max_regresion = 25.0  # Usar datos hasta los 25 segundos
        
        idx_min_regresion = np.searchsorted(Tiempo, tiempo_min_regresion)
        idx_max_regresion = np.searchsorted(Tiempo, tiempo_max_regresion)
        
        lyapunov_exponent = np.nan
        lyapunov_error = np.nan
        
        if idx_min_regresion < idx_max_regresion and idx_max_regresion > idx_min_regresion + 1:
            t_reg = Tiempo[idx_min_regresion:idx_max_regresion]
            log_disc_reg = log_discrepancia[idx_min_regresion:idx_max_regresion]
            
            # Filtrar valores no válidos (NaN o Inf)
            valid_indices = np.isfinite(log_disc_reg) & np.isfinite(t_reg)
            t_reg_clean = t_reg[valid_indices]
            log_disc_reg_clean = log_disc_reg[valid_indices]
            
            if len(t_reg_clean) > 1:
                # Realizar regresión lineal y obtener parámetros de ajuste y matriz de covarianza
                coef, cov = np.polyfit(t_reg_clean, log_disc_reg_clean, 1, cov=True)
                lyapunov_exponent = coef[0]  # Pendiente = exponente de Lyapunov
                intercept = coef[1]  # Intercepto
                
                # Calcular errores estándar del exponente e intercepto
                lyapunov_error = np.sqrt(cov[0, 0])  # Error de la pendiente
                intercept_error = np.sqrt(cov[1, 1])  # Error del intercepto
                
                poly1d_fn = np.poly1d(coef)
                y_predicted = poly1d_fn(t_reg_clean)
                
                # Calcular R²
                correlation_matrix = np.corrcoef(log_disc_reg_clean, y_predicted)
                correlation_xy = correlation_matrix[0, 1]
                r_squared = correlation_xy**2
                
                # Plotear la línea de ajuste sin label para la leyenda
                plt.plot(t_reg_clean, y_predicted, 'k--', label='Ajuste lineal')
                
                # Añadir ecuación de la regresión con errores y R² en un cuadro de texto
                ecuacion_text = f'$\\ln(Discrepancia) = ({lyapunov_exponent:.4f} \\pm {lyapunov_error:.4f}) \\cdot Tiempo + ({intercept:.4f} \\pm {intercept_error:.4f})$\n$R^2 = {r_squared:.4f}$'
                plt.text(0.02, 0.98, ecuacion_text, transform=plt.gca().transAxes, 
                        fontsize=14, verticalalignment='top', horizontalalignment='left',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.xlabel('Tiempo (s)')
        plt.ylabel('ln(Discrepancia)')
        plt.title('Estimación del Exponente de Lyapunov')
        plt.legend()
        plt.grid(True)
        plt.axis('tight')
        plt.tight_layout()
        plt.savefig(f'{carpetaDatos}/PenduloForzado_ExponenteLyapunov.png')
        plt.close()
        
       # Análisis de comparación de trayectorias
        # Normalizar la discrepancia respecto al máximo valor posible
        max_theta_range = 2*np.pi  # El rango máximo de theta es 2π
        max_omega_range = np.max(np.abs(Variables_norm[:, 1])) * 2  # Estimación del rango máximo de omega
        
        # Normalización de la discrepancia para determinar puntos significativos
        discrepancia_normalizada = discrepancia / np.sqrt(max_theta_range**2 + max_omega_range**2)
        
        # Puntos de divergencia significativa (donde la discrepancia supera cierto umbral)
        umbral_divergencia = 0.01  # 1% de la amplitud máxima
        puntos_divergencia = np.where(discrepancia_normalizada > umbral_divergencia)[0]
        
        if len(puntos_divergencia) > 0:
            tiempo_primera_divergencia = Tiempo[puntos_divergencia[0]]
        else:
            tiempo_primera_divergencia = float('inf')
            
        # Estadísticas de comparación de trayectorias
        print("\n----- Comparación de Trayectorias -----")
        print(f"Diferencia final en θ: {abs(Variables_norm[-1, 0] - Variables_pert_norm[-1, 0]):.6f} rad")
        print(f"Diferencia final en θ̇: {abs(Variables_norm[-1, 1] - Variables_pert_norm[-1, 1]):.6f} rad/s")
        
        if not np.isinf(tiempo_primera_divergencia):
            print(f"Tiempo de primera divergencia significativa: {tiempo_primera_divergencia:.4f} s")
        else:
            print("Las trayectorias no divergen significativamente durante la simulación")
            
        # Calcular correlación entre las trayectorias
        corr_theta = np.corrcoef(Variables_norm[:, 0], Variables_pert_norm[:, 0])[0, 1]
        corr_omega = np.corrcoef(Variables_norm[:, 1], Variables_pert_norm[:, 1])[0, 1]
        print(f"Correlación en θ: {corr_theta:.6f}")
        print(f"Correlación en θ̇: {corr_omega:.6f}")

    # Visualización del diagrama de bifurcaciones (si está habilitado)
    if realizar_diagrama_bifurcaciones and 'A_bif' in locals() and len(A_bif) > 0:
        print("\nVisualizando diagrama de bifurcaciones...")
        
        plt.figure(figsize=(12, 8))
        plt.scatter(A_bif, theta_bif, s=0.1, color='black')
        plt.xlabel('Amplitud de la fuerza externa (A)')
        plt.ylabel('$\\theta$ (rad)')
        plt.title('Diagrama de Bifurcaciones del Péndulo Forzado')
        plt.grid(True)
        plt.axis('tight')
        plt.tight_layout()
        plt.savefig(f'{carpetaDatos}/PenduloForzado_Bifurcaciones.png', dpi=300)
        plt.close()
        
        print(f"Diagrama de bifurcaciones guardado en: {carpetaDatos}/PenduloForzado_Bifurcaciones.png")

print(f"\n¡Simulación completada exitosamente!")
print(f"Todos los gráficos se han guardado en la carpeta: {carpetaDatos}")
print(f"Número de puntos de Poincaré encontrados: {len(puntos_poincare)}")
if calcular_errores:
    print(f"Error RK4 estimado (norma):")
    print(f"  - theta: {norma_error_theta:.8e} rad")
    print(f"  - dtheta/dt: {norma_error_dtheta_dt:.8e} rad/s")

if realizar_analisis_discrepancia and 'lyapunov_exponent' in locals() and not np.isnan(lyapunov_exponent):
    print(f"Exponente de Lyapunov estimado: {lyapunov_exponent:.4f} ± {lyapunov_error:.4f}")
    if lyapunov_exponent > 0:
        print("El sistema muestra comportamiento caótico (exponente de Lyapunov positivo)")
    else:
        print("El sistema no muestra comportamiento caótico (exponente de Lyapunov no positivo)")