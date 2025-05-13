import numpy as np
import matplotlib.pyplot as plt
import os
import Funciones as Funciones
import subprocess
subprocess.run('clear', shell=True)

"""
PÉNDULO DOBLE
============
Simulación completa del péndulo doble, que incluye:
- Integración numérica mediante método Runge-Kutta de 4º orden
- Cálculo de puntos de Poincaré
- Análisis de errores numéricos
- Análisis de discrepancia (sensibilidad a condiciones iniciales)
- Análisis de conservación de energía

Las ecuaciones del movimiento utilizadas corresponden exactamente a las derivadas
en el archivo EcuacionesMovimiento.py, que implementa el formalismo Hamiltoniano
según la documentación en docs/csm_ch06.pdf.

Sistema de unidades consistente:
- Tiempo: segundos (s)
- Longitud: metros (m)
- Masa: kilogramos (kg)
- Ángulos: radianes (rad)
- Momento angular: kg·m²/s
- Energía: Julios (J)
"""

# =====================================================================
#           SECCIÓN 1: DEFINICIÓN DE FUNCIONES
# =====================================================================

def PenduloDoble(Tiempo, Variables, m=1.0, L=1.0, Gravedad=9.8):
    """
    Define las ecuaciones de Hamilton del péndulo doble.
    
    IMPORTANTES:
    - Esta implementación sigue exactamente las ecuaciones derivadas en EcuacionesMovimiento.py
    - Se basa en el formalismo Hamiltoniano documentado en docs/csm_ch06.pdf
    
    Parámetros:
    -----------
    Tiempo: float
        Tiempo actual (s).
    Variables: array
        Estado actual del sistema [theta1, theta2, p1, p2].
    m: float
        Masa de cada partícula (kg).
    L: float
        Longitud de cada varilla (m).
    Gravedad: float
        Aceleración de la gravedad (m/s²).
        
    Retorna:
    --------
    array
        Las derivadas [dtheta1/dt, dtheta2/dt, dp1/dt, dp2/dt].
    """
    theta1, theta2, p1, p2 = Variables
    
    # Calcular diferencia de ángulos y sus funciones trigonométricas
    delta = theta1 - theta2
    sin_delta = np.sin(delta)
    cos_delta = np.cos(delta)
    
    # Denominador común D = 1 + sin^2(delta)
    D_denom = 1 + sin_delta**2
    
    # Ecuaciones derivadas en EcuacionesMovimiento.py:
    
    # Derivada de theta1 respecto al tiempo (dH/dp1)
    dtheta1_dt = (2*p1 - 2*p2*cos_delta) / (2*L**2*m*D_denom)
    
    # Derivada de theta2 respecto al tiempo (dH/dp2)
    dtheta2_dt = (-2*p1*cos_delta + 4*p2) / (2*L**2*m*D_denom)
    
    # Derivada de p1 respecto al tiempo (-dH/dtheta1)
    # Esta expresión fue simplificada parcialmente del original
    term1 = -2*L*Gravedad*m*np.sin(theta1)
    term2 = -(p1*p2*sin_delta) / (L**2*m*D_denom)
    term3 = ((p1**2 - 2*p1*p2*cos_delta + 2*p2**2) * sin_delta * cos_delta) / (L**2*m*D_denom**2)
    dp1_dt = term1 + term2 + term3
    
    # Derivada de p2 respecto al tiempo (-dH/dtheta2)
    # Esta expresión fue simplificada parcialmente del original
    term1 = -L*Gravedad*m*np.sin(theta2)
    term2 = (p1*p2*sin_delta) / (L**2*m*D_denom)
    term3 = -((p1**2 - 2*p1*p2*cos_delta + 2*p2**2) * sin_delta * cos_delta) / (L**2*m*D_denom**2)
    dp2_dt = term1 + term2 + term3
    
    return np.array([dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt])

def HamiltonPenduloDoble(Variables, m=1.0, L=1.0, Gravedad=9.8):
    """
    Calcula la energía Hamiltoniana del péndulo doble.
    
    IMPORTANTE: Esta implementación sigue exactamente el Hamiltoniano definido en EcuacionesMovimiento.py
    
    Parámetros:
    -----------
    Variables: array
        Estado actual del sistema [theta1, theta2, p1, p2].
    m: float
        Masa de cada partícula (kg).
    L: float
        Longitud de cada varilla (m).
    Gravedad: float
        Aceleración de la gravedad (m/s²).
    
    Retorna:
    --------
    float
        La energía Hamiltoniana total del sistema (J).
    """
    theta1, theta2, p1, p2 = Variables
    
    # Calcular diferencia de ángulos y su seno cuadrado
    delta = theta1 - theta2
    sin_delta_sq = np.sin(delta)**2
    cos_delta = np.cos(delta)
    
    # Denominador común D = 1 + sin_delta_sq
    D_denom = 1 + sin_delta_sq
    
    # Energía cinética T
    T_num = p1**2 + 2*p2**2 - 2*p1*p2*cos_delta
    T_den = 2*m*L**2 * D_denom
    T = T_num / T_den
    
    # Energía potencial U
    U = m*Gravedad*L*(3 - 2*np.cos(theta1) - np.cos(theta2))
    
    # Hamiltoniano total H = T + U
    H = T + U
    
    return H

def EcuacionSegundoGrado(EnergiaDeseada, CondicionesInicialesParciales, m, L, Gravedad, signo='+'):
    """
    Calcula el valor de p2 para una energía Hamiltoniana deseada.
    
    IMPORTANTE: Esta implementación es consistente con el Hamiltoniano
    definido en EcuacionesMovimiento.py
    
    Parámetros:
    -----------
    EnergiaDeseada: float
        Energía Hamiltoniana deseada para el sistema (J).
    CondicionesInicialesParciales: array
        Valores conocidos de [theta1, theta2, p1].
    m: float
        Masa de cada partícula (kg).
    L: float
        Longitud de cada varilla (m).
    Gravedad: float
        Aceleración de la gravedad (m/s²).
    signo: str
        '+' o '-' para seleccionar la solución deseada.
        
    Retorna:
    --------
    float
        Valor de p2 que cumple con la energía deseada.
    """
    theta1, theta2, p1 = CondicionesInicialesParciales
    
    # Calcular diferencia de ángulos y denominador D
    delta = theta1 - theta2
    sin_delta_sq = np.sin(delta)**2
    cos_delta = np.cos(delta)
    D_denom = 1 + sin_delta_sq
    
    # Energía potencial U según EcuacionesMovimiento.py
    U = m * Gravedad * L * (3 - 2 * np.cos(theta1) - np.cos(theta2))
    
    # La energía cinética tiene la forma: T = (a*p2^2 + b*p2 + c) / D_denom
    # Donde:
    # a = 1 / (m*L^2)
    # b = -p1*cos_delta / (m*L^2)
    # c = p1^2 / (2*m*L^2)
    
    # Coeficientes para la ecuación cuadrática
    a = 1.0 / (m * L**2 * D_denom)
    b = -p1 * cos_delta / (m * L**2 * D_denom)
    c = p1**2 / (2 * m * L**2 * D_denom) + U - EnergiaDeseada
    
    # Discriminante para verificar soluciones reales
    discriminante = b**2 - 4.0 * a * c
    
    if discriminante < 0:
        print(f"No hay soluciones reales. Discriminante = {discriminante}")
        return np.nan
    
    # Cálculo de las raíces
    if signo == '+':
        p2 = (-b + np.sqrt(discriminante)) / (2.0 * a)
    else:
        p2 = (-b - np.sqrt(discriminante)) / (2.0 * a)
    
    return p2

def PoincarePenduloDoble(Tiempo, Variables, Derivadas, nombreArchivoDatos=None, nombreArchivoPoincare=None, criterio_theta1=False, t_min=3.0):
    """
    Genera puntos de Poincaré para el péndulo doble de forma más robusta,
    utilizando interpolación lineal para mayor precisión en la sección de cruce.
    Se toman muestras cuando theta1=0 (con dtheta1/dt > 0) o theta2=0 (con dtheta2/dt > 0),
    según el criterio seleccionado, y evitando artefactos debidos a la normalización en +/-pi.

    Parámetros:
    -----------
    Tiempo: array
        Vector de tiempos de la simulación (s).
    Variables: ndarray
        Array con las variables (idealmente normalizadas a [-pi, pi] para los ángulos ANTES de llamar a esta función)
        del sistema en cada paso de tiempo. Columnas: [theta1, theta2, p1, p2]
    Derivadas: ndarray
        Array con las derivadas de las variables en cada paso de tiempo.
        Columnas: [dtheta1/dt, dtheta2/dt, dp1/dt, dp2/dt]
    nombreArchivoDatos: str, opcional
        Nombre del archivo para guardar todos los datos de la trayectoria.
    nombreArchivoPoincare: str, opcional
        Nombre del archivo para guardar los puntos de Poincaré (tiempo_cruce, theta1_cruce, theta2_cruce, p1_cruce, p2_cruce).
    criterio_theta1: bool
        Si es True, usa theta1=0 como criterio de sección; si es False, usa theta2=0.
    t_min: float
        Tiempo mínimo para comenzar a registrar puntos de Poincaré (s), para descartar el transitorio.

    Retorna:
    --------
    datos_completos: ndarray
        Datos completos de la trayectoria guardados en el formato [t, theta1, theta2, p1, p2].
    puntos_poincare_para_retorno: ndarray
        Puntos de Poincaré guardados en el archivo, formato: [t_cruce, th1_interp, th2_interp, p1_interp, p2_interp].
        Si no se guardan en archivo, este array puede estar vacío o contener los puntos.
        La función devuelve los puntos interpolados (solo estado) para uso directo si es necesario: [th1, th2, p1, p2]
    """
    # Índice de la columna de ángulo que usaremos como criterio (0 para theta1, 1 para theta2)
    idx_criterio_angulo = 0 if criterio_theta1 else 1
    # Índice de la columna de velocidad angular correspondiente para la condición de cruce
    idx_criterio_velocidad = idx_criterio_angulo

    PuntosPoincare_lista = [] # Lista para almacenar los puntos de Poincaré interpolados [t_cruce, th1, th2, p1, p2]

    for i in range(1, len(Tiempo)):
        if Tiempo[i-1] >= t_min: # Empezar a buscar cruces después del tiempo mínimo
            # Variables en el paso anterior (i-1) y actual (i)
            t_prev, t_curr = Tiempo[i-1], Tiempo[i]
            vars_prev, vars_curr = Variables[i-1, :], Variables[i, :]
            
            # Ángulo de criterio en el paso anterior y actual
            theta_criterio_prev = vars_prev[idx_criterio_angulo]
            theta_criterio_curr = vars_curr[idx_criterio_angulo]

            # Condición 1: Verificar si el ángulo de criterio cruza por cero (cambio de signo)
            # Y que el ángulo no sea ya cero en el punto previo para evitar dobles conteos si Δt es muy grande
            if theta_criterio_prev * theta_criterio_curr <= 0 and theta_criterio_prev != 0:
                # Condición 2: Verificar que el salto angular es pequeño (< pi)
                # Esto ayuda a distinguir un cruce genuino por cero de un "salto"
                # debido a la normalización en los límites +/- pi.
                if abs(theta_criterio_curr - theta_criterio_prev) < np.pi:
                    # Condición 3: Verificar que la velocidad angular sea positiva (o la dirección deseada)
                    # Se usa la velocidad en el punto 'i' (final del intervalo) o se podría interpolar también.
                    # Usar Derivadas[i, idx_criterio_velocidad] es más simple.
                    if Derivadas[i, idx_criterio_velocidad] > 0:
                        
                        # --- Inicio de la Interpolación Lineal ---
                        # Calcular el factor de interpolación 's'
                        # s = (target_value - value_prev) / (value_curr - value_prev)
                        # target_value para el ángulo de criterio es 0.
                        denominador_s = theta_criterio_curr - theta_criterio_prev
                        if abs(denominador_s) < 1e-12: # Evitar división por cero si los ángulos son idénticos
                            s = 0.5 # O tomar el promedio, o el punto actual. Si hay cruce, no deberían ser idénticos.
                        else:
                            s = -theta_criterio_prev / denominador_s
                        
                        # Asegurar que s esté en [0, 1] para una interpolación válida dentro del intervalo
                        s = np.clip(s, 0.0, 1.0)

                        # Interpolar todas las variables de estado y el tiempo
                        t_cruce = t_prev + s * (t_curr - t_prev)
                        vars_interp = vars_prev + s * (vars_curr - vars_prev)
                        
                        # Añadir el punto interpolado (tiempo y estado completo)
                        PuntosPoincare_lista.append(np.concatenate(([t_cruce], vars_interp)))
                        # --- Fin de la Interpolación Lineal ---

    puntos_poincare_array = np.array(PuntosPoincare_lista)

    # Preparar datos completos para guardar (Tiempo y Variables)
    datos_completos = np.column_stack((Tiempo, Variables))

    if nombreArchivoDatos:
        try:
            np.savetxt(nombreArchivoDatos, datos_completos, delimiter='  ',
                       header='Tiempo(s) Theta1(rad) Theta2(rad) P1(kg·m²/s) P2(kg·m²/s)')
        except Exception as e:
            print(f"Error al guardar datos completos en {nombreArchivoDatos}: {e}")

    puntos_poincare_para_retorno = np.array([])
    if puntos_poincare_array.size > 0:
        # Para el archivo, guardamos [t_cruce, th1_interp, th2_interp, p1_interp, p2_interp]
        if nombreArchivoPoincare:
            try:
                np.savetxt(nombreArchivoPoincare, puntos_poincare_array, delimiter='  ',
                           header='Tiempo_Cruce(s) Theta1_interp(rad) Theta2_interp(rad) P1_interp(kg·m²/s) P2_interp(kg·m²/s)')
            except Exception as e:
                 print(f"Error al guardar puntos de Poincaré en {nombreArchivoPoincare}: {e}")
        # Para el retorno, usualmente solo se necesitan las variables de estado [th1, th2, p1, p2] del cruce
        puntos_poincare_para_retorno = puntos_poincare_array[:, 1:] 
    else:
        if nombreArchivoPoincare:
             try:
                 with open(nombreArchivoPoincare, 'w') as f:
                     f.write('# No se encontraron puntos de Poincaré con la lógica de interpolación.\n')
             except Exception as e:
                 print(f"Error al crear archivo Poincaré en {nombreArchivoPoincare}: {e}")
        print("Advertencia: No se encontraron puntos de Poincaré con la lógica de interpolación.")

    return datos_completos, puntos_poincare_para_retorno

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
        
    Retorna:
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
        
    Retorna:
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

# =====================================================================
#           SECCIÓN 2: PARÁMETROS DEL SISTEMA
# =====================================================================

# Parámetros físicos del sistema:
m = 1.0                 # Masa de las partículas (kg)
L = 1.0                 # Longitud de las varillas (m)
Gravedad = 9.8          # Aceleración de la gravedad (m/s²)

# Parámetros de simulación:
# Intervalo de Tiempo:
TInicial = 0                 # Tiempo inicial de la simulación (s)
TFinal = 120                 # Tiempo final de la simulación (s)
TamañoPaso = 0.001           # Tamaño de paso para integración numérica (s)

# Nombres de los archivos CSV de salida
carpetaDatos = "data/pendulodoble"
nombreArchivoDatos = f'{carpetaDatos}/PenduloDoble_Datos.csv'
nombreArchivoPoincare = f'{carpetaDatos}/PenduloDoble_Poincare.csv'
nombreArchivoEnergia = f'{carpetaDatos}/PenduloDoble_Energia.csv'
nombreArchivoError = f'{carpetaDatos}/PenduloDoble_ErrorRK4.csv'
nombreArchivoDiscrepancia = f'{carpetaDatos}/PenduloDoble_Discrepancia.csv'

# Parámetros para la función PoincarePenduloDoble
criterio_theta1 = False       # Si es True usa theta1=0 como criterio; si es False usa theta2=0
t_min = 3.0                  # Tiempo mínimo para registrar puntos de Poincaré (s)

# Asegurarse de que la carpeta existe
os.makedirs(carpetaDatos, exist_ok=True)

# Energía del sistema para encontrar los puntos en el plano de Poincaré
EnergiaDeseada = 1          # Energía deseada para las condiciones iniciales (J)

# Establezco las condiciones inciales de la posición y la velocidad en un vector fila.
# CondicionesIniciales = [theta1, theta2, p1, p2]
# donde: theta1, theta2 están en radianes (rad)
#        p1, p2 son los momentos angulares generalizados (kg·m²/s)
CondicionesIniciales = np.array([0.2, 0.0, 0.0, 0.0]) 

# Calcular p2 para obtener la energía deseada
CondicionesIniciales[3] = EcuacionSegundoGrado(
    EnergiaDeseada, CondicionesIniciales[:3],
    m, L, Gravedad, signo='+'
)
print("\n===== PÉNDULO DOBLE (E={:.1f} J) =====".format(EnergiaDeseada))
print(f"Condiciones iniciales: θ₁={CondicionesIniciales[0]:.3f} rad, θ₂={CondicionesIniciales[1]:.3f} rad, p₁={CondicionesIniciales[2]:.3f} kg·m²/s, p₂={CondicionesIniciales[3]:.3f} kg·m²/s")

if np.isnan(CondicionesIniciales[3]):
    print("¡ADVERTENCIA! No se encontró solución real que cumpla con dot_theta2 > 0.")
    print("Verifique que la energía deseada y las condiciones iniciales son físicamente consistentes.")
    exit(1)

# Parámetros para análisis de discrepancia
realizar_analisis_discrepancia = True  # Bandera para activar/desactivar el análisis
perturbacion = 1e-10  # Tamaño de la perturbación para el análisis de discrepancia
CondicionesInicialesPerturbadas = CondicionesIniciales + np.array([perturbacion, 0, 0, 0])  # Perturbación en theta1
lyapunov_exponent_pd = np.nan  # Inicializar variable para el exponente de Lyapunov

# Flags para controlar la ejecución
calcular_errores = True        # Si es False, omite el cálculo de errores para una ejecución más rápida
generar_graficas = True        # Si es False, no generará gráficas (solo datos en archivos CSV)

# =====================================================================
#           SECCIÓN 3: SIMULACIÓN Y ANÁLISIS
# =====================================================================

print(f"1. Integrando ecuaciones de movimiento (t={TInicial}s → {TFinal}s, h={TamañoPaso}s)...")

# Aplico el método de Runge-Kutta pasando las constantes definidas:
Tiempo, Variables, Derivadas = Funciones.RungeKutta(
    lambda Tiempo, Variables: PenduloDoble(Tiempo, Variables, m, L, Gravedad), 
    CondicionesIniciales, 
    TInicial, 
    TFinal, 
    TamañoPaso
)

if calcular_errores:
    # Calcular errores del método Runge-Kutta en puntos seleccionados
    print(f"2. Calculando estimación de error RK4...")
    # Usar todos los puntos para el cálculo de error
    num_puntos_error = len(Tiempo)  # Usar todos los puntos disponibles
    indices_error = np.linspace(0, len(Tiempo)-1, num_puntos_error, dtype=int)
    errores_locales = np.zeros((num_puntos_error, 4))  # 4 componentes del sistema
    tiempos_error = np.zeros(num_puntos_error)
    
    for i, idx in enumerate(indices_error):
        t = Tiempo[idx]
        y = Variables[idx]
        
        # Calcular error local del método RK4
        error_local = calcular_error_rk4(
            lambda t, y: PenduloDoble(t, y, m, L, Gravedad),
            y, t, TamañoPaso
        )
        
        errores_locales[i] = error_local
        tiempos_error[i] = t

    # Calcular normas de los errores para cada componente
    norma_error_theta1 = np.linalg.norm(errores_locales[:, 0])
    norma_error_theta2 = np.linalg.norm(errores_locales[:, 1])
    norma_error_p1 = np.linalg.norm(errores_locales[:, 2])
    norma_error_p2 = np.linalg.norm(errores_locales[:, 3])

# Normalizar los ángulos después de la integración y antes de la visualización
Variables_norm = Funciones.NormalizarAngulos(Variables)

# Calcular la energía Hamiltoniana en cada paso de tiempo
print("3. Calculando energía y análisis de conservación...")
EnergiaHamiltoniana = np.array([HamiltonPenduloDoble(Variables[i, :], m, L, Gravedad) for i in range(len(Tiempo))])

# Error en la conservación de energía como medida adicional de precisión
error_energia = EnergiaHamiltoniana - EnergiaHamiltoniana[0]
error_energia_relativo = error_energia / EnergiaHamiltoniana[0]

# Análisis de discrepancia (si está habilitado)
if realizar_analisis_discrepancia:
    print("4. Realizando análisis de sensibilidad a condiciones iniciales...")
    
    # Integrar trayectoria con condición inicial perturbada
    Tiempo_pert, Variables_pert, Derivadas_pert = Funciones.RungeKutta(
        lambda Tiempo, Variables: PenduloDoble(Tiempo, Variables, m, L, Gravedad), 
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
    
    # Calcular el exponente de Lyapunov mediante regresión lineal
    # Usar solo una parte específica de la curva, entre t=0.5s y t=25s
    tiempo_min_regresion = 0.5  # Omitir los primeros 0.5 segundos (transitorio)
    tiempo_max_regresion = 25.0  # Usar datos hasta los 25 segundos
    
    idx_min_regresion = np.searchsorted(Tiempo, tiempo_min_regresion)
    idx_max_regresion = np.searchsorted(Tiempo, tiempo_max_regresion)
    
    lyapunov_exponent_pd = np.nan
    lyapunov_error_pd = np.nan
    
    if idx_min_regresion < idx_max_regresion and idx_max_regresion > idx_min_regresion + 1:  # Necesitamos al menos 2 puntos
        t_reg = Tiempo[idx_min_regresion:idx_max_regresion]
        log_disc_reg = log_discrepancia[idx_min_regresion:idx_max_regresion]
        
        # Filtrar NaNs o Infs que puedan haber en log_disc_reg
        valid_indices = np.isfinite(log_disc_reg) & np.isfinite(t_reg)
        t_reg_clean = t_reg[valid_indices]
        log_disc_reg_clean = log_disc_reg[valid_indices]
        
        if len(t_reg_clean) > 1:
            # Realizar regresión lineal y obtener parámetros de ajuste
            coef, cov = np.polyfit(t_reg_clean, log_disc_reg_clean, 1, cov=True)
            lyapunov_exponent_pd = coef[0]  # Pendiente = exponente de Lyapunov
            
            # Calcular error estándar del exponente (raíz de la varianza de la pendiente)
            lyapunov_error_pd = np.sqrt(cov[0, 0])
    
    # Guardar datos de discrepancia
    datos_discrepancia = np.column_stack((Tiempo, discrepancia, log_discrepancia))
    np.savetxt(nombreArchivoDiscrepancia, datos_discrepancia, delimiter='  ', 
               header='Tiempo(s) Discrepancia log(Discrepancia)')

# =====================================================================
#           SECCIÓN 4: GUARDAR DATOS Y CÁLCULO DE PUNTOS DE POINCARÉ
# =====================================================================

print("5. Guardando datos y calculando puntos de Poincaré...")

# Guardar datos de error RK4 si se calcularon
if calcular_errores:
    error_data = np.column_stack((
        tiempos_error, 
        errores_locales[:, 0], errores_locales[:, 1],
        errores_locales[:, 2], errores_locales[:, 3]
    ))
    np.savetxt(nombreArchivoError, error_data, delimiter='  ', 
               header='Tiempo(s) Error_theta1(rad) Error_theta2(rad) Error_p1(kg·m²/s) Error_p2(kg·m²/s)')

# Guardar datos de energía
energia_data = np.column_stack((Tiempo, EnergiaHamiltoniana, error_energia, error_energia_relativo))
np.savetxt(nombreArchivoEnergia, energia_data, delimiter='  ', 
           header='Tiempo(s) Energia(J) Error_Absoluto(J) Error_Relativo')

# Cálculo de puntos de Poincaré
datos_completos, puntos_poincare = PoincarePenduloDoble(
    Tiempo, 
    Variables_norm, 
    Derivadas, 
    nombreArchivoDatos,
    nombreArchivoPoincare,
    criterio_theta1,
    t_min
)

# =====================================================================
#           SECCIÓN 5: VISUALIZACIÓN
# =====================================================================

# Si se solicitan gráficas, generarlas
if generar_graficas:
    # Configuración para todas las gráficas
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.titlesize': 16,
        'figure.figsize': (8, 6),
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })

    print("6. Generando gráficas seleccionadas...")

    # 1. Gráfica combinada de todas las variables vs tiempo con subplots separados
    plt.figure(figsize=(12, 10))
    fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    axs[0].plot(Tiempo, Variables_norm[:, 0], color='blue')
    axs[0].set_ylabel('$\\theta_1$ (rad)')
    axs[0].set_title('Evolución Temporal de Variables del Péndulo Doble')
    axs[0].grid(True)
    axs[1].plot(Tiempo, Variables_norm[:, 1], color='red')
    axs[1].set_ylabel('$\\theta_2$ (rad)')
    axs[1].grid(True)
    axs[2].plot(Tiempo, Variables[:, 2], color='green')
    axs[2].set_ylabel('$p_1$ (kg·m²/s)')
    axs[2].grid(True)
    axs[3].plot(Tiempo, Variables[:, 3], color='purple')
    axs[3].set_ylabel('$p_2$ (kg·m²/s)')
    axs[3].set_xlabel('Tiempo (s)')
    axs[3].grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)
    plt.savefig(f'{carpetaDatos}/PenduloDoble_TodasVariables_vs_Tiempo.png')
    plt.close()

    # 1. Diagrama de Fase: Momento1 vs Theta1
    plt.figure()
    plt.scatter(Variables_norm[:, 0], Variables[:, 2], s=1, color='blue', label='Trayectoria')
    plt.xlabel('$\\theta_1$ (rad)')
    plt.ylabel('$p_1$ (kg·m²/s)')
    plt.title('Diagrama de Fase ($p_1$ vs $\\theta_1$)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloDoble_Fase_p1_theta1.png')
    plt.close()

    # 2. Diagrama de Fase: Momento2 vs Theta2
    plt.figure()
    plt.scatter(Variables_norm[:, 1], Variables[:, 3], s=1, color='blue', label='Trayectoria')
    plt.xlabel('$\\theta_2$ (rad)')
    plt.ylabel('$p_2$ (kg·m²/s)')
    plt.title('Diagrama de Fase ($p_2$ vs $\\theta_2$)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloDoble_Fase_p2_theta2.png')
    plt.close()

    # 5. Mapa de Poincaré (theta2 vs p2)
    plt.figure()
    if puntos_poincare.size > 0:
        theta2_poincare = puntos_poincare[:, 1]
        p2_poincare = puntos_poincare[:, 3]
        plt.scatter(theta2_poincare, p2_poincare, s=10, color='red', label=f'Sección $\\theta_2=0, \\dot{{\\theta}}_2>0$')
    else:
        plt.text(0.5, 0.5, 'No se encontraron puntos de Poincaré', 
                 horizontalalignment='center', verticalalignment='center', 
                 transform=plt.gca().transAxes)
    plt.xlabel('$\\theta_2$ (rad)')
    plt.ylabel('$p_2$ (kg·m²/s)')
    plt.title(f'Mapa de Poincaré (E ≈ {EnergiaHamiltoniana[0]:.2f} J)')
    plt.legend()
    plt.grid(True)
    plt.xlim(-np.pi, np.pi)
    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloDoble_Poincare_theta2_p2.png')
    plt.close()

    # 6. Mapa de Poincaré (theta1 vs p1)
    plt.figure()
    if puntos_poincare.size > 0:
        theta1_poincare = puntos_poincare[:, 0]
        p1_poincare = puntos_poincare[:, 2]
        plt.scatter(theta1_poincare, p1_poincare, s=10, color='red', label=f'Sección $\\theta_2=0, \\dot{{\\theta}}_2>0$')
    else:
        plt.text(0.5, 0.5, 'No se encontraron puntos de Poincaré', 
                horizontalalignment='center', verticalalignment='center', 
                transform=plt.gca().transAxes)
    plt.xlabel('$\\theta_1$ (rad)')
    plt.ylabel('$p_1$ (kg·m²/s)')
    plt.title(f'Mapa de Poincaré (E ≈ {EnergiaHamiltoniana[0]:.2f} J)')
    plt.legend()
    plt.grid(True)
    plt.xlim(-np.pi, np.pi)
    plt.tight_layout()
    plt.savefig(f'{carpetaDatos}/PenduloDoble_Poincare_theta1_p1.png')
    plt.close()

    # Gráfica Lyapunov
    if realizar_analisis_discrepancia and not np.isnan(lyapunov_exponent_pd):
        plt.figure(figsize=(10, 8))
        plt.plot(Tiempo, log_discrepancia, 'r-')
        
        # Usar los valores ya calculados en la sección 3
        if not np.isnan(lyapunov_exponent_pd):
            # Obtener los datos para la visualización de la línea de ajuste
            idx_max_regresion = np.argmax(Tiempo > tiempo_max_regresion)
            if idx_max_regresion == 0 and Tiempo[-1] <= tiempo_max_regresion:
                idx_max_regresion = len(Tiempo)
            
            if idx_max_regresion > 1:
                t_reg = Tiempo[:idx_max_regresion]
                log_disc_reg = log_discrepancia[:idx_max_regresion]
                
                valid_indices = np.isfinite(log_disc_reg) & np.isfinite(t_reg)
                t_reg_clean = t_reg[valid_indices]
                
                if len(t_reg_clean) > 1:
                    # Crear la línea de tendencia con el exponente ya calculado
                    poly1d_fn = np.poly1d([lyapunov_exponent_pd, log_disc_reg_clean[0] - lyapunov_exponent_pd * t_reg_clean[0]])
                    plt.plot(t_reg_clean, poly1d_fn(t_reg_clean), 'k--', 
                            label=f'$\\lambda_{{max}} = {lyapunov_exponent_pd:.4f} \\pm {lyapunov_error_pd:.4f}$')
        
        plt.xlabel('Tiempo (s)')
        plt.ylabel('ln(Discrepancia)')
        plt.title('Estimación del Exponente de Lyapunov (Péndulo Doble)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'{carpetaDatos}/PenduloDoble_ExponenteLyapunov.png')
        plt.close()

# =====================================================================
#           SECCIÓN 6: RESUMEN DE RESULTADOS
# =====================================================================

print("\n======== RESULTADOS DE LA SIMULACIÓN ========")
print(f"Puntos de Poincaré: {len(puntos_poincare) if puntos_poincare.size > 0 else 0}")

# Información de conservación de energía
print("\n----- Conservación de Energía -----")
print(f"Energía inicial: {EnergiaHamiltoniana[0]:.8f} J")
print(f"Energía final: {EnergiaHamiltoniana[-1]:.8f} J")
print(f"Variación absoluta: {np.abs(EnergiaHamiltoniana[-1] - EnergiaHamiltoniana[0]):.8e} J")

# Error en conservación de energía
error_energia_relativo_max = np.max(np.abs(error_energia_relativo)) * 100
error_energia_relativo_avg = np.mean(np.abs(error_energia_relativo)) * 100

print(f"Error relativo máximo: {error_energia_relativo_max:.8e} %")
print(f"Error relativo promedio: {error_energia_relativo_avg:.8e} %")

# RK4 error analysis
if calcular_errores:
    print("\n----- Análisis de Error del Método RK4 -----")
    max_error_theta1 = np.max(errores_locales[:, 0])
    max_error_theta2 = np.max(errores_locales[:, 1])
    max_error_p1 = np.max(errores_locales[:, 2])
    max_error_p2 = np.max(errores_locales[:, 3])
    
    avg_error_theta1 = np.mean(errores_locales[:, 0])
    avg_error_theta2 = np.mean(errores_locales[:, 1])
    avg_error_p1 = np.mean(errores_locales[:, 2])
    avg_error_p2 = np.mean(errores_locales[:, 3])
    
    print(f"Error máximo (theta1, theta2): ({max_error_theta1:.2e}, {max_error_theta2:.2e}) rad")
    print(f"Error máximo (p1, p2): ({max_error_p1:.2e}, {max_error_p2:.2e}) kg·m²/s")
    
    print(f"Error promedio (theta1, theta2): ({avg_error_theta1:.2e}, {avg_error_theta2:.2e}) rad")
    print(f"Error promedio (p1, p2): ({avg_error_p1:.2e}, {avg_error_p2:.2e}) kg·m²/s")
    
    print(f"Norma RMS (theta1, theta2): ({norma_error_theta1:.2e}, {norma_error_theta2:.2e}) rad")
    print(f"Norma RMS (p1, p2): ({norma_error_p1:.2e}, {norma_error_p2:.2e}) kg·m²/s")

# Análisis de sensibilidad a condiciones iniciales
if realizar_analisis_discrepancia:
    # Calcular estadísticas de discrepancia
    discrepancia_max = np.max(discrepancia)
    discrepancia_final = discrepancia[-1]
    tiempo_duplicacion = np.nan
    
    # Estimar tiempo de duplicación
    if discrepancia[0] > 0:
        idx_duplicacion = np.where(discrepancia > 2*discrepancia[0])[0]
        if len(idx_duplicacion) > 0:
            tiempo_duplicacion = Tiempo[idx_duplicacion[0]] - Tiempo[0]
    
    print("\n----- Sensibilidad a Condiciones Iniciales -----")
    print(f"Perturbación inicial: {perturbacion:.2e} rad en θ₁")
    print(f"Discrepancia máxima: {discrepancia_max:.4e}")
    print(f"Discrepancia final (t={Tiempo[-1]}s): {discrepancia_final:.4e}")
    
    if not np.isnan(tiempo_duplicacion):
        print(f"Tiempo de duplicación: {tiempo_duplicacion:.4f} s")
    
    if not np.isnan(lyapunov_exponent_pd):
        print(f"Exponente de Lyapunov: {lyapunov_exponent_pd:.4f} ± {lyapunov_error_pd:.4f}")
        estado = "caótico" if lyapunov_exponent_pd > 0.01 else "no caótico"
        print(f"Sistema detectado como: {estado}")
    
    # Análisis estadístico de diferencias entre trayectorias
    diferencia_theta1 = np.abs(Variables_norm[:, 0] - Variables_pert_norm[:, 0])
    diferencia_theta2 = np.abs(Variables_norm[:, 1] - Variables_pert_norm[:, 1])
    
    max_dif_theta1 = np.max(diferencia_theta1)
    max_dif_theta2 = np.max(diferencia_theta2)
    avg_dif_theta1 = np.mean(diferencia_theta1)
    avg_dif_theta2 = np.mean(diferencia_theta2)
    
    print("\nDiferencias entre trayectorias:")
    print(f"Máxima diferencia: (θ₁={max_dif_theta1:.2e}, θ₂={max_dif_theta2:.2e}) rad")
    print(f"Diferencia promedio: (θ₁={avg_dif_theta1:.2e}, θ₂={avg_dif_theta2:.2e}) rad")

print(f"\nDatos guardados en: {carpetaDatos}")
if generar_graficas:
    print(f"Gráficas guardadas en: {carpetaDatos}")
