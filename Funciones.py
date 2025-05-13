import numpy as np
"""
Definición de las funciones comunes utilizadas en las simulaciones de sistemas dinámicos.
Incluye el método Runge-Kutta de 4º orden y funciones de normalización de ángulos.

Sistema de unidades consistente:
- Tiempo: segundos (s)
- Longitud: metros (m)
- Masa: kilogramos (kg)
- Ángulos: radianes (rad)
- Momento angular: kg·m²/s
- Energía: Julios (J)
"""

#Función que implementa el método Runge-Kutta:
def RungeKutta(Funcion, CondicionesIniciales, TInicial, TFinal, TamañoPaso):
    """
    Entradas:

    Función -> Función de primer orden a la que se le aplicará el método.
    CondicionesIniciales -> Incluye la posición y velocidad inicial de nuestro sistema.
    TInicial -> Punto inicial, donde empezamos a evaluar la función (s)
    TFinal -> Punto final, hasta donde evaluamos la función (s)
    TamañoPaso -> Es el incremento de tiempo, la 'cantidad' de tiempo que avanzamos en cada paso (s).
                 A partir de este valor y de los tiempos inicial y final es que vamos a obterner 
                 el número de pasos

    Salidas: 
    Tiempo -> Vector fila que agrupa todos los tiempos en los que se ha evaluado la función (s), 
              va desde TInicial hata TFinal (equiespaciados) con una dimensión igual al numero de pasos.
    Variables -> Array con las variables del sistema en cada paso de tiempo.
    Derivadas -> Array con las derivadas de las variables del sistema en cada paso de tiempo.
    """

    #Definimos el número de pasos:
    NumeroPasos = int((TFinal - TInicial) / TamañoPaso) + 1  # +1 para incluir exactamente el punto final

    #Incializamos las variables: 
    Tiempo = np.linspace(TInicial, TFinal, NumeroPasos)  # Usar linspace directamente para garantizar TFinal exacto
    Variables = np.zeros((NumeroPasos, len(CondicionesIniciales)))
    Derivadas = np.zeros((NumeroPasos, len(CondicionesIniciales))) # Array para guardar derivadas

    #Aplico las condiciones iniciales:
    Variables[0] = CondicionesIniciales
    Derivadas[0] = np.array(Funcion(Tiempo[0], Variables[0])) # Calcular derivada inicial

    #Aplico el método Runge-Kutta:
    for i in range(NumeroPasos - 1): # Iterar hasta NumeroPasos - 1
        K1 = np.array(Funcion(Tiempo[i], Variables[i]))
        K2 = np.array(Funcion(Tiempo[i] + TamañoPaso / 2, Variables[i] + (TamañoPaso / 2) * K1))
        K3 = np.array(Funcion(Tiempo[i] + TamañoPaso / 2, Variables[i] + (TamañoPaso / 2) * K2))
        K4 = np.array(Funcion(Tiempo[i] + TamañoPaso, Variables[i] + TamañoPaso * K3))
        
        Variables[i+1,:] = Variables[i,:] + (TamañoPaso / 6) * (K1 + 2 * K2 + 2 * K3 + K4)     
        Derivadas[i+1,:] = np.array(Funcion(Tiempo[i+1], Variables[i+1,:])) # Calcular derivada en el nuevo punto
        
        # Eliminado el bloque de normalización de ángulos para preservar la precisión del integrador RK4
        # y la conservación de energía del sistema Hamiltoniano

    return Tiempo, Variables, Derivadas # Devolver Derivadas

#Función para normalizar ángulos:
def NormalizarAngulos(Variables):
    """
    Normaliza todos los ángulos para que estén en el rango [-pi, pi].
    Esta función se aplica DESPUÉS de terminar la integración numérica,
    no durante el proceso de integración.
    
    Parámetros:
    -----------
    Variables: ndarray
        Array con las variables del sistema. Se asume que las dos primeras
        columnas son ángulos (rad) que deben ser normalizados.
        
    Retorna:
    --------
    Variables_norm: ndarray
        Array con las mismas dimensiones que el original pero con los ángulos
        normalizados al rango [-pi, pi] (rad).
    """
    Variables_norm = Variables.copy()
    
    # Asumimos que las dos primeras columnas son ángulos
    # (como en el péndulo simple, forzado o doble)
    n_cols = min(2, Variables.shape[1])
    
    for j in range(n_cols):
        Variables_norm[:, j] = (Variables[:, j] + np.pi) % (2 * np.pi) - np.pi
    
    return Variables_norm