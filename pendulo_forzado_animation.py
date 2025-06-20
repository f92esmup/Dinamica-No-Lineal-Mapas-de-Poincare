import numpy as np
from manim import *
import sys
import os

# Añadir el directorio actual al path para importar Funciones
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import Funciones

class SimplePenduloForzado(Scene):
    def construct(self):
        # Parámetros del péndulo forzado amortiguado
        self.coef_amortiguamiento = 0.2     # Coeficiente de amortiguamiento (1/s)
        self.frecuencia_natural = 1.0       # Frecuencia natural (rad/s)
        self.amplitud = 0.5                 # Amplitud de la fuerza externa (1/s²)
        self.frecuencia_externa = 2.0       # Frecuencia de la fuerza externa (rad/s)
        
        # Parámetros de simulación
        self.t_inicial = 0.0
        self.t_final = 30.0      # 30 segundos de simulación para ver mejor la dinámica
        self.dt = 0.02           # Paso de tiempo pequeño para capturar la dinámica
        
        # Escala para la visualización
        self.scale = 2.5
        self.L = 1.0  # Longitud visual del péndulo
        
        # Condiciones iniciales [theta, dtheta/dt]
        condiciones_iniciales = np.array([0.1, 0.0])  # Pequeña perturbación inicial
        
        print(f"Péndulo forzado amortiguado con Amplitud = {self.amplitud}")
        print(f"Condiciones iniciales: θ={condiciones_iniciales[0]:.3f} rad, dθ/dt={condiciones_iniciales[1]:.3f} rad/s")
        print(f"Parámetros: γ={self.coef_amortiguamiento}, ω₀={self.frecuencia_natural}, F={self.amplitud}, Ω={self.frecuencia_externa}")
        
        # Integrar las ecuaciones de movimiento
        print("Integrando ecuaciones de movimiento...")
        tiempo, variables, _ = Funciones.RungeKutta(
            self.pendulo_forzado_equations,
            condiciones_iniciales,
            self.t_inicial,
            self.t_final,
            self.dt
        )
        
        # Convertir ángulos a posiciones cartesianas
        x_data, y_data = self.angle_to_position(variables[:, 0])
        
        # Crear elementos visuales
        # Punto de suspensión (fijo)
        pivot = Dot(ORIGIN, color=WHITE, radius=0.1)
        
        # Masa del péndulo
        mass = Dot(color=YELLOW, radius=0.2)
        
        # Varilla del péndulo
        rod = Line(ORIGIN, [x_data[0] * self.scale, -y_data[0] * self.scale, 0], 
                  color=WHITE, stroke_width=4)
        
        # Posicionar la masa inicial
        mass.move_to([x_data[0] * self.scale, -y_data[0] * self.scale, 0])
        
        # Añadir elementos a la escena
        self.add(pivot, rod, mass)
        
        # Función para actualizar el péndulo
        def update_pendulum(mob, alpha):
            # Calcular índice de tiempo basado en alpha
            idx = int(alpha * (len(tiempo) - 1))
            if idx >= len(tiempo):
                idx = len(tiempo) - 1
                
            # Posiciones actuales
            x = x_data[idx] * self.scale
            y = -y_data[idx] * self.scale  # Negativo porque Y aumenta hacia arriba en Manim
            
            # Actualizar varilla
            rod.become(Line(ORIGIN, [x, y, 0], color=WHITE, stroke_width=4))
            
            # Actualizar masa
            mass.move_to([x, y, 0])
            
            return mob
        
        # Crear un grupo que contenga todos los elementos móviles
        pendulum_group = Group(rod, mass)
        
        # Animar el péndulo durante el tiempo especificado
        self.play(
            UpdateFromAlphaFunc(pendulum_group, update_pendulum),
            run_time=15,  # 15 segundos de animación
            rate_func=linear
        )
        
        # Pausa al final
        self.wait(1)
    
    def pendulo_forzado_equations(self, t, variables):
        """
        Ecuaciones de movimiento del péndulo forzado amortiguado.
        
        Variables: [theta, dtheta/dt]
        
        La ecuación diferencial es:
        d²θ/dt² = -γ(dθ/dt) - (ω₀² + 2F cos(Ωt)) sin(θ)
        
        donde:
        γ = coeficiente de amortiguamiento
        ω₀ = frecuencia natural
        F = amplitud de la fuerza externa
        Ω = frecuencia de la fuerza externa
        """
        theta, dtheta_dt = variables
        
        # Primera ecuación: dθ/dt
        du1 = dtheta_dt
        
        # Segunda ecuación: d²θ/dt²
        du2 = (-self.coef_amortiguamiento * dtheta_dt - 
               (self.frecuencia_natural**2 + 2 * self.amplitud * np.cos(self.frecuencia_externa * t)) * 
               np.sin(theta))
        
        return np.array([du1, du2])
    
    def angle_to_position(self, theta):
        """Convierte ángulo a posición cartesiana"""
        # Posición de la masa
        x = self.L * np.sin(theta)
        y = self.L * np.cos(theta)
        
        return x, y


class PenduloForzadoConTraza(Scene):
    def construct(self):
        # Parámetros del péndulo forzado amortiguado
        self.coef_amortiguamiento = 0.2     
        self.frecuencia_natural = 1.0       
        self.amplitud = 0.5                 
        self.frecuencia_externa = 2.0       
        
        # Parámetros de simulación
        self.t_inicial = 0.0
        self.t_final = 30.0      # Más tiempo para ver la evolución
        self.dt = 0.02           
        
        # Escala para la visualización
        self.scale = 2.5
        self.L = 1.0  
        
        # Condiciones iniciales 
        condiciones_iniciales = np.array([0.1, 0.0])  
        
        # Integrar las ecuaciones de movimiento
        tiempo, variables, _ = Funciones.RungeKutta(
            self.pendulo_forzado_equations,
            condiciones_iniciales,
            self.t_inicial,
            self.t_final,
            self.dt
        )
        
        # Convertir ángulos a posiciones cartesianas
        x_data, y_data = self.angle_to_position(variables[:, 0])
        
        # Crear elementos visuales
        pivot = Dot(ORIGIN, color=WHITE, radius=0.1)
        mass = Dot(color=YELLOW, radius=0.15)
        rod = Line(ORIGIN, [x_data[0] * self.scale, -y_data[0] * self.scale, 0], 
                  color=WHITE, stroke_width=3)
        
        # Crear traza de la masa
        trace = VMobject(color=ORANGE, stroke_width=2, stroke_opacity=0.7)
        trace.set_points_as_corners([])
        
        # Posicionar la masa inicial
        mass.move_to([x_data[0] * self.scale, -y_data[0] * self.scale, 0])
        
        # Añadir elementos a la escena
        self.add(pivot, rod, mass, trace)
        
        # Variables para la traza
        trace_points = []
        max_trace_points = 100  # Limitar puntos para rendimiento
        
        # Función para actualizar el péndulo con traza
        def update_pendulum_with_trace(mob, alpha):
            idx = int(alpha * (len(tiempo) - 1))
            if idx >= len(tiempo):
                idx = len(tiempo) - 1
                
            x = x_data[idx] * self.scale
            y = -y_data[idx] * self.scale
            
            # Actualizar varilla y masa
            rod.become(Line(ORIGIN, [x, y, 0], color=WHITE, stroke_width=3))
            mass.move_to([x, y, 0])
            
            # Actualizar traza
            trace_points.append([x, y, 0])
            if len(trace_points) > max_trace_points:
                trace_points.pop(0)
            
            if len(trace_points) > 1:
                trace.set_points_as_corners(trace_points)
            
            return mob
        
        pendulum_group = Group(rod, mass)
        
        # Animar el péndulo con traza
        self.play(
            UpdateFromAlphaFunc(pendulum_group, update_pendulum_with_trace),
            run_time=15,  
            rate_func=linear
        )
        
        self.wait(2)
    
    def pendulo_forzado_equations(self, t, variables):
        """Ecuaciones de movimiento del péndulo forzado amortiguado"""
        theta, dtheta_dt = variables
        
        du1 = dtheta_dt
        du2 = (-self.coef_amortiguamiento * dtheta_dt - 
               (self.frecuencia_natural**2 + 2 * self.amplitud * np.cos(self.frecuencia_externa * t)) * 
               np.sin(theta))
        
        return np.array([du1, du2])
    
    def angle_to_position(self, theta):
        """Convierte ángulo a posición cartesiana"""
        x = self.L * np.sin(theta)
        y = self.L * np.cos(theta)
        return x, y
