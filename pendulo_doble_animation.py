import numpy as np
from manim import *
import sys
import os

# Añadir el directorio actual al path para importar Funciones
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import Funciones

class SimplePenduloDoble(Scene):
    def construct(self):
        # Parámetros del péndulo doble
        self.m = 1.0     # Masa de cada partícula (kg)
        self.L = 1.0     # Longitud de cada varilla (m)
        self.g = 9.8     # Gravedad (m/s²)
        
        # Parámetros de simulación
        self.t_inicial = 0.0
        self.t_final = 10.0      # 10 segundos de simulación (reducido)
        self.dt = 0.05           # Paso de tiempo mayor para eficiencia
        
        # Escala para la visualización (convertir metros a unidades de Manim)
        self.scale = 2.0
        
        # Configurar energía de 15J
        energia_deseada = 15.0
        condiciones_iniciales = np.array([1.1, 0.0, 0.0, 0.0])
        
        # Calcular p2 para obtener la energía deseada
        condiciones_iniciales[3] = self.ecuacion_segundo_grado(
            energia_deseada, condiciones_iniciales[:3], '+'
        )
        
        print(f"Condiciones iniciales para E={energia_deseada}J:")
        print(f"θ₁={condiciones_iniciales[0]:.3f} rad, θ₂={condiciones_iniciales[1]:.3f} rad")
        print(f"p₁={condiciones_iniciales[2]:.3f} kg·m²/s, p₂={condiciones_iniciales[3]:.3f} kg·m²/s")
        
        # Integrar las ecuaciones de movimiento
        print("Integrando ecuaciones de movimiento...")
        tiempo, variables, _ = Funciones.RungeKutta(
            self.pendulo_doble_equations,
            condiciones_iniciales,
            self.t_inicial,
            self.t_final,
            self.dt
        )
        
        # Convertir ángulos a posiciones cartesianas
        x1_data, y1_data, x2_data, y2_data = self.angles_to_positions(variables)
        
        # Crear elementos visuales
        # Punto de suspensión (fijo)
        pivot = Dot(ORIGIN, color=WHITE, radius=0.1)
        
        # Masas del péndulo
        mass1 = Dot(color=RED, radius=0.15)
        mass2 = Dot(color=BLUE, radius=0.15)
        
        # Varillas del péndulo
        rod1 = Line(ORIGIN, [x1_data[0] * self.scale, -y1_data[0] * self.scale, 0], color=WHITE, stroke_width=3)
        rod2 = Line([x1_data[0] * self.scale, -y1_data[0] * self.scale, 0], 
                   [x2_data[0] * self.scale, -y2_data[0] * self.scale, 0], color=WHITE, stroke_width=3)
        
        # Trazas de las masas (opcional - comentado para simplicidad)
        # trace1 = VMobject(color=RED, stroke_width=1, stroke_opacity=0.5)
        # trace2 = VMobject(color=BLUE, stroke_width=1, stroke_opacity=0.5)
        
        # Posicionar las masas iniciales
        mass1.move_to([x1_data[0] * self.scale, -y1_data[0] * self.scale, 0])
        mass2.move_to([x2_data[0] * self.scale, -y2_data[0] * self.scale, 0])
        
        # Añadir elementos a la escena
        self.add(pivot, rod1, rod2, mass1, mass2)
        
        # Función para actualizar el péndulo
        def update_pendulum(mob, alpha):
            # Calcular índice de tiempo basado en alpha
            idx = int(alpha * (len(tiempo) - 1))
            if idx >= len(tiempo):
                idx = len(tiempo) - 1
                
            # Posiciones actuales
            x1 = x1_data[idx] * self.scale
            y1 = -y1_data[idx] * self.scale  # Negativo porque Y aumenta hacia arriba en Manim
            x2 = x2_data[idx] * self.scale
            y2 = -y2_data[idx] * self.scale
            
            # Actualizar varillas
            rod1.become(Line(ORIGIN, [x1, y1, 0], color=WHITE, stroke_width=3))
            rod2.become(Line([x1, y1, 0], [x2, y2, 0], color=WHITE, stroke_width=3))
            
            # Actualizar masas
            mass1.move_to([x1, y1, 0])
            mass2.move_to([x2, y2, 0])
            
            return mob
        
        # Crear un grupo que contenga todos los elementos móviles
        pendulum_group = Group(rod1, rod2, mass1, mass2)
        
        # Animar el péndulo durante el tiempo especificado
        self.play(
            UpdateFromAlphaFunc(pendulum_group, update_pendulum),
            run_time=8,  # 8 segundos de animación
            rate_func=linear
        )
        
        # Pausa al final
        self.wait(1)
    
    def pendulo_doble_equations(self, t, variables):
        """Ecuaciones de Hamilton del péndulo doble"""
        theta1, theta2, p1, p2 = variables
        
        # Calcular diferencia de ángulos
        delta = theta1 - theta2
        sin_delta = np.sin(delta)
        cos_delta = np.cos(delta)
        
        # Denominador común
        D_denom = 1 + sin_delta**2
        
        # Derivadas
        dtheta1_dt = (2*p1 - 2*p2*cos_delta) / (2*self.L**2*self.m*D_denom)
        dtheta2_dt = (-2*p1*cos_delta + 4*p2) / (2*self.L**2*self.m*D_denom)
        
        # Términos para dp1/dt
        term1 = -2*self.L*self.g*self.m*np.sin(theta1)
        term2 = -(p1*p2*sin_delta) / (self.L**2*self.m*D_denom)
        term3 = ((p1**2 - 2*p1*p2*cos_delta + 2*p2**2) * sin_delta * cos_delta) / (self.L**2*self.m*D_denom**2)
        dp1_dt = term1 + term2 + term3
        
        # Términos para dp2/dt
        term1 = -self.L*self.g*self.m*np.sin(theta2)
        term2 = (p1*p2*sin_delta) / (self.L**2*self.m*D_denom)
        term3 = -((p1**2 - 2*p1*p2*cos_delta + 2*p2**2) * sin_delta * cos_delta) / (self.L**2*self.m*D_denom**2)
        dp2_dt = term1 + term2 + term3
        
        return np.array([dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt])
    
    def ecuacion_segundo_grado(self, energia_deseada, condiciones_parciales, signo='+'):
        """Calcula p2 para una energía dada"""
        theta1, theta2, p1 = condiciones_parciales
        
        delta = theta1 - theta2
        sin_delta_sq = np.sin(delta)**2
        cos_delta = np.cos(delta)
        
        D_denom = 1 + sin_delta_sq
        
        # Energía potencial
        U = self.m*self.g*self.L*(3 - 2*np.cos(theta1) - np.cos(theta2))
        
        # Energía cinética objetivo
        T_objetivo = energia_deseada - U
        
        # Coeficientes de la ecuación cuadrática en p2
        a = 2 / (2*self.m*self.L**2 * D_denom)
        b = -2*p1*cos_delta / (2*self.m*self.L**2 * D_denom)
        c = p1**2 / (2*self.m*self.L**2 * D_denom) - T_objetivo
        
        # Resolver ecuación cuadrática
        discriminante = b**2 - 4*a*c
        
        if discriminante < 0:
            return np.nan
        
        if signo == '+':
            return (-b + np.sqrt(discriminante)) / (2*a)
        else:
            return (-b - np.sqrt(discriminante)) / (2*a)
    
    def angles_to_positions(self, variables):
        """Convierte ángulos a posiciones cartesianas"""
        theta1 = variables[:, 0]
        theta2 = variables[:, 1]
        
        # Posición de la primera masa
        x1 = self.L * np.sin(theta1)
        y1 = self.L * np.cos(theta1)
        
        # Posición de la segunda masa
        x2 = x1 + self.L * np.sin(theta2)
        y2 = y1 + self.L * np.cos(theta2)
        
        return x1, y1, x2, y2
