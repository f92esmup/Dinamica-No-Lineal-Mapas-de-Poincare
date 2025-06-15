from manim import *
import numpy as np
from scipy.integrate import odeint

class LorenzAttractor(Scene):
    def construct(self):
        # Parámetros del atractor de Lorenz
        sigma = 10.0
        rho = 28.0
        beta = 8.0/3.0
        
        # Sistema de ecuaciones diferenciales de Lorenz
        def lorenz_system(state, t):
            x, y, z = state
            dx_dt = sigma * (y - x)
            dy_dt = x * (rho - z) - y
            dz_dt = x * y - beta * z
            return [dx_dt, dy_dt, dz_dt]
        
        # Condiciones iniciales para dos partículas
        initial_state_1 = [1.0, 1.0, 1.0]
        initial_state_2 = [1.01, 1.0, 1.0]  # Ligeramente diferente para mostrar divergencia
        
        # Tiempo de integración
        t = np.linspace(0, 20, 2000)
        
        # Resolver las ecuaciones para ambas partículas
        trajectory_1 = odeint(lorenz_system, initial_state_1, t)
        trajectory_2 = odeint(lorenz_system, initial_state_2, t)
        
        # Calcular el centro del atractor para centrarlo en pantalla
        all_x = np.concatenate([trajectory_1[:, 0], trajectory_2[:, 0]])
        all_z = np.concatenate([trajectory_1[:, 2], trajectory_2[:, 2]])
        center_x = (np.max(all_x) + np.min(all_x)) / 2
        center_z = (np.max(all_z) + np.min(all_z)) / 2
        
        # Escalar las trayectorias para que se vean bien en pantalla
        scale = 0.15
        
        # Crear las trayectorias como funciones paramétricas
        def trajectory_func_1(t_param):
            idx = int(t_param * (len(trajectory_1) - 1))
            if idx >= len(trajectory_1):
                idx = len(trajectory_1) - 1
            x, y, z = trajectory_1[idx]
            return np.array([(x - center_x) * scale, (z - center_z) * scale, 0])  # Proyección x-z centrada
        
        def trajectory_func_2(t_param):
            idx = int(t_param * (len(trajectory_2) - 1))
            if idx >= len(trajectory_2):
                idx = len(trajectory_2) - 1
            x, y, z = trajectory_2[idx]
            return np.array([(x - center_x) * scale, (z - center_z) * scale, 0])  # Proyección x-z centrada
        
        # Crear las trayectorias
        trajectory_curve_1 = ParametricFunction(
            trajectory_func_1,
            t_range=[0, 1],
            color=BLUE,
            stroke_width=2
        )
        
        trajectory_curve_2 = ParametricFunction(
            trajectory_func_2,
            t_range=[0, 1],
            color=RED,
            stroke_width=2
        )
        
        # Crear las partículas (puntos móviles)
        particle_1 = Dot(radius=0.08, color=BLUE)
        particle_2 = Dot(radius=0.08, color=RED)
        
        # Crear las trayectorias que se irán dibujando dinámicamente
        trail_1 = VMobject(color=BLUE, stroke_width=2, stroke_opacity=0.7)
        trail_2 = VMobject(color=RED, stroke_width=2, stroke_opacity=0.7)
        
        # Posiciones iniciales
        particle_1.move_to(trajectory_func_1(0))
        particle_2.move_to(trajectory_func_2(0))
        
        # Agregar las partículas y las trayectorias vacías
        self.add(particle_1, particle_2, trail_1, trail_2)
        
        # Listas para almacenar los puntos de las trayectorias
        trail_points_1 = [trajectory_func_1(0)]
        trail_points_2 = [trajectory_func_2(0)]
        
        # Función para actualizar las partículas
        def update_particle_1(mob, alpha):
            pos = trajectory_func_1(alpha)
            mob.move_to(pos)
            # Agregar punto a la trayectoria
            trail_points_1.append(pos)
            if len(trail_points_1) > 1:
                trail_1.clear_points()
                trail_1.set_points_smoothly(trail_points_1)  # Usar suavizado
            return mob
            
        def update_particle_2(mob, alpha):
            pos = trajectory_func_2(alpha)
            mob.move_to(pos)
            # Agregar punto a la trayectoria
            trail_points_2.append(pos)
            if len(trail_points_2) > 1:
                trail_2.clear_points()
                trail_2.set_points_smoothly(trail_points_2)  # Usar suavizado
            return mob
        
        # Animar las partículas y sus trayectorias
        self.play(
            UpdateFromAlphaFunc(particle_1, update_particle_1),
            UpdateFromAlphaFunc(particle_2, update_particle_2),
            run_time=15,
            rate_func=linear
        )
        
        self.wait(2)
