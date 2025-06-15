from manim import *
import numpy as np

class SimplePoincareMapExplanation(Scene):
    def construct(self):
        # 1. Setup the Scene
        # Create coordinate axes without grid
        x_axis = Line(start=LEFT*4, end=RIGHT*4, color=WHITE, stroke_width=2)
        y_axis = Line(start=DOWN*3, end=UP*3, color=WHITE, stroke_width=2)
        self.add(x_axis, y_axis)

        # 2. Define the objects
        # Define a logarithmic spiral using a ParametricFunction. It should start far from the center.
        # This will be our system's "trajectory".
        trajectory = ParametricFunction(
            lambda t: np.array([2.8 * np.exp(-0.08*t) * np.cos(t), 2.8 * np.exp(-0.08*t) * np.sin(t), 0]),
            t_range=np.array([0, 25]),
            color=BLUE,
            stroke_width=4
        )
        
        # Define a diagonal line as our "Poincaré Section" (bisectriz entre 2° y 4° cuadrante)
        # La línea y = -x pasa por la bisectriz del 2° y 4° cuadrante
        poincare_section = Line(UP*3 + LEFT*3, DOWN*3 + RIGHT*3, color=RED, stroke_width=8)

        # 3. Animate the process step-by-step
        # First show the Poincaré section
        self.play(Create(poincare_section), run_time=2)
        self.wait(1)
        
        # Then animate the creation of the trajectory
        self.play(Create(trajectory), run_time=4)
        self.wait(2)
        
        # Calculate intersection points with the diagonal line y = -x
        # For the spiral: x = 2.8*exp(-0.08*t)*cos(t), y = 2.8*exp(-0.08*t)*sin(t)
        # Intersection when y = -x, so: 2.8*exp(-0.08*t)*sin(t) = -2.8*exp(-0.08*t)*cos(t)
        # This gives: tan(t) = -1, so t = -π/4 + nπ
        
        intersection_points = []
        for n in range(8):  # Más puntos para mejor visualización
            t_intersect = -np.pi/4 + n * np.pi
            if 0 <= t_intersect <= 25:
                x_pos = 2.8 * np.exp(-0.08 * t_intersect) * np.cos(t_intersect)
                y_pos = 2.8 * np.exp(-0.08 * t_intersect) * np.sin(t_intersect)
                # Verificar que realmente esté cerca de la línea y = -x
                if abs(y_pos + x_pos) < 0.1:  # Tolerancia pequeña
                    intersection_points.append(np.array([x_pos, y_pos, 0]))
        
        # Si no hay suficientes puntos calculados, usar puntos aproximados
        if len(intersection_points) < 4:
            intersection_points = [
                np.array([1.9, -1.9, 0]),
                np.array([0.7, -0.7, 0]),
                np.array([-0.2, 0.2, 0]),
                np.array([-0.9, 0.9, 0]),
                np.array([-0.35, 0.35, 0]),
                np.array([-0.12, 0.12, 0])
            ]

        # Loop through the intersection points and show the "recording" mechanism.
        poincare_dots = VGroup()
        
        for i, point in enumerate(intersection_points[:6]):  # Tomar solo los primeros 6
            # Create a dot at the intersection.
            dot = Dot(point, color=YELLOW, radius=0.12)
            
            # Add a small arrow showing the direction of crossing
            # Dirección perpendicular a la línea diagonal
            arrow_direction = np.array([1, 1, 0]) / np.sqrt(2)  # Perpendicular a y=-x
            arrow = Arrow(
                start=point + arrow_direction*0.4,
                end=point - arrow_direction*0.4,
                color=GREEN,
                stroke_width=4,
                max_tip_length_to_length_ratio=0.3
            )
            
            # Animate a Flash and the creation of the dot.
            self.play(
                Flash(dot, color=YELLOW, flash_radius=0.4),
                Create(arrow),
                run_time=0.8
            )
            self.play(FadeIn(dot), FadeOut(arrow), run_time=0.5)
            poincare_dots.add(dot)
            
            # Add label for each point
            point_label = Text(f"P{i+1}", font_size=20, color=YELLOW).next_to(dot, UP*0.5 + RIGHT*0.5, buff=0.1)
            self.play(Write(point_label), run_time=0.3)
            poincare_dots.add(point_label)

        self.wait(2)

        # 4. Show the final result
        # Fade out the trajectory and the Poincaré section
        self.play(FadeOut(trajectory), FadeOut(poincare_section), run_time=2)
        
        # Emphasize the final points.
        dots_only = VGroup(*[mob for mob in poincare_dots if isinstance(mob, Dot)])
        self.play(dots_only.animate.set_color(RED).scale(1.3), run_time=1.5)
        
        self.wait(3)
