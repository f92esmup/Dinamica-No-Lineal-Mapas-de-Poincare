#!/usr/bin/env python3
"""
Manim animation showing sensitivity to initial conditions using the Hénon Map.
Two particles start at nearly identical positions and diverge dramatically over time.
"""

from manim import *
import numpy as np

class HenonMapDivergence(Scene):
    def construct(self):
        # 1. Define the Hénon map function
        def henon_map(x, y, a=1.4, b=0.3):
            """
            Hénon map function: transforms (x, y) to the next point
            Standard parameters: a=1.4, b=0.3
            """
            x_next = 1 - a * x**2 + y
            y_next = b * x
            return x_next, y_next

        # 2. Set up the scene
        # Create 2D axes covering the typical range of the Hénon attractor
        axes = Axes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-0.5, 0.5, 0.25],
            x_length=10,
            y_length=6,
            axis_config={"color": GRAY}
        )
        
        # Add labels to axes using Text instead of LaTeX
        x_label = Text("x", font_size=24).next_to(axes.get_x_axis().get_end(), RIGHT)
        y_label = Text("y", font_size=24).next_to(axes.get_y_axis().get_end(), UP)
        
        # Create title
        title = Text("Hénon Map: Sensitivity to Initial Conditions", font_size=36)
        title.to_edge(UP)
        
        # Add axes and title to scene
        self.add(axes, x_label, y_label, title)
        
        # 3. Create the initial points
        # Two nearly identical initial positions
        p1_start = np.array([0.1, 0.1])
        p2_start = np.array([0.1, 0.10001])  # Tiny difference: 0.00001
        
        # Convert to scene coordinates
        p1_pos = axes.coords_to_point(p1_start[0], p1_start[1])
        p2_pos = axes.coords_to_point(p2_start[0], p2_start[1])
        
        # Create dots
        dot1 = Dot(p1_pos, color=BLUE, radius=0.05)
        dot2 = Dot(p2_pos, color=RED, radius=0.05)
        
        # Add labels for the dots
        label1 = Text("Particle 1", font_size=20, color=BLUE).next_to(dot1, UP)
        label2 = Text("Particle 2", font_size=20, color=RED).next_to(dot2, DOWN)
        
        # 4. Create TracedPath objects for both dots
        # This will leave a trail behind the dots to draw the attractor
        trace1 = TracedPath(dot1.get_center, stroke_color=BLUE, stroke_width=2, stroke_opacity=0.7)
        trace2 = TracedPath(dot2.get_center, stroke_color=RED, stroke_width=2, stroke_opacity=0.7)
        
        # Add everything to the scene
        self.add(dot1, dot2, label1, label2, trace1, trace2)
        
        # Show initial setup
        self.wait(1)
        
        # Add initial condition text
        initial_text = Text(
            "Initial positions differ by only 0.00001",
            font_size=24,
            color=YELLOW
        ).to_edge(DOWN)
        self.play(Write(initial_text))
        self.wait(1)
        
        # 5. Animate the iterations
        current_p1 = p1_start.copy()
        current_p2 = p2_start.copy()
        
        # Number of iterations
        n_iterations = 150
        
        # Create a simple counter using Text instead of Variable to avoid LaTeX
        iteration_counter = Text("Iteration: 0", font_size=24)
        iteration_counter.to_corner(UL)
        self.add(iteration_counter)
        
        for i in range(n_iterations):
            # Calculate next positions using Hénon map
            next_p1 = henon_map(current_p1[0], current_p1[1])
            next_p2 = henon_map(current_p2[0], current_p2[1])
            
            # Convert to scene coordinates
            next_pos1 = axes.coords_to_point(next_p1[0], next_p1[1])
            next_pos2 = axes.coords_to_point(next_p2[0], next_p2[1])
            
            # Create animations for smooth movement
            new_counter_text = Text(f"Iteration: {i + 1}", font_size=24).to_corner(UL)
            animations = [
                dot1.animate.move_to(next_pos1),
                dot2.animate.move_to(next_pos2),
                Transform(iteration_counter, new_counter_text)
            ]
            
            # Adjust animation speed - faster in the beginning, slower later to show divergence
            if i < 20:
                run_time = 0.1  # Fast initial iterations
            elif i < 50:
                run_time = 0.08
            else:
                run_time = 0.05  # Slower to appreciate the chaotic behavior
            
            # Play the animation
            self.play(*animations, run_time=run_time)
            
            # Update current positions
            current_p1 = np.array(next_p1)
            current_p2 = np.array(next_p2)
            
            # Update labels to follow dots (only for first few iterations)
            if i < 10:
                label1.next_to(dot1, UP)
                label2.next_to(dot2, DOWN)
            elif i == 10:
                # Fade out labels to avoid clutter
                self.play(FadeOut(label1), FadeOut(label2), run_time=0.3)
            
            # Add divergence text after some iterations
            if i == 30:
                self.play(FadeOut(initial_text))
                divergence_text = Text(
                    "The particles are now following completely different paths!",
                    font_size=24,
                    color=YELLOW
                ).to_edge(DOWN)
                self.play(Write(divergence_text))
            elif i == 60:
                self.play(FadeOut(divergence_text))
                chaos_text = Text(
                    "This is the essence of chaos: sensitive dependence on initial conditions",
                    font_size=22,
                    color=YELLOW
                ).to_edge(DOWN)
                self.play(Write(chaos_text))
        
        # Final pause to appreciate the full attractor
        self.wait(3)
        
        # Add final explanatory text
        final_text = VGroup(
            Text("The Hénon Map creates a strange attractor", font_size=24),
            Text("Tiny differences lead to completely different trajectories", font_size=20),
            Text("This demonstrates the butterfly effect in dynamical systems", font_size=20)
        ).arrange(DOWN, buff=0.3).to_edge(DOWN)
        
        self.play(FadeOut(chaos_text))
        self.play(Write(final_text))
        self.wait(4)


if __name__ == "__main__":
    # This allows running the script directly
    scene = HenonMapDivergence()
    scene.render()
