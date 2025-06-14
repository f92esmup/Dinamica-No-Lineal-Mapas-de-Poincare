from manim import *

class SpiralAndColorClock(Scene):
    def construct(self):
        # 1. SETUP THE SCENE
        # Create simple coordinate axes (just lines, no arrows or markers)
        x_axis = Line(start=[-3, 0, 0], end=[3, 0, 0], color=WHITE, stroke_width=2)
        y_axis = Line(start=[0, -3, 0], end=[0, 3, 0], color=WHITE, stroke_width=2)
        
        trajectory = ParametricFunction(
            lambda t: np.array([
                2.5 * np.exp(-t/8) * np.cos(2*t),
                2.5 * np.exp(-t/8) * np.sin(2*t),
                0
            ]), t_range=[0, 9], color=BLUE
        )
        
        # Create progress bar at the bottom
        progress_bar_border = Rectangle(width=8, height=0.5, color=WHITE, stroke_width=2).to_edge(DOWN, buff=0.5)
        progress_bar_fill = Rectangle(width=0.1, height=0.4, color=YELLOW, fill_opacity=0.8)
        progress_bar_fill.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
        
        # Create time label
        time_label = Text("0T", font_size=36).next_to(progress_bar_border, LEFT, buff=0.3)
        
        self.add(x_axis, y_axis, progress_bar_border, time_label)

        # 2. ANIMATE THE CONTINUOUS TRAJECTORY
        # Animate the creation of the full blue spiral over 1.5 seconds (much faster).
        # This shows the entire path the system will take.
        self.play(Create(trajectory), run_time=1.5, rate_func=linear)
        self.wait(0.5)

        # 3. ANIMATE THE SAMPLING PROCESS WITH PROGRESS BAR
        # Now, we will highlight the sampling using the progress bar
        
        # Period 1 (0s to 3s)
        # Fill the progress bar over 3 seconds, then mark Poincaré point
        self.add(progress_bar_fill)
        
        def update_progress_bar_1(mob, alpha):
            new_width = alpha * 7.9  # 7.9 to leave small margin
            mob.become(Rectangle(width=new_width, height=0.4, color=YELLOW, fill_opacity=0.8))
            mob.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
            return mob
            
        self.play(
            UpdateFromAlphaFunc(progress_bar_fill, update_progress_bar_1),
            run_time=3,
            rate_func=linear
        )
        
        dot1 = Dot(trajectory.point_from_proportion(3/9), color=RED, radius=0.08)
        self.play(Flash(dot1, flash_radius=0.3), FadeIn(dot1))
        
        # Update time label to 1T
        self.play(Transform(time_label, Text("1T", font_size=36).next_to(progress_bar_border, LEFT, buff=0.3)))
        
        # Reset progress bar for period 2
        progress_bar_fill.become(Rectangle(width=0.1, height=0.4, color=ORANGE, fill_opacity=0.8))
        progress_bar_fill.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
        self.wait(0.2)
        
        # Period 2 (3s to 6s)
        def update_progress_bar_2(mob, alpha):
            new_width = alpha * 7.9
            mob.become(Rectangle(width=new_width, height=0.4, color=ORANGE, fill_opacity=0.8))
            mob.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
            return mob
            
        self.play(
            UpdateFromAlphaFunc(progress_bar_fill, update_progress_bar_2),
            run_time=3,
            rate_func=linear
        )

        dot2 = Dot(trajectory.point_from_proportion(6/9), color=RED, radius=0.08)
        self.play(Flash(dot2, flash_radius=0.3), FadeIn(dot2))
        
        # Update time label to 2T
        self.play(Transform(time_label, Text("2T", font_size=36).next_to(progress_bar_border, LEFT, buff=0.3)))
        
        # Reset progress bar for period 3
        progress_bar_fill.become(Rectangle(width=0.1, height=0.4, color=GREEN, fill_opacity=0.8))
        progress_bar_fill.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
        self.wait(0.2)

        # Period 3 (6s to 9s)
        def update_progress_bar_3(mob, alpha):
            new_width = alpha * 7.9
            mob.become(Rectangle(width=new_width, height=0.4, color=GREEN, fill_opacity=0.8))
            mob.move_to(progress_bar_border.get_center()).align_to(progress_bar_border, LEFT).shift(RIGHT * 0.05)
            return mob
            
        self.play(
            UpdateFromAlphaFunc(progress_bar_fill, update_progress_bar_3),
            run_time=3,
            rate_func=linear
        )
        
        dot3 = Dot(trajectory.point_from_proportion(9/9), color=RED, radius=0.08)
        self.play(Flash(dot3, flash_radius=0.3), FadeIn(dot3))
        
        # Update time label to 3T
        self.play(Transform(time_label, Text("3T", font_size=36).next_to(progress_bar_border, LEFT, buff=0.3)))
        
        # Make the blue trajectory disappear to show only the Poincaré points
        self.play(FadeOut(trajectory), run_time=1)
        
        self.wait(2)
