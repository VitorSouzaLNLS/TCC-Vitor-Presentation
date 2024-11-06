from manim import *
from manim_slides import Slide, ThreeDSlide
import numpy as np
from numpy import pi, cos, sin
from numpy.linalg import norm

stdcolor = BLACK
background_color = WHITE
if stdcolor != BLACK:
    background_color = BLACK
slidetitlepos = 3*UP + 2*LEFT
FSDOCTITLE = 1.2
FSTITLE = 1.05
FSTEXT = 0.95

def doslidetitle(string, section=None):
            if section is None:
                return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)
            else:
                return Tex("{\\small "+section+":} \\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)
def doslideframe(string):
    return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSDOCTITLE*1.2)

sN = 0

class SceneTCC(Slide):
    def construct(self):
        self.camera.background_color = background_color

        slideNUMBER = Integer(sN).scale(0.6).to_corner(DOWN + RIGHT).set_color(background_color)
        slideNUMBER.add_updater(
            lambda mob: mob.set_value(sN)
        )

        #!########### TITLE PAGE ############
        title = Tex("\\centering \\textbf{Improvement on Particle Tracking Tool\\\\ for Accelerator Simulations at LNLS}", color=stdcolor).scale(FSDOCTITLE)
        self.play(Write(title), run_time=2)

        vitor = Tex("\\centering Vitor Davi de Souza", color=stdcolor).scale(FSTITLE)
        self.play(title.animate.shift(1.4*UP), run_time=0.5)
        self.play(FadeIn(vitor), run_time=1)
        self.next_slide()


        #!########### CONTENTS ############
        contents = doslidetitle("Contents").to_corner(UP)

        self.play(FadeOut(vitor), ReplacementTransform(title, contents))

        itensc = Tex("""
        \\begin{itemize}
            \\item Theoretical Review
            \\item Path Length Limitation
            \\item Path Length Modification
            \\item Results
            \\item Conclusion
        \\end{itemize}""", color=stdcolor).scale(FSTEXT*0.7)
        # itensc.next_to(contents, DOWN)
        self.play(Write(itensc))

        #!########### THEORETICAL REVIEW ############

        #* Section Frame
        self.next_slide()
        ftitle = "Theoretical Review"
        capframe = doslideframe(ftitle)
        self.play(FadeOut(itensc), ReplacementTransform(contents, capframe))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.set_color(stdcolor).update()
        slidetitle = doslidetitle("Synchrotrons", section=ftitle).to_corner(UP+LEFT)
        self.play(FadeIn(slideNUMBER), FadeTransform(capframe, slidetitle))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.update()
        temp = slidetitle
        slidetitle = doslidetitle("Particle Dynamics", section=ftitle).to_corner(UP+LEFT)
        self.play(FadeTransform(temp, slidetitle))

        #!########### PATH LENGHT X ############

        #* Section Frame
        self.next_slide()
        ftitle = "Path Length X"
        capframe = doslideframe(ftitle)
        slideNUMBER.set_color(background_color)
        self.play(ReplacementTransform(slidetitle, capframe))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.set_color(stdcolor).update()
        slidetitle = doslidetitle("Issue", section=ftitle).to_corner(UP+LEFT)
        self.play(FadeTransform(capframe, slidetitle))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.update()
        temp = slidetitle
        slidetitle = doslidetitle("Particle Dynamics", section=ftitle).to_corner(UP+LEFT)
        self.play(FadeTransform(temp, slidetitle))


        #!########### END ############
        self.next_slide()

        end = doslidetitle("End").scale(FSDOCTITLE)
        self.play(FadeOut(slideNUMBER), ReplacementTransform(slidetitle, end))


class Orbiting(ThreeDSlide):
    def construct(self):
        self.camera.background_color = background_color

        runtime = 7
        fps = 30
        dtk = 1/fps

        # axes = ThreeDAxes()
        # self.add(axes)

        title = doslidetitle("Teste 1").to_corner(UP+LEFT)
        self.add_fixed_in_frame_mobjects(title)

        R = 3.0
        orbit = Circle(radius=R, color=stdcolor, stroke_width=0.8)
        particle = Dot3D(radius=0.05, color=RED)
        t = 0

        def ideal_p_updater(mob, dt):
            nonlocal t, R, runtime
            t -= dt
            tt = t
            mob.move_to([R*cos(2*pi*tt/runtime), R*sin(2*pi*tt/runtime), 0])

        def real_p_updater(mob, dt):
            nonlocal t, R, runtime
            t -= dt
            ampx = R * (1 + 0.1*cos(2*pi * 0.7 * t + np.random.rand()/10))
            ampy = 0.1*cos(2*pi * 0.4 * t + np.random.rand()/10)
            tt = t -2*dt - 0.01*np.random.rand()
            mob.move_to([ampx*cos(2*pi*tt/runtime), ampx*sin(2*pi*tt/runtime), ampy])


        #* criando a particula e a orbita ideal
        self.set_camera_orientation(phi=70*DEGREES)
        self.play(GrowFromCenter(orbit))
        particle.update(0)
        self.play(FadeIn(particle))

        #* run da part. na orbita ideal
        self.next_slide(loop=True)
        particle.add_updater(ideal_p_updater)
        path = TracedPath(particle.get_center, dissipating_time=0.3, stroke_opacity=[0, 1], stroke_color=RED, stroke_width=3)
        self.add(path)
        self.play(particle.animate, run_time=runtime)

        #* run da part numa orbita "real"
        self.next_slide(loop=True)
        # self.remove(path)
        particle.remove_updater(ideal_p_updater)
        particle.add_updater(real_p_updater)
        # path = TracedPath(particle.get_center, dissipating_time=0.3, stroke_opacity=[0, 1], stroke_color=RED, stroke_width=2)
        # self.add(path)
        self.play(particle.animate, run_time=runtime)

        #* definicao dos sistema de coordenadas
        self.next_slide()

        def shat_updater():
            nonlocal t, R
            start = np.array([R*cos(2*pi*t/runtime), R*sin(2*pi*t/runtime), 0])
            tnext = t - dtk
            end = np.array([R*cos(2*pi*tnext/runtime), R*sin(2*pi*tnext/runtime), 0])
            # return Arrow(start=start, end=end, color=stdcolor).scale(3)
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start,start+direction)
        def xhat_updater():
            nonlocal t, R
            start = np.array([R*cos(2*pi*t/runtime), R*sin(2*pi*t/runtime), 0])
            Rplus = R + 0.8
            end = np.array([Rplus*cos(2*pi*t/runtime), Rplus*sin(2*pi*t/runtime), 0])
            # return Arrow(start=end, end=start, color=stdcolor, buff=0.5)
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start, start+direction)
        def yhat_updater():
            nonlocal t, R
            start = np.array([R*cos(2*pi*t/runtime), R*sin(2*pi*t/runtime), 0])
            dy = 0.8
            end = np.array([R*cos(2*pi*t/runtime), R*sin(2*pi*t/runtime), dy])
            # return Arrow(start=end, end=start, color=stdcolor)
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start,start+direction)

        shat = always_redraw(shat_updater)
        xhat = always_redraw(xhat_updater)
        yhat = always_redraw(yhat_updater)

        def sstr_updater():
            end = shat.get_end()
            stt = shat.get_start()
            di = (end-stt)
            return MathTex("\\hat{s}", color=stdcolor,).scale(0.4).rotate(pi/2, axis=RIGHT).move_to(end + 0.3*di)
        def xstr_updater():
            end = xhat.get_end()
            stt = xhat.get_start()
            di = (end-stt)
            return MathTex("\\hat{x}", color=stdcolor,).scale(0.4).rotate(pi/2, axis=RIGHT).move_to(end + 0.3*di)
        def ystr_updater():
            end = yhat.get_end()
            stt = yhat.get_start()
            di = (end-stt)
            return MathTex("\\hat{y}", color=stdcolor,).scale(0.4).rotate(pi/2, axis=RIGHT).move_to(end + 0.3*di)

        shatstr = always_redraw(sstr_updater)
        xhatstr = always_redraw(xstr_updater)
        yhatstr = always_redraw(ystr_updater)

        particle.remove_updater(real_p_updater)

        self.play(FadeIn(shat), FadeIn(xhat), FadeIn(yhat))

        self.play(Write(shatstr),Write(xhatstr),Write(yhatstr))

        #* run particula + orbita real + sistema de coordenadas
        self.next_slide(loop=True)
        particle.add_updater(real_p_updater)

        self.play(particle.animate, run_time=runtime)
