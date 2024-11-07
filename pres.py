from manim import *
from manim_slides import Slide, ThreeDSlide
import numpy as np
from numpy import cos, sin
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
        return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)
        # return Tex("{\\small "+section+":} \\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)

def doslideframe(string):
    return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSDOCTITLE*1.2)


class SceneTCC1(Slide):
    def construct(self):
        sN = 0
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
            \\item Particle Tracking
            \\item Path Length
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

class SceneTCC2(ThreeDSlide):
    def construct(self):
        sN = 2
        self.camera.background_color = background_color
        self.set_camera_orientation(phi=70*DEGREES)


        #! parametros

        R = 3.0
        t = 0.05
        runtime = 7
        fps = 30
        dtk = 1/fps
        dtsk = -2*dtk


        #! titulo

        #* Slide

        # slideNUMBER = Integer(sN).scale(0.6).set_color(BLACK).to_corner(DOWN + RIGHT)
        # slideNUMBER.add_updater(
        #     lambda mob: mob.set_value(sN)
        # )
        slideNUMBER = MathTex(str(sN)).scale(0.6).set_color(BLACK).to_corner(DOWN + RIGHT)

        self.next_slide()
        slideNUMBER.update()
        self.add_fixed_in_frame_mobjects(slideNUMBER)
        ftitle = "Theoretical Review"
        slidetitle = doslidetitle("Synchrotrons", section=ftitle).to_corner(UP+LEFT)
        temp = slidetitle
        self.add_fixed_in_frame_mobjects(temp)
        # self.remove(temp)
        # self.add_fixed_in_frame_mobjects(temp)
        slidetitle = doslidetitle("Coordinate System", section=ftitle).to_corner(UP+LEFT)
        self.add_fixed_in_frame_mobjects(slidetitle)
        self.remove(slidetitle)
        self.play(Unwrite(temp), Write(slidetitle, reverse=True))
        # self.add_fixed_in_frame_mobjects(slidetitle)


        # title = doslidetitle("Particle Dynamics", "Theoretical Review").to_corner(UP+LEFT)
        # self.add_fixed_in_frame_mobjects(title)
        # self.remove(title)
        # self.play(Write(title))


        #! objetos principais

        orbit = Circle(radius=R, color=stdcolor, stroke_width=0.8)
        orbdraw = Dot(radius=0.001, color=background_color)
        particle = Dot3D(radius=0.05, color=RED)

        def ideal_p_updater(mob, dt):
            nonlocal t, R, runtime
            t -= dt
            tt = t + 0.0
            mob.move_to(np.array([R*cos(2*PI*tt/runtime), R*sin(2*PI*tt/runtime), 0]))

        def real_p_updater(mob, dt):
            nonlocal t, R, runtime
            np.random.seed(int(abs(t*1000)))
            rn = np.random.rand(2)
            ampx = R * (1 + 0.1*cos(2*PI * t + rn[0]/10))
            ampy = 0.1*cos(2*PI * t + rn[1]/10 - 0.3)
            tt = t -1*dt - 0.01*np.random.rand()
            mob.move_to(np.array([ampx*cos(2*PI*tt/runtime), ampx*sin(2*PI*tt/runtime), ampy]))


        #! desenhando a orbita ideal


        orbpath = TracedPath(orbdraw.get_center, stroke_width=0.8, stroke_color=stdcolor)

        # self.add(orbdraw)
        orbdraw.add_updater(ideal_p_updater)
        orbdraw.update()

        def orbvec_updater():
            start = np.array([0,0,0])
            end = orbdraw.get_center()
            direction = (end - start)
            mag = norm(direction)
            # if mag == 0:
            #     direction = [1, 0, 0]
            # else:
            direction /= mag
            return Vector(direction=direction, buff = 0.2, color=GRAY, stroke_width=0.8).put_start_and_end_on(start, end)

        orbvec = always_redraw(orbvec_updater)

        def orbstr_updater():
            end = orbvec.get_end()
            stt = orbvec.get_start()
            di = (end-stt)
            return MathTex("\\vec{r}_0(s)", color=GRAY).scale(0.5).rotate(PI/2, axis=RIGHT).move_to(end - di/2 + 0.5*UP)

        orbstr = always_redraw(orbstr_updater)

        #* criando vetor da orbita
        self.next_slide()
        self.add(slidetitle)
        orbdraw.remove_updater(ideal_p_updater)
        self.play(Create(orbvec), run_time=0.7, rate_func=linear)
        self.play(Write(orbstr), runtime=0.7)
        # self.remove(orbdraw)

        #* desenhando a orbita
        self.next_slide()
        # runtime = 5
        # self.add(orbdraw)
        orbdraw.add_updater(ideal_p_updater)
        self.add(orbpath)
        self.play(orbdraw.animate, run_time=runtime/3, rate_func=linear)
        orbdraw.remove_updater(ideal_p_updater)

        #* nomeando a orbita ideal
        # self.next_slide()
        self.remove(orbpath)

        self.add(orbit) # adiciona a órbita estática
        orbittext = Tex("design orbit", color=stdcolor).scale(0.5).rotate(PI/2, axis=RIGHT).move_to([-R, -R, 0])
        self.play(Write(orbittext))

        # #* retirando o texto da orbita ideal
        # self.next_slide()
        # self.play(FadeOut(orbittext))


        #! definindo frenet-serret

        #* nomeando a orbita ideal
        self.next_slide()
        frenet = Tex("Frenet-Serret\\\\ Coordinate System", color=stdcolor).scale(0.65)
        self.add_fixed_in_frame_mobjects(frenet)
        self.remove(frenet)
        frenet.shift(2.5*DOWN)
        self.play(Write(frenet), FadeOut(orbittext))

        #* updaters para os versores
        def shat_updater():
            nonlocal t, R, runtime
            start = np.array([R*cos(2*PI*t/runtime), R*sin(2*PI*t/runtime), 0])
            tnext = t - dtk
            end = np.array([R*cos(2*PI*tnext/runtime), R*sin(2*PI*tnext/runtime), 0])
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start,start+direction)

        def xhat_updater():
            nonlocal t, R, runtime
            start = np.array([R*cos(2*PI*t/runtime), R*sin(2*PI*t/runtime), 0])
            Rplus = R + 0.8
            end = np.array([Rplus*cos(2*PI*t/runtime), Rplus*sin(2*PI*t/runtime), 0])
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start, start+direction)

        def yhat_updater():
            nonlocal t, R, runtime
            start = np.array([R*cos(2*PI*t/runtime), R*sin(2*PI*t/runtime), 0])
            dy = 0.8
            end = np.array([R*cos(2*PI*t/runtime), R*sin(2*PI*t/runtime), dy])
            direction = (end - start)
            direction /= norm(direction)
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start,start+direction)

        #* versores
        shat = always_redraw(shat_updater)
        xhat = always_redraw(xhat_updater)
        yhat = always_redraw(yhat_updater)

        #* updaters para labels dos versores
        def sstr_updater():
            end = shat.get_end()
            stt = shat.get_start()
            di = (end-stt)
            return MathTex("\\hat{s}", color=stdcolor,).scale(0.4).rotate(PI/2, axis=RIGHT).move_to(end + 0.3*di)
        def xstr_updater():
            end = xhat.get_end()
            stt = xhat.get_start()
            di = (end-stt)
            return MathTex("\\hat{x}", color=stdcolor,).scale(0.4).rotate(PI/2, axis=RIGHT).move_to(end + 0.3*di)
        def ystr_updater():
            end = yhat.get_end()
            stt = yhat.get_start()
            di = (end-stt)
            return MathTex("\\hat{y}", color=stdcolor,).scale(0.4).rotate(PI/2, axis=RIGHT).move_to(end + 0.3*di)

        #* labels dos versores
        shatstr = always_redraw(sstr_updater)
        xhatstr = always_redraw(xstr_updater)
        yhatstr = always_redraw(ystr_updater)

        #* definindo s
        self.next_slide()
        self.play(Create(shat))
        self.play(Write(shatstr))

        shattext = MathTex("\\hat{s} = \\dfrac{d\\vec{r}_0}{ds}", color=stdcolor).scale(0.6).to_corner(RIGHT).shift(UP)
        self.add_fixed_in_frame_mobjects(shattext)
        self.remove(shattext)
        self.play(Write(shattext))

        #* definindo x
        self.next_slide()
        self.play(Create(xhat))
        self.play(Write(xhatstr))

        xhattext = MathTex("\\hat{x} = -\\rho\\dfrac{d\\hat{s}}{ds}", color=stdcolor).scale(0.6).to_corner(RIGHT).next_to(shattext, DOWN)
        self.add_fixed_in_frame_mobjects(xhattext)
        self.remove(xhattext)
        self.play(Write(xhattext))

        #* definindo y
        self.next_slide()
        self.play(Create(yhat))
        self.play(Write(yhatstr))

        yhattext = MathTex("\\hat{y} = \\hat{s}\\times\\hat{x}", color=stdcolor).scale(0.6).to_corner(RIGHT).next_to(xhattext, DOWN)
        self.add_fixed_in_frame_mobjects(yhattext)
        self.remove(yhattext)
        self.play(Write(yhattext))


        #! simulação

        #* addiciona partícula
        self.next_slide()
        particle.add_updater(real_p_updater)
        particle.update()
        particle.remove_updater(real_p_updater)
        self.play(FadeIn(particle))

        #* run
        # self.next_slide(loop=True)
        # orbdraw.add_updater(ideal_p_updater)
        # particle.add_updater(real_p_updater)
        # particle.update()
        # orbdraw.update()
        path = TracedPath(particle.get_center, dissipating_time=1.2, stroke_opacity=[0, 1], stroke_color=RED, stroke_width=2)
        self.add(path)
        # self.play(particle.animate, run_time=1*runtime, rate_func=linear)


        self.next_slide()
        orbdraw.remove_updater(ideal_p_updater)
        particle.remove_updater(real_p_updater)
        self.remove(path)

        self.play(particle.animate, FadeOut(particle),FadeOut(orbit),
        FadeOut(orbdraw),FadeOut(orbvec), FadeOut(orbstr),
        FadeOut(shat),FadeOut(shatstr),FadeOut(shattext),
        FadeOut(xhat),FadeOut(xhatstr),FadeOut(xhattext),
        FadeOut(yhat),FadeOut(yhatstr),FadeOut(yhattext),
        FadeOut(frenet))#, FadeOut(slidetitle), FadeOut(slideNUMBER))

        # self.next_slide()

        sN += 1
        self.remove(slideNUMBER)
        slideNUMBER = MathTex(str(sN)).scale(0.6).set_color(BLACK).to_corner(DOWN + RIGHT)
        self.add_fixed_in_frame_mobjects(slideNUMBER)
        self.remove(slidetitle)
        temp = slidetitle
        slidetitle = doslidetitle("Longitudinal Motion").to_corner(UP+LEFT)
        self.add_fixed_in_frame_mobjects(temp)
        self.add_fixed_in_frame_mobjects(slidetitle)
        self.remove(slidetitle)
        self.play(Unwrite(temp), Write(slidetitle, reverse=True))



class SceneTCC3(Slide):
    def construct(self):
        sN = 3
        self.camera.background_color = background_color

        slideNUMBER = Integer(sN).scale(0.6).to_corner(DOWN + RIGHT).set_color(stdcolor)
        slideNUMBER.add_updater(
            lambda mob: mob.set_value(sN)
        )
        #!########### PATH LENGHT X ############

        # #* Section Frame
        # self.next_slide()
        sN += 1
        slideNUMBER.update()
        self.add(slideNUMBER)
        temp = doslidetitle("Longitudinal Motion", section="Theoretical Review").to_corner(UP+LEFT)
        self.add(temp)
        slidetitle = doslidetitle("Synchrotron Oscillations", section="Theoretical Review").to_corner(UP+LEFT)
        self.play(Unwrite(temp), Write(slidetitle, reverse=True))

        self.next_slide()
        self.add(slidetitle)
        testtext = MathTex("f(x) = \\dfrac{x^2}{2\\!}", color=stdcolor).scale(0.7)
        self.play(Write(testtext))

        self.next_slide()
        self.play(FadeOut(testtext))

        # #* Section Frame
        self.next_slide()
        # slidetitle = doslidetitle("Particle Dynamics", section="Theoretical Review").to_corner(UP+LEFT)
        # # self.add(slidetitle)
        ftitle = "Particle Tracking"
        capframe = doslideframe(ftitle)
        slideNUMBER.set_color(background_color)
        self.play(ReplacementTransform(slidetitle, capframe))
        # self.play(Write(capframe))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.set_color(stdcolor).update()
        self.add(slideNUMBER)
        slidetitle = doslidetitle("Maps", section=ftitle).to_corner(UP+LEFT)
        self.play(ReplacementTransform(capframe, slidetitle))

        #* Slide
        sN += 1
        self.next_slide()
        slideNUMBER.update()
        temp = slidetitle
        self.remove(slidetitle)
        # self.add(temp)
        slidetitle = doslidetitle("1-turn Map", section=ftitle).to_corner(UP+LEFT)
        self.play(Unwrite(temp), Write(slidetitle, reverse=True))

        sN += 1
        self.next_slide()
        slideNUMBER.update()
        temp = slidetitle
        self.remove(slidetitle)
        # self.add(temp)
        slidetitle = doslidetitle("Closed Orbit", section=ftitle).to_corner(UP+LEFT)
        self.play(Unwrite(temp), Write(slidetitle, reverse=True))

        #!########### END ############
        self.next_slide()

        end = doslidetitle("End").scale(FSDOCTITLE)
        self.play(FadeOut(slideNUMBER), ReplacementTransform(slidetitle, end))
