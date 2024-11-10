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
FSTEXT = 0.95*0.8

# np.random.seed(int(abs(t*1000)))
np.random.seed(1297482)

def dosectitle(string, section=None, align=None):
    if align is None:
        ali = UP+LEFT
    else:
        ali = align
    return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSTITLE).to_corner(ali)
    # else:
    #     return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSTITLE).to_corner(ali)
        # return Tex("{\\small "+section+":} \\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)

def docaptitle(string):
    return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSDOCTITLE*1.2)

def doslidenumber(n):
    slideNUMBER = Integer(n).scale(0.6).to_corner(DOWN + RIGHT).set_color(stdcolor)
    # def slideupdater(mob, _n):
    #     mob.set_value(_n)
    # slideNUMBER.add_updater(slideupdater)
    return slideNUMBER

def dosec2sec(scene, curr_sec_mob, next_sec_title, proj=None):
    nextsec = dosectitle(next_sec_title)
    if proj == '3d':
        scene.add_fixed_in_frame_mobjects(nextsec)
        scene.remove(nextsec)
    scene.play(Unwrite(curr_sec_mob, reverse=True), run_time=0.6)
    scene.play(Write(nextsec))
    # return nextsec

def docap2sec(scene, curr_cap_mob, next_sec_title):
    nextsec = dosectitle(next_sec_title)
    scene.play(FadeTransform(curr_cap_mob, nextsec))
    # return nextsec

class S1(Slide):
    def construct(self):
        self.camera.background_color = background_color

        #!########### TITLE PAGE ############
        self.next_slide()
        title = Tex("\\centering \\textbf{Improvement on particle tracking Tool\\\\ for accelerator simulations at LNLS}", color=stdcolor).scale(FSDOCTITLE)
        self.play(Write(title), run_time=2)

        vitor = Tex("\\centering Vitor Davi de Souza", color=stdcolor).scale(FSTITLE)
        self.play(title.animate.shift(1.4*UP), run_time=0.5)
        self.play(FadeIn(vitor), run_time=1)

        data = Tex("\\centering November 10, 2024", color=stdcolor).scale(FSTEXT).next_to(vitor, DOWN).shift(DOWN)
        self.play(Write(data), run_time=0.5)

        #!########### CONTENTS ############
        self.next_slide()
        contents = dosectitle("Contents", align=UP)

        self.play(FadeOut(vitor), FadeOut(data), ReplacementTransform(title, contents))

        itensc = Tex("""
        \\begin{itemize}
            \\item Theoretical Review
            \\item Particle Tracking
            \\item Path Length
            \\item Results
            \\item Conclusion
        \\end{itemize}""", color=stdcolor).scale(FSTEXT)
        # itensc.next_to(contents, DOWN)
        self.play(Write(itensc))

        #!########### THEORETICAL REVIEW ############

        #* Section Frame
        self.next_slide()
        captitle = docaptitle("Theoretical Review")
        self.play(FadeOut(itensc), ReplacementTransform(contents, captitle))

        #* Slide
        self.next_slide()
        slideNUMBER = doslidenumber(1)
        docap2sec(self, captitle, "Synchrotrons")
        self.play(FadeIn(slideNUMBER), run_time=0.5)


class S2(Slide):
    def construct(self):
        self.camera.background_color = background_color
        slideNUMBER = doslidenumber(1)
        self.add(slideNUMBER)
        slidetitle = dosectitle("Synchrotrons")
        self.add(slidetitle)

        #* END this
        self.next_slide()
        slideNUMBER.set_value(2)
        dosec2sec(self, slidetitle, "Coordinate System")


class S3(ThreeDSlide):
    def construct(self):
        self.camera.background_color = background_color
        self.set_camera_orientation(phi=70*DEGREES)

        #! titulo
        #* Slide

        slideNUMBER = doslidenumber(2)
        self.add_fixed_in_frame_mobjects(slideNUMBER)
        slidetitle = dosectitle("Coordinate System")
        self.add_fixed_in_frame_mobjects(slidetitle)

        #! parametros

        R = 3.0
        t = 0.05
        runtime = 7
        fps = 30
        dtk = 1/fps
        dtsk = -2*dtk

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
            rn = np.random.rand(2)
            ampx = R * (1 + 0.1*cos(2*PI * t + rn[0]/10))
            ampy = 0.1*cos(4*PI * t + rn[1]/10 - 0.3)
            tt = t -1*dt - 0.02*np.random.rand()
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
            return Vector(direction=direction, buff = 0.2, color=stdcolor, stroke_width=0.8).put_start_and_end_on(start,start+direction).rotate(PI/2, OUT)

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
        shattext = MathTex("\\hat{s} = \\dfrac{d\\vec{r}_0}{ds}", color=stdcolor).scale(0.6).to_corner(RIGHT).shift(UP)
        self.add_fixed_in_frame_mobjects(shattext)
        self.remove(shattext)

        self.play(Create(shat))
        self.play(Write(shatstr), Write(shattext))


        #* definindo x
        self.next_slide()
        xhattext = MathTex("\\hat{x} = -\\rho\\dfrac{d\\hat{s}}{ds}", color=stdcolor).scale(0.6).to_corner(RIGHT).next_to(shattext, DOWN)
        self.add_fixed_in_frame_mobjects(xhattext)
        self.remove(xhattext)
        self.play(Create(xhat))
        self.play(Write(xhatstr), Write(xhattext))

        #* definindo y
        self.next_slide()
        yhattext = MathTex("\\hat{y} = \\hat{s}\\times\\hat{x}", color=stdcolor).scale(0.6).to_corner(RIGHT).next_to(xhattext, DOWN)
        self.add_fixed_in_frame_mobjects(yhattext)
        self.remove(yhattext)
        self.play(Create(yhat))
        self.play(Write(yhatstr), Write(yhattext))

        print(yhat.get_start(), yhat.get_end())

        return
        #! simulação

        #* addiciona partícula
        self.next_slide()
        particle.add_updater(real_p_updater)
        particle.update()
        particle.remove_updater(real_p_updater)
        self.play(FadeIn(particle))

        #* run
        self.next_slide(loop=True)
        orbdraw.add_updater(ideal_p_updater)
        particle.add_updater(real_p_updater)
        particle.update()
        orbdraw.update()
        path = TracedPath(particle.get_center, dissipating_time=1.2*runtime, stroke_opacity=[0, 1], stroke_color=RED, stroke_width=2)
        self.add(path)
        self.play(particle.animate, run_time=1*runtime, rate_func=linear)


        self.next_slide()
        orbdraw.remove_updater(ideal_p_updater)
        particle.remove_updater(real_p_updater)
        # self.remove(path)

        ds, dx, dy = 0, 0.25, 0.2
        p0 = np.array([R+dx, 0, dy])
        self.play(
            # particle.animate, FadeOut(particle), FadeOut(orbit),
            FadeOut(path),
            FadeOut(orbdraw),FadeOut(orbvec), FadeOut(orbstr),
            # FadeOut(shat),FadeOut(shatstr),
            FadeOut(shattext),
            # FadeOut(xhat),FadeOut(xhatstr),
            FadeOut(xhattext),
            # FadeOut(yhat),FadeOut(yhatstr),
            FadeOut(yhattext),
            FadeOut(frenet),
            #, FadeOut(slidetitle), FadeOut(slideNUMBER))
            particle.animate.move_to(p0)
        )

class S4(ThreeDSlide):
    def construct(self):
        self.camera.background_color = background_color
        self.set_camera_orientation(phi=70*DEGREES)

        #! titulo
        #* Slide

        slideNUMBER = doslidenumber(3)
        self.add_fixed_in_frame_mobjects(slideNUMBER)
        slidetitle = dosectitle("Coordinate System")
        self.add_fixed_in_frame_mobjects(slidetitle)

        #! parametros

        R = 3.0
        t = 1.0
        runtime = 7
        fps = 30
        dtk = 1/fps
        dtsk = -2*dtk

        #! objetos principais
        orbitfull = Circle(radius=R, color=stdcolor, stroke_width=0.8)
        orbit = Arc(radius=R, start_angle=-PI/4, angle=PI/2, color=stdcolor, stroke_width=0.8)
        ds, dx, dy = 0, 0.25, 0.2
        p0 = np.array([R+dx, 0, dy])
        particle = Dot3D(radius=0.05, color=RED).move_to(p0)
        start = np.array([R, 0, 0])
        d1, d2, d3 = np.array([1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, 1])
        vdx = np.array([-dx,0,0])
        vds = np.array([0,+dx,0])
        vdy = np.array([0,0,-dx])
        ascale=1
        vx = Arrow(start=start+vdx, end=start+d1*ascale, stroke_width=0.8, max_tip_length_to_length_ratio=0.1, color=stdcolor)
        vs = Arrow(start=start+vds, end=start+d2*ascale, stroke_width=0.8, max_tip_length_to_length_ratio=0.1, color=stdcolor)
        vy = Arrow(start=start+vdy, end=start+d3*ascale, stroke_width=0.8, max_tip_length_to_length_ratio=0.1, color=stdcolor).rotate(PI/2, OUT)
        vxstr, vsstr, vystr = [MathTex("\\hat{"+a+"}", color=stdcolor).scale(0.4).rotate(PI/2,RIGHT) for a in ["x", "s", "y"]]
        vxstr.move_to(vx.get_end()+0.3*vx.get_unit_vector())
        vsstr.move_to(vs.get_end()+0.3*vs.get_unit_vector())
        vystr.move_to(vy.get_end()+0.3*vy.get_unit_vector())
        grou = VGroup(orbit, particle, vs, vx, vy, vxstr, vsstr, vystr)
        # self.play(Create(orbit, reverse=True))
        # self.play(Create(vx), Create(vy), Create(vs), Create(vxstr), Create(vystr), Create(vsstr))
        # self.play(FadeIn(particle))
        self.add(orbit)
        self.add(particle)
        self.add(vx, vy, vs, vxstr, vystr, vsstr)
        self.play(FadeOut(orbitfull))

        # self.next_slide()
        # self.move_camera(theta=-110*DEGREES, phi=70*DEGREES)
        # self.play(grou.animate.scale(4).shift(4.5*LEFT))
        # , vxstr.animate.scale(0.6)
        # , vsstr.animate.scale(0.6)
        # , vystr.animate.scale(0.6))

        self.next_slide()
        self.move_camera(theta=-110*DEGREES, phi=70*DEGREES, zoom=3, frame_center=p0)
        self.play(vsstr.animate.scale(0.9), vxstr.animate.scale(0.9), vystr.animate.scale(0.9))

        self.next_slide()
        p0 = particle.get_center()
        oS = vs.get_start()
        proy = DashedLine(start=p0, end=[p0[0], p0[1], oS[2]], color=GRAY, stroke_width=3)
        prox = DashedLine(start=p0, end=[oS[0], oS[1], p0[2]], color=GRAY, stroke_width=3)
        valpy = Line(start=oS, end=[oS[0], oS[1], p0[2]], stroke_width=3, color=GREEN_D)
        valpx = Line(start=oS, end=[p0[0], p0[1], oS[2]], stroke_width=3, color=BLUE_D)
        yval = MathTex("y", color=GREEN_D).scale(0.3).rotate(PI/2, RIGHT).next_to(valpy.get_center(), 0.3*LEFT)#, LEFT)
        xval = MathTex("x", color=BLUE_D).scale(0.3).rotate(PI/2, RIGHT).next_to(valpx.get_center(), 0.3*IN)
        self.play(Create(proy), Create(prox))
        self.play(Create(valpy), Create(valpx), Write(yval), Write(xval))
        self.remove(particle)
        self.add(particle)

        self.next_slide()
        p = np.array([0.07, -0.4, 0.07])*4
        p0 = particle.get_center()
        pvec = Arrow(start=p0, end=p0+p, buff=0, color=RED, stroke_width=3, max_tip_length_to_length_ratio=0.03).rotate(PI/2, p)
        p1 = pvec.get_end()
        pstr = MathTex("p", color=RED).scale(0.3).rotate(PI/2,RIGHT).next_to(p1, OUT)
        self.play(Create(pvec), Write(pstr))
        self.remove(particle)
        self.add(particle)

        self.next_slide()
        sparallel = DashedLine(start=p0, end=[p0[0], p1[1], p0[2]], color=GRAY, stroke_width=2)
        proj = np.array([p1[0], p1[1], p0[2]])
        pparallel = DashedLine(start=p0, end=proj, color=GRAY, stroke_width=2)
        self.play(Create(sparallel), Create(pparallel))
        self.remove(particle)
        self.add(particle)

        self.next_slide()
        xprime = Arc(radius=0.9, start_angle=-PI/2, angle=-p[0]/p[1], arc_center=p0, color=BLUE_D, stroke_width=2)
        xprimestr = MathTex("x^\\prime", color=BLUE_D).scale(0.2).rotate(PI/2, RIGHT).next_to(xprime.get_center(), 0.9*DOWN)
        self.play(Create(xprime), Write(xprimestr))


        self.next_slide()
        pax = pvec.get_unit_vector()
        pax[-1] = 0.0
        pax /= norm(pax)
        ryp=0.95
        pk = particle.get_center() + ryp*pax
        # self.add(Dot3D(pk, radius=0.03, color=BLACK))
        yprime = Arc(radius=ryp, start_angle=-PI/2, angle=-p[2]/p[1], arc_center=p0, color=GREEN_D, stroke_width=2).rotate(-PI/2, pax, about_point=pk)
        yprimestr = MathTex("y^\\prime", color=GREEN_D).scale(0.2).rotate(PI/2, RIGHT).next_to(yprime.get_center(), 0.12*RIGHT+0.25*DOWN)
        self.play(Create(yprime), Write(yprimestr))

        # self.next_slide(loop=True)
        # self.begin_ambient_camera_rotation(rate=0.3, about='theta')
        # self.wait(10)
        # self.next_slide()
        # self.stop_ambient_camera_rotation()

        # dosec2sec(self, slidetitle, "Longitudinal Motion", proj='3d')
        # slideNUMBER.set_value(100)
