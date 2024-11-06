from manim import *
from manim_slides import Slide, ThreeDSlide
from manim.opengl import *

stdcolor = BLACK
background_color = WHITE
if stdcolor != BLACK:
    background_color = BLACK
slidetitlepos = 3*UP + 2*LEFT
FSDOCTITLE = 1.2
FSTITLE = 1.05
FSTEXT = 0.95

class SceneTCC(Slide):
    def construct(self):
        self.camera.background_color = background_color
        def doslidetitle(string, section=None):
            if section is None:
                return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)
            else:
                return Tex("{\\small "+section+":} \\textbf{"+string+"}", color=stdcolor).scale(FSTITLE)
        def doslideframe(string):
            return Tex("\\textbf{"+string+"}", color=stdcolor).scale(FSDOCTITLE*1.2)

        sN = 0
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


# class Orbiting(ThreeDSlide):
#     def construct(self):
