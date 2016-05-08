#DCVIZ

from numpy import sin, cos,  array
from os.path import join, split
import re

from DCViz_sup import DCVizPlotter


class GeneticAlorithmScala(DCVizPlotter):

    nametag = "genetic(\d+)\.dat"
    isFamilyMember = True
    transpose = True
    loadSequential = True
    #loadLatest = True
    ziggyMagicNumber = 1

    figMap = {"figure" : "subfigure", "figure2": "subfigure2"}

    hugifyFonts = True

    fig_size = [8, 6]

    @staticmethod
    def eval_fourier_series(x, a0, cos_coeffs, sin_coeffs):
        f = a0
        for n, (c, s) in enumerate(zip(cos_coeffs, sin_coeffs)):
            f += s*sin((n+1)*x) + c*cos((n+1)*x)

        return f

    def plotsingle(self, x, entry):
        a0 = float(self.loader.get_metadata()[1][0].strip())
        cos_coeffs, sin_coeffs = entry
        self.subfigure.plot(x, self.eval_fourier_series(x, a0, cos_coeffs, sin_coeffs))

        self.subfigure2.plot([0, 0], [0, a0], "k-", linewidth=10)

        for n, coeff in enumerate(sin_coeffs):
            self.subfigure2.plot([n+1, n+1], [0, coeff], "r-", linewidth=10)

        for n, coeff in enumerate(cos_coeffs):
            self.subfigure2.plot([n+0.5, n+0.5], [0, coeff], "b-", linewidth=10)


    def plot(self, data):

        x = []
        f = []

        dir = split(self.filepath)[0]
        with open(join(dir, "target_data.dat"), 'r') as _file:
            for line in _file:
                xi, fi = line.split()
                x.append(float(xi))
                f.append(float(fi))

        x = array(x)
        f = array(f)

        if not (self.loadSequential or self.loadLatest):
            for entry in data:
                self.plotsingle(x, entry)
        else:
            self.plotsingle(x, data)

        n = re.findall(self.nametag, self.filename)[0]
        self.subfigure.set_title("Generation %s" % n)

        self.subfigure.plot(x, f)

        span = f.max() - f.min()
        self.subfigure.set_ylim(f.min() - span/3, f.max() + span/3)
        self.subfigure.set_xlim(x.min(), x.max())

        self.subfigure.set_ylabel("Function value")
        self.subfigure.set_xlabel("x")

        self.subfigure2.set_xlabel("n")
        self.subfigure2.set_ylabel("Wave coefficient")

        self.subfigure2.set_xlim(-0.15, self.subfigure2.get_xlim()[1]+0.15)
        self.subfigure2.set_ybound(0)
