#DCVIZ

from numpy import sin, cos, linspace, pi

from DCViz_sup import DCVizPlotter


class GeneticAlorithmScala(DCVizPlotter):

    nametag = "genetic(\d+)\.dat"
    isFamilyMember = True
    transpose = True
    loadLatest = True
    ziggyMagicNumber = 1

    xmax = 1
    #xmax = 2*pi

    def analytical(self, x):
        return 0 if x < 0.5 else 1
        #return sin(x) + 5*cos(2*x) + 2*sin(5*x)

    @staticmethod
    def eval_fourier_series(x, cos_coeffs, sin_coeffs):
        f = 0
        for n, (c, s) in enumerate(zip(cos_coeffs, sin_coeffs)):
            f += s*sin((n+1)*x) + c*cos((n+1)*x)

        return f

    def plotsingle(self, x, entry):
        cos_coeffs, sin_coeffs = entry
        self.subfigure.plot(x, self.eval_fourier_series(x, cos_coeffs, sin_coeffs))

    def plot(self, data):

        x = linspace(0, self.xmax, 1000)

        if not (self.loadSequential or self.loadLatest):
            for entry in data:
                self.plotsingle(x, entry)
        else:
            self.plotsingle(x, data)
            self.subfigure.set_title(self.filename)

        self.subfigure.plot(x, [self.analytical(xi) for xi in x])

        self.subfigure.set_ylim(-0.5, 1.5)

