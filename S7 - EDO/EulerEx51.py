import numpy as np

Xstart = 0; Xend = 1
Ustart = 0; Uend = 0.5*(np.exp(Xend) - np.sin(Xend) - np.cos(Xend))

Nbh = 16

Eexpl = np.zeros(Nbh)
Eimpl = np.zeros(Nbh)

for i in range(Nbh):
    n= pow(2, i+1)
    h = (Xend - Xstart)/n
    X = np.linspace(Xstart, Xend, n+1)
    Uexpl = np.zeros(n+1)   ; Uexpl[0] = Ustart
    Uimpl = np.zeros(n+1)   ; Uimpl[0] = Ustart

    for j in range(n):
        Uexpl[j+1] = Uexpl[j] * (1 + h) + h * np.sin(X[j])
        Uimpl[j+1] = (Uimpl[j] + h* np.sin(X[j+1]))/(1-h)
    Eexpl[i] = abs(Uend-Uexpl[-1])
    Eimpl[i] = abs(Uend-Uimpl[-1])

LogExpl = np.log(abs(Eexpl[:-1]/Eexpl[1:]))/np.log(2)
LogImpl = np.log(abs(Eimpl[:-1]/Eimpl[1:]))/np.log(2)

print("Echelle logarithmiques Erreur explicite", '\n',*['%.4f' % val for val in LogExpl])
print("Echelle logarithmiques Erreur implicite", '\n',*['%.4f' % val for val in LogImpl])