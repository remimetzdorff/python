import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import datetime
from celestialbody import CelestialBody

def prepare_data(names, ref=None, start=None, stop=None, step=None):
    Xs, Ys, Zs = [], [], []
    for name in names:
        body = CelestialBody(name)
        X, Y, Z, _ = body.trajectory(start=start, stop=stop, step=step)
        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
    X_sun, Y_sun, Z_sun = np.zeros(len(X)), np.zeros(len(X)), np.zeros(len(X))
    Xs.append(X_sun)
    Ys.append(Y_sun)
    Zs.append(Z_sun)
    names.append("Sun")
    if ref not in names:
        print("ref must be in names!")
        return
    else:
        for i, name in enumerate(names):
            if name == ref:
                Xs.append(Xs[i])
                Ys.append(Ys[i])
                Zs.append(Zs[i])
        names.append(ref)

    Xs_ref, Ys_ref, Zs_ref = [], [], []
    for X, Y, Z in zip(Xs, Ys, Zs):
        Xs_ref.append(X - Xs[-1])
        Ys_ref.append(Y - Ys[-1])
        Zs_ref.append(Z - Zs[-1])

    coord_sun = [Xs, Ys, Zs]
    coord_ref = [Xs_ref, Ys_ref, Zs_ref]
    return coord_sun, coord_ref

def frame_of_reference(names, ref=None, start = datetime.datetime(2021, 1, 1), stop = datetime.datetime(2029, 1, 1), step = 2):
    """
    Affiche deux graphes contenant les trajectoires des objets contenus dans la liste names :
        - le premier dans le référentiel écliptique héliocentrique
        - le deuxième dans le référentiel écliptique lié à l'objet nommé ref
    """
    coord_sun, coord_ref = prepare_data(names, ref=ref, start = start, stop = stop, step = step)

    fig = plt.figure(figsize=(16, 8))
    sps = (1, 2)
    ax1 = plt.subplot2grid(sps, (0, 0))
    ax2 = plt.subplot2grid(sps, (0, 1))

    for coord, ax in zip([coord_sun, coord_ref], [ax1, ax2]):
        Xs, Ys, Zs = coord
        i = 0
        for name, X, Y, Z in zip(names[:-1], Xs[:-1], Ys[:-1], Zs[:-1]):
            if name == "Sun":
                color = "gold"
            else:
                color = "C"+str(i)
            position,   = ax.plot(X[0], Y[0], "o", color=color, label=name)
            trajectory, = ax.plot(X, Y, "-", color=color, linewidth=.5)
            i += 1

        ax.set_xlabel("X (au)")
        ax.set_ylabel("Y (au)")
        ax.set_aspect("equal")
        ax.legend()
    ax1.set_xlim(ax2.get_xlim())
    ax1.set_ylim(ax2.get_ylim())
    return fig