import matplotlib.pyplot as plt
import numpy as np


def l2_score(x, y):
    return np.sqrt(np.sum(np.square(x - y), axis=1))


def plot_side(arr, ylim=[-2, 2.5], yticks=[0, 2], pic_name=None):
    """
    Adapted from APARENT code (Bogard et al. 2019)
    """
    assert arr.shape[0] % 2 == 0, "arr must have even length."
    midpoint = int(arr.shape[0] / 2)
    pl = arr[:midpoint]
    mn = arr[midpoint:]
    plt.bar(
        range(pl.shape[0]),
        pl,
        width=-2,
        color="r",
    )
    plt.bar(range(mn.shape[0]), -mn, width=-2, color="b")
    axes = plt.gca()
    axes.set_ylim(ylim)
    axes.set_yticks(yticks)
    axes.set_xticks([])
    axes.spines[["right", "top", "bottom"]].set_visible(False)
    plt.xlim(-0.5, pl.shape[0] - 0.5)
    axes.tick_params(labelleft=False)

    if pic_name is None:
        plt.show()
    else:
        plt.savefig(pic_name, transparent=True)
        plt.close()
