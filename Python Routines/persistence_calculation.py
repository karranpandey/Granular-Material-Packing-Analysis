import pyms3d
import numpy as np
import matplotlib.pyplot as plt


def compute_pers_diagm(data_file_name, dim):
    """Comput the persistence diagram

    Args:
        data_file_name (str): raw file with distance field
        dim (tuple): dimensions of distance field

    Returns:
        None: None
    """
    # Comput msc
    msc = pyms3d.mscomplex()
    msc.compute_bin(data_file_name, dim)
    # simplify for base case
    msc.simplify_pers(thresh=0.0, is_nrm=True)
    # get the critical points
    cps_max = msc.cps(3)
    # cps_min, cps_1sad, cps_2sad = msc.cps(0), msc.cps(1), msc.cps(2)
    # value of scalar function at the critical point
    cps_fun_vals = msc.cps_func()

    # simplify for highest persistence (normalized value = 1)
    msc.simplify_pers(thresh=1, is_nrm=True)
    # get the critical point pairs that got cancelled (all got cancelled)
    # saddle--maximum pairs cancelled in simplification
    cp_pairs = msc.cps_pairid()

    # persistence diagram between birth and death
    sad = cp_pairs[cps_max, np.newaxis]
    b_value = cps_fun_vals[sad]
    d_value = cps_fun_vals[cps_max, np.newaxis]
    pers = d_value - b_value
    # p_diagm_list = np.concatenate((b_value, d_value, pers, sad,
    #                               cps_max[:, np.newaxis]), axis=1)

    print('Saddle values (quantiles - 0.25 spacing) ',
          np.quantile(b_value, [0, 0.25, 0.50, 0.75, 1.0]))
    print('Max values (quantiles - 0.25 spacing) ',
          np.quantile(d_value, [0, 0.25, 0.50, 0.75, 1.0]))

    plt.figure()
    plt.hist(b_value, stacked=True, label="saddle value", rwidth=0.1,
             cumulative=True, density=True, histtype='bar')
    plt.hist(d_value, stacked=True, label="Max value", rwidth=0.1,
             cumulative=True, density=True, histtype='bar')
    plt.xlabel('func value')
    plt.legend()
    plt.show()

    print("Take a note of the persistence value at the knee!")
    plt.figure()
    plt.plot(np.sort(pers, axis=0)[::-1], np.arange(pers.shape[0]))
    plt.xlabel("Persistence")
    plt.ylabel("Survived critical points")
    plt.title("Persistence curve")
    plt.show()

    plt.figure()
    plt.plot(b_value, d_value, 'r.')
    plt.plot([0, max(max(b_value), max(d_value))],
             [0, max(max(b_value), max(d_value))])
    plt.xlabel("Birth")
    plt.ylabel("Death")
    plt.title("Persistence diagram")
    plt.show()
    return None
