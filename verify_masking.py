import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from elmo import get_simulation


def run_simulation(classname, challenges):
    print("Getting Simulation!)")
    Simulation = get_simulation(classname, repository='elmo/projects')
    print("Got Simulation!")
    simulation = Simulation()
    simulation.set_challenges(challenges)
    print("Running Simulation!")
    res = simulation.run()
    print("Simulation completed!")
    assert not res['error'], res['error']
    return simulation


def printed_coefficients(pd):
    """Rebuild the printed uint16 coefficients from the per-byte output.

    print2bytes writes the two bytes of each uint16 in big-endian order
    (high byte first), so a coefficient is pd[2k] * 256 + pd[2k + 1].
    """
    pd = np.asarray(pd, dtype=np.int64)
    return pd[:, 0::2] * 256 + pd[:, 1::2]


def hw(x):
    """Hamming weight of a 16-bit value."""
    x = int(x) & 0xFFFF
    c = 0
    while x:
        c += x & 1
        x >>= 1
    return c


def max_corr_vs_secret(traces, hws):
    """Max |Pearson correlation| between any trace column and HW(secret).

    Columns with (near-)zero variance carry no information and are skipped,
    otherwise a constant column would produce a spurious huge correlation.
    """
    n, npoints = traces.shape
    t = traces - traces.mean(axis=0)
    t_std = t.std(axis=0)
    valid = t_std > 1e-12
    worst = 0.0
    for c in range(hws.shape[1]):
        h = hws[:, c]
        h_c = h - h.mean()
        denom = n * t_std * h_c.std()
        with np.errstate(invalid='ignore', divide='ignore'):
            corr = np.nan_to_num((t * h_c[:, None]).sum(axis=0) / denom)
        corr[~valid] = 0.0
        m = np.abs(corr).max()
        if m > worst:
            worst = m
    return worst


def main():
    np.random.seed(0)
    # Enough traces to push the max|corr| noise floor well below the 0.5
    # threshold: for a masked trace, max|corr| ~ sqrt(2*ln(P)/N) (pure noise),
    # while the unmasked leak stays at ~1.0.
    nb_challenges = 256
    challenges = [np.random.choice(
        [-2, -1, 0, 1, 2],
        (2, 256),
        p=[1/16, 4/16, 6/16, 4/16, 1/16],
    ) for _ in range(nb_challenges)]

    # NOTE: the two simulations share the same ELMO output directory, so every
    # result (printed data, traces) must be read right after the run that
    # produced it, before the other simulation overwrites it.
    print('Running unmasked KyberNTT...')
    sim_unmasked = run_simulation('KyberNTTSimulation', challenges)
    pd_unmasked = np.array(sim_unmasked.get_printed_data())
    traces_unmasked = sim_unmasked.get_traces()

    print('Running masked KyberNTT...')
    sim_masked = run_simulation('KyberNTTMaskedSimulation', challenges)
    pd_masked = np.array(sim_masked.get_printed_data())
    traces_masked = sim_masked.get_traces()

    # --- Correctness: both versions must compute the same NTT(s) mod q ---
    cu = printed_coefficients(pd_unmasked)
    cm = printed_coefficients(pd_masked)
    assert cu.shape == cm.shape, (cu.shape, cm.shape)
    # barrett_reduce maps zero-residue coefficients to either 0 or 3329
    # depending on the sign of the input, so compare modulo q.
    n_diff = np.count_nonzero((cu - cm) % 3329)
    print('Correctness: coefficients differ in {} of {}'.format(n_diff, cu.size))
    assert n_diff == 0, 'Masked and unmasked NTT results differ!'

    print('Correctness check passed!')

    # --- Defense: CPA-style correlation of traces with HW(s) ---
    challenges_flat = np.array([ch.ravel() for ch in challenges])
    hws = np.array([[hw(v) for v in row] for row in challenges_flat])

    corr_unmasked = max_corr_vs_secret(traces_unmasked, hws)
    corr_masked = max_corr_vs_secret(traces_masked, hws)
    print('Unmasked max |corr| vs HW(s): {:.3f}'.format(corr_unmasked))
    print('Masked   max |corr| vs HW(s): {:.3f}'.format(corr_masked))

    assert corr_unmasked > 0.4, 'Unmasked version should leak the secret!'
    assert corr_masked < corr_unmasked * 0.5, 'Masked version still leaks!'
    assert corr_masked < 0.5, 'Masked max correlation too high!'

    # --- Plots: the two traces must differ ---
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    axes[0].plot(np.arange(traces_unmasked.shape[1]), traces_unmasked[0], color='C0')
    axes[0].set_title('Unmasked KyberNTT - power trace (trace 0)')
    axes[1].plot(np.arange(traces_masked.shape[1]), traces_masked[0], color='C1')
    axes[1].set_title('Masked KyberNTT - power trace (trace 0)')
    for ax in axes:
        ax.set_xlabel('instruction index')
        ax.set_ylabel('power')
    fig.tight_layout()
    fig.savefig('traces_comparison.png', dpi=120)

    print('Saved plots to traces_comparison.png')
    print('All checks passed!')


if __name__ == '__main__':
    main()
