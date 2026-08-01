import numpy as np

from elmo import get_simulation


def run_simulation(classname, challenges):
    Simulation = get_simulation(classname, repository='elmo/projects')
    simulation = Simulation()
    simulation.set_challenges(challenges)
    res = simulation.run()
    assert not res['error'], res['error']
    return simulation


def main():
    np.random.seed(0)
    nb_challenges = 128
    challenges = [np.random.choice(
        [-2, -1, 0, 1, 2],
        (2, 256),
        p=[1/16, 4/16, 6/16, 4/16, 1/16],
    ) for _ in range(nb_challenges)]

    print('Running unmasked KyberNTT...')
    sim_unmasked = run_simulation('KyberNTTSimulation', challenges)
    print('Running masked KyberNTT...')
    sim_masked = run_simulation('KyberNTTMaskedSimulation', challenges)

    # --- Correctness: both versions must compute the same NTT(s) ---
    pd_unmasked = np.array(sim_unmasked.get_printed_data())
    pd_masked = np.array(sim_masked.get_printed_data())
    assert pd_unmasked.shape == pd_masked.shape, (pd_unmasked.shape, pd_masked.shape)
    # barrett_reduce maps zero-residue coefficients to 0 (non-negative inputs)
    # or 3329 (negative multiples of q), so compare modulo q.
    n_diff = np.count_nonzero((pd_unmasked - pd_masked) % 3329)
    print('Correctness: printed values differ in {} of {}'.format(n_diff, pd_unmasked.size))
    assert n_diff == 0, 'Masked and unmasked NTT results differ!'

    print('Correctness check passed!')


if __name__ == '__main__':
    main()
