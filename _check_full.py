import numpy as np
from elmo import get_simulation

Sim = get_simulation('KyberFullSimulation', repository='elmo/projects')
sim = Sim()
sim.set_challenges(sim.get_random_challenges(3))
res = sim.run()
assert not res['error'], res['error']

traces = sim.get_traces()
print('traces shape:', traces.shape)
n_win = sim.N_TRACES_PER_CHALLENGE
assert traces.shape[0] == 3 * n_win, traces.shape

n = sim.get_number_of_challenges()
assert n == 3, n
for ci in range(n):
    d = sim.get_printed_data_per_challenge(ci)
    n_labels = len(d['labels'])
    total = n_labels * sim.LABEL_BYTES + sim.DUMP_BYTES
    print('challenge {}: {} labels, per-challenge bytes = {}'.format(ci, n_labels, total))
    assert n_labels == sim.N_WINDOWS
    assert total == sim.PER_CHALLENGE_BYTES
    assert d['ss_enc'].tobytes() == d['ss_dec'].tobytes(), 'ss mismatch in challenge {}'.format(ci)
    assert len(d['pk']) == 800 and len(d['sk']) == 1632 and len(d['ct']) == 768

print('ss_enc == ss_dec for all challenges: OK')
print('ALL FULL-SIM CHECKS PASSED')