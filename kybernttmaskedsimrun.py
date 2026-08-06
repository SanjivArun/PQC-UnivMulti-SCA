from elmo import get_simulation

Sim = get_simulation('KyberNTTMaskedSimulation', repository='elmo/projects')
sim = Sim()
sim.set_challenges(sim.get_random_challenges(10))
sim.run()
print(sim.get_traces().shape)          # (10, nb_instructions) — note: longer than unmasked
print(sim.get_printed_data())