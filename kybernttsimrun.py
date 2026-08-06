from elmo import get_simulation

KyberNTTSimulation = get_simulation('KyberNTTSimulation')

simulation = KyberNTTSimulation()
challenges = simulation.get_random_challenges(10)
simulation.set_challenges(challenges)

simulation.run()
traces = simulation.get_traces()
print('traces shape:', traces.shape)
