from elmo import get_simulation
Sim = get_simulation("KyberFullSimulation", repository="elmo/projects")
sim = Sim(); sim.set_challenges(sim.get_random_challenges(1)); r=sim.run()
import numpy as np; print("SHAPE", sim.get_traces().shape); print("NTRACES", r.get("nb_traces"), "ERR", r.get("error"))
