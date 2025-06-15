import graph_tool.all as gt #from Tiago Peixoto
import numpy as np
from numba import jit,njit # to do just in time compilation at machine code
from concurrent.futures import ProcessPoolExecutor  # to do multiprocessing across graphs realizations
import pickle,glob
from dataclasses import dataclass
from typing   import Dict, Any, List
import networkx as nx #to generate Watts-Strogatz faster
import time
import datetime #to print the runtime of intermediate chunks

def graph_to_csr(G):
    """
    Given an undirected graph-tool Graph G,
    returns (offs, nbrs) in CSR form, such that the neighbors of node i are nbrs[offs[i]:offs[i+1]]
    """
    degs = np.array( G.get_total_degrees(G.get_vertices()),    dtype=np.uint16)
    offs = np.r_[0, np.cumsum(degs, dtype=np.uint32)]          # (0, k_0, k_0+k_1, ... , k_0+...+k_{N-2}, k_0+...+k_{N-1})
    nbrs = np.empty(offs[-1], dtype=np.uint16)

    aux_offs = np.delete(offs,-1)
    for v in G.vertices():
        i = int(v)
        edgelist = G.get_out_edges(v)        
        targets   = edgelist[:, 1]
        
        for u in targets:
            nbrs[aux_offs[i]] = int(u)
            aux_offs[i] += 1
    
    return nbrs,offs

def log_schedule(start, stop, factor=1.5):
    """
    this is just to create a list of times at wich we sample the order parameter
    """
    out = [start]
    val = start
    while True:
        val = int(round(val * factor))
        if val > stop: break
        out.append(val)
    return out

def bits_needed(N): #to exctract a sample in 0,N
    b = 1
    while (1 << b) < N:
        b += 1
    return b
@njit
def xorshift128p_next_float(seed_array):
    """
    In-place update of the 2-element seed_array [s0, s1] using the Xorshift128+ algorithm.
    Returns a random float64 in [0,1).
    """
    s1 = seed_array[0]
    s0 = seed_array[1]
    
    # Rotate
    seed_array[0] = s0
    
    # Xorshift steps
    s1 ^= (s1 << np.uint64(23))
    s1 ^= (s1 >> np.uint64(17))
    s1 ^= s0
    s1 ^= (s0 >> np.uint64(26))
    
    # Write back to seed_array
    seed_array[1] = s1
    
    # Output is s0 + s1
    out = s0 + s1  # still a uint64

    # Convert to float in [0,1).
    # Mask out 53 bits (the mantissa of double) for uniform distribution:
    return (out & np.uint64(0x1FFFFFFFFFFFFF)) / float(1 << 53)
@njit
def xorshift128p_next_uint(seed): #for the choiche of nodes to try to update

    """
    In‐place update of the 2‐element seed [s0, s1].
    Returns a uint64 = s0 + s1 (the next state).
    """
    s1 = seed[0]
    s0 = seed[1]
    seed[0] = s0
    # the core xorshift:
    s1 ^= (s1 << np.uint64(23))
    s1 ^= (s1 >> np.uint64(17))
    s1 ^= s0
    s1 ^= (s0 >> np.uint64(26))
    seed[1] = s1
    return s0 + s1

@njit(inline='always')
def uint64_to_unit_float(u): #for the flip trial
    # 53 random mantissa bits → float in [0,1)
    return (u & np.uint64(0x1FFFFFFFFFFFFF)) * (1.0 / (1 << 53))

@njit
def init_sigma_xorshift(N, seed):  #for the initialization of the configuration
    sigma = np.empty(N, np.int8)
    for i in range(N):
        # take the least-significant bit of a 64-bit PRNG step
        sigma[i] = 1 if (xorshift128p_next_uint(seed) & 1) else -1
    return sigma

@njit 
def ρ(σ,E): #This is the order parameter
    count = 0
    for e in E:
        count+= (σ[e[0]] !=σ[e[1]])

    return count / len(E)

@njit 
def do_1_step(N, σ, seed, nn,cs_k, inv_deg, a):

    for _ in range(N):
        u = xorshift128p_next_uint(seed)    # ONE RNG call for the choiche of nodes and for the update trial
        i = int(u % N)                       # there are some avoidable correlations 
                                              # but not a tragedy for the scope of the project

        s_nb = 0
        for p in range(cs_k[i], cs_k[i+1]):
            s_nb += σ[nn[p]]  # Σ_{j ∈ nn(i)} σ_j       with csr (no slicing needed)

        p_flip = a + (1 - 2*a) * (0.5 * (1.0 - σ[i] * s_nb * inv_deg[i]))

        if p_flip == 0.0: #can be useful to speed up the non-noise part
            continue
        if p_flip == 1.0 or uint64_to_unit_float(u) < p_flip:  # use low bits for probability (see RNG)
            σ[i] = -σ[i]

@njit 
def do_T_steps_ρ_sample (N,σ, ρ_t, sample_times , seed, E,nn,cs_k, inv_deg, a):

    for t in range (0,sample_times[0]):  #first i sample at each timestep

        ρ_t[t] = ρ(σ,E)

        do_1_step (N,σ, seed, nn,cs_k, inv_deg, a)
    
    sample_idx = 0                       # as time passes i switch to a log sampling (see auxiliary fuction)
    for t in range (sample_times[0],sample_times[-1]+1):

        if t == sample_times[sample_idx]:
            ρ_t[sample_times[0]+sample_idx]  = ρ(σ,E)
            sample_idx+=1

        do_1_step (N,σ, seed, nn,cs_k, inv_deg, a)

@njit (fastmath=True)
def simulate_Voter_dynamics(N,E, nn,cs_k, inv_deg, num_realiz, sample_times, seed, a):

    """
    I define a single array ρ(t) and for each realization the code writes over the previous one, 
    but I keep track of the mean and variance (over realization of the dynamics) using a Welford update.
    All quantities finish with _t to signal that they are time dependent trajectories, ρ_t[t] is a number.
    """
    sample_size = sample_times[0]+len(sample_times)

    ρ_t     = np.zeros(sample_size, dtype = np.float64)
    ρ_avg_t = np.zeros(sample_size, dtype = np.float64)
    ρ_Var_t = np.zeros(sample_size, dtype = np.float64)

    for r in range (0,num_realiz):

        σ = init_sigma_xorshift(len(cs_k)-1, seed)

        do_T_steps_ρ_sample (N,σ, ρ_t, sample_times, seed, E,nn,cs_k, inv_deg, a)

        Δρ_t     = ρ_t - ρ_avg_t
        ρ_avg_t += Δρ_t / (r+1)
        ρ_Var_t  += Δρ_t * (ρ_t - ρ_avg_t ) # at the end of the loop this is num_realiz * Var
    
    ρ_Var_t = ρ_Var_t/num_realiz

    return ρ_avg_t,ρ_Var_t

@dataclass
class GraphConfig: #instanciated many times inside a particular experiment 

    N:           int
    model:       str          # “BA”, “ER”, “CM”, etc.
    params:      Dict[str,Any]   #  {"m":4}, or {"p":0.05}, or {"gamma":2.5}
    graph_seed:  int

@dataclass
class DynConfig:   #instanciated one time in main

    num_realiz:   int
    sample_times: np.ndarray
    seed0:        np.ndarray
    noise_values: list[float]        

from numpy.random import default_rng
from numpy.random import SeedSequence

def Klemm_Eguiluz_ssf(N, m, rng, a=None): # SSF stands for structured scale free 
    #undirected version of the process explained in Phys. Rev. E 65, 036123 (2002) (chapter 3)
    if a==None: a=m

    G = gt.complete_graph(m)                             # the process starts with an m-clique
    active_nodes = list(G.get_vertices())                # all starting nodes are active

    for t in range(m,N):

        G.add_vertex()                                   # a node gets added    

        new_edges = [(t, j) for j in active_nodes]       # it forms a link with the active nodes

        
        G.add_edge_list(new_edges)

        active_nodes.append(t)                           # it gets activated

        deactiv_weights = 1 / (a+G.get_total_degrees(active_nodes))
        deactiv_probs = deactiv_weights/deactiv_weights.sum()

        dactivated_node = rng.choice(active_nodes,p=deactiv_probs)   ##SEEEEEEEEEED
        active_nodes.remove(dactivated_node)             # one active node gets deactivated with some P

    return G

def cut_ring (G,N,m):  #circular network with a cut to compare with Klemm Eguiluz wich is not closed

    for i in range(m):
        for k in range(m - i):
            e = G.edge(i, N - k - 1)   
            if e is not None: G.remove_edge(e)

    return G


def make_graph(cfg,μ): #takes a GraphConfig instance and returns a Graph object
    """
    Valid models: {"BA", "SSF", "SWSF", "String", "SW", "EN", "SBM"}
    """
    gt.seed_rng(cfg.graph_seed + μ)   # reproducibility

    p     = cfg.params            # shorthand
    model = cfg.model.upper()

    if model == "BA":                                     # Barabasi Albert
        m = p["m"]
        return gt.price_network(cfg.N, m=m, directed=False, gamma=1,c=0) #put c=0 to get the original BA

    elif model == "SSF":                                  
        m = p.get("m", 3)

        seed = cfg.graph_seed + μ
        rng  = default_rng(seed)

        return Klemm_Eguiluz_ssf(cfg.N, m=m, rng=rng)

    elif model == "SWSF":                                 # Small world scale free
        m = p.get("m", 3);  prob = p["p"]                 # for p = 1 you get random scale free

        seed = cfg.graph_seed + μ
        rng  = default_rng(seed)

        G = Klemm_Eguiluz_ssf(cfg.N, m=m, rng=rng)
        E = G.num_edges() 
        n_iter = int(-0.5 * E * np.log(1-prob)) 
        #correct number of swap trials (at least for SSF rewire) such that prob is a true rewire_p

                                            # preserves degree distribution but not correlations
        n_failed = gt.random_rewire(G, model="configuration", edge_sweep=False, n_iter=n_iter) #modifies g in place
        return G

    elif model == "STRING":
        m = p.get("m", 3)

        G = gt.circular_graph(cfg.N, m)

        return cut_ring (G,cfg.N,m)

    elif model == "SW":                                   # Wats-Strogatz cut (with high p the cut is useless)
        m = p.get("m", 3);  prob = p["p"]                 # for p = 1 you get random network
        
        seed = cfg.graph_seed + μ
        H = nx.watts_strogatz_graph(n=cfg.N, k=2*m, p=prob, seed=seed)
        #      ↑ k = numero totale di vicini (m a destra + m a sinistra)

        # 2) converti in graph-tool
        G = gt.Graph(directed=False)
        G.add_vertex(cfg.N)
        G.add_edge_list(list(H.edges()))

        return cut_ring(G,cfg.N,m)

    elif model == "EN":                                    # exponential network, no bias
        m = p.get("m", 3)
        return gt.price_network(cfg.N, m=m, directed=False, gamma=0)  

    elif model == "SBM":
        sizes = p["sizes"]             # list of block sizes
        Pmat  = p["P"]                 # 2-D list / np.array of probs
        G, _  = gt.random_graph(cfg.N, model="blockmodel",
                                bg_probs=Pmat, sizes=sizes)
        return G

    else:
        raise ValueError("Unknown model "+model)
    

def worker(params):
    """
    params is a 2-tuple: (graph_cfg, dyn_cfg, mu_index) 
    everything inside here must be picklable/well‐defined at top level.
    """
    graph_cfg, dyn_cfg, μ = params

    G = make_graph(graph_cfg,μ) 

    # extract CSR + inv_degrees + Edge list
    nn,cs_k = graph_to_csr(G) # see auxiliary functions
    inv_deg = (1.0/G.degree_property_map("total").a).astype(np.float32)
    E       = np.array(G.get_edges(), dtype=np.uint32)

    # 2) for each noise value in dyn_cfg.noises
    out = []
    for idx,a in enumerate(dyn_cfg.noise_values):

        seed_sim = dyn_cfg.seed0.copy()
        # mix in the graph index
        seed_sim[0] ^= np.uint64(μ)
        seed_sim[1] ^= np.uint64(μ << 32)
        # mix in the noise‐level index
        seed_sim[0] ^= np.uint64(idx << 16)
        seed_sim[1] ^= np.uint64(idx << 48)

        # run the Numba‐JIT’d dynamics
        ρ_avg_t, ρ_Var_t = simulate_Voter_dynamics(
            graph_cfg.N,
            E, nn, cs_k, inv_deg,
            dyn_cfg.num_realiz,
            dyn_cfg.sample_times,
            seed_sim,
            a                       
        )
        
        params = dict(graph_cfg.params)     #params here is just m for BA, p for ER ecc..(could be multiple)
        params["noise"] = a                 # we add noise info

        out.append( (graph_cfg.model, params, μ, ρ_avg_t, ρ_Var_t) )
        
    return out            # list-of-tuples for this graph

class VoterExperiments:
    """
    A thin wrapper that:
      • builds GraphConfig lists for each experiment type
      • optionally chunks them
      • runs them in a fixed-size ProcessPool
      • checkpoints each chunk to disk if needed
    """
    def __init__(self, dyn_cfg: DynConfig, num_graphs: int = 12):
        self.dyn_cfg    = dyn_cfg
        self.num_graphs = num_graphs

    # ---------- internal helpers ---------------------------------

    def _pool_run(self, graph_cfgs):
        tasks = [(cfg, self.dyn_cfg, μ) 
                for cfg in graph_cfgs  for μ in range(self.num_graphs)]
        
        all_results = []
        with ProcessPoolExecutor(max_workers=self.num_graphs) as ex:
            for chunk in ex.map(worker, tasks):
                all_results.extend(chunk)   # flatten the inner lists
                
        return all_results

    def _run_with_checkpoint(self, graph_cfgs: List[GraphConfig], *,
                             tag: str,
                             checkpoint: bool,
                             num_chunks: int = 4,
                             out_pattern: str = "{tag}_chunk_{i}.pkl"):
        """
          * if checkpoint==False  → one shot
          * else                  → split cfgs into chunks, run+pickle each
        """
        if not checkpoint:
            return self._pool_run(graph_cfgs)

        all_results = []
        chunks = np.array_split(graph_cfgs, num_chunks)

        start_time = time.time()

        for i, sub_cfgs in enumerate(chunks):
            sub_cfgs = list(sub_cfgs)              # numpy -> python list

            elapsed   = time.time() - start_time
            elapsed_str = str(datetime.timedelta(seconds=int(elapsed)))
            print(f"[{tag}] chunk {i+1}/{num_chunks} : {len(sub_cfgs)} sweep_param values"
                                                            f" (elapsed {elapsed_str})")
            
            res = self._pool_run(sub_cfgs)

            fname = out_pattern.format(tag=tag, i=i)
            with open(fname, "wb") as f:
                pickle.dump(res, f)
            print(f"   → saved to {fname}")

            all_results.extend(res)
            del res                                    # free RAM
        return all_results

    # ---------- public sweep methods, add experiments here! ------

    def BA_sweep_m(self, N, m_values, *, base_seed, checkpoint=False, num_chunks=4,tag=None):
        """Fix N, sweep BA parameter m."""
        cfgs = [GraphConfig(N=N, model="BA",
                            params={"m": int(m)},
                            graph_seed=base_seed + int(m))
                for m in m_values]
        
        return self._run_with_checkpoint(cfgs, tag=tag, checkpoint=checkpoint, num_chunks=num_chunks)

    def BA_sweep_N(self, m, N_values, *, base_seed,  checkpoint=False, num_chunks=4):
        """Fix m, sweep graph size N."""
        cfgs = [GraphConfig(N=int(N), model="BA",
                            params={"m": int(m)},
                            graph_seed=base_seed + i)
                for i, N in enumerate(N_values)]
        
        return self._run_with_checkpoint(cfgs, tag="BA_N_m{m}", checkpoint=checkpoint, num_chunks=num_chunks)
    
    def SSF_sweep_m(self, N, m_values, *, base_seed, checkpoint=False, num_chunks=4,tag=None):
        """Fix N, sweep Klemm_Eguiluz parameter m."""

        cfgs = [GraphConfig(N=N, model="SSF",
                            params={"m": int(m)},
                            graph_seed=base_seed + int(m))
                for m in m_values]
        
        return self._run_with_checkpoint(cfgs, tag=tag, checkpoint=checkpoint, num_chunks=num_chunks)
    
    def SWSF_sweep_p(self, N, p_values, *, base_seed, checkpoint=False, num_chunks=4,tag=None):
        """Fix N, fix m=3(basta non fare nulla),sweep on rewiring p."""

        cfgs = [GraphConfig(N=N, model="SWSF",
                            params={"p": p},
                            graph_seed=base_seed + int(p*200000))
                for p in p_values]
        
        return self._run_with_checkpoint(cfgs, tag=tag, checkpoint=checkpoint, num_chunks=num_chunks)
    
    def SW_sweep_p(self, N, p_values, *, base_seed, checkpoint=False, num_chunks=4,tag=None):
        """Fix N, fix m=3(basta non fare nulla),sweep on rewiring p."""

        cfgs = [GraphConfig(N=N, model="SW",
                            params={"p": p},
                            graph_seed=base_seed + int(p*200000))
                for p in p_values]
        
        return self._run_with_checkpoint(cfgs, tag=tag, checkpoint=checkpoint, num_chunks=num_chunks)
    
if __name__ == "__main__":
    dyn_cfg = DynConfig(
        num_realiz     = 50,
        sample_times   = log_schedule(200, 10000, 1.1),
        seed0          = np.array([167, 390], np.uint64),        #used only for the dynamics
        noise_values   = [1E-3]
    )

    experiment = VoterExperiments(dyn_cfg, num_graphs=8)

    results = experiment.SWSF_sweep_p(
        N=200,
        p_values=[0,0.01],
        base_seed=12025,                                         #used only for the structure
        checkpoint=False,       # this line makes it so the trajectories are saved and divided in chunks 
        num_chunks=1,
        tag = 'SWSF_m3_p0_01_noisy_N200'         
    )
#the following is to store single trajectories
@njit
def do_1_step_and_update_upcount(N, σ, up, seed, nn, cs_k, inv_deg, a):
    """
    One lattice sweep *and* incremental update of the running
    number-of-up-spins counter `up`.  Returns the new up count.
    """
    for _ in range(N):
        u = xorshift128p_next_uint(seed)
        i = int(u % N)

        s_nb = 0
        for p in range(cs_k[i], cs_k[i+1]):
            s_nb += σ[nn[p]]

        p_flip = a + (1 - 2*a) * (0.5 * (1.0 - σ[i] * s_nb * inv_deg[i]))

        # test for a flip
        if p_flip != 0.0 and (p_flip == 1.0 or uint64_to_unit_float(u) < p_flip):
            σ[i] = -σ[i]
            # update the running counter in O(1)
            up += 1 if σ[i] == 1 else -1
    return up


@njit
def simulate_single_run_uptrace(N, T, seed, nn, cs_k, inv_deg, a):
    """
    One realisation, T full sweeps, sampling at every sweep.
    Returns a uint16 array of length T with the up-spin counts.
    """
    σ = init_sigma_xorshift(len(cs_k) - 1, seed)

    up_counts = np.empty(T, dtype=np.uint16)

    # initial count
    up = 0
    for s in σ:
        if s == 1:
            up += 1
    up_counts[0] = up

    # main loop
    for t in range(1, T):
        up = do_1_step_and_update_upcount(
            N, σ, up, seed, nn, cs_k, inv_deg, a
        )
        up_counts[t] = up

    return up_counts

def run_and_store_number_of_upspins_trajectory(a, T, graph_seed,seed, tag):
    graph_cfg = GraphConfig(N=4000, model="BA", params={"m": 2}, graph_seed=graph_seed)
    G = make_graph(graph_cfg,0)

    # extract CSR + inv_degrees + Edge list
    nn,cs_k = graph_to_csr(G) # see auxiliary functions
    inv_deg = (1.0/G.degree_property_map("total").a).astype(np.float32)

    trace = simulate_single_run_uptrace(graph_cfg.N, T, seed, nn, cs_k, inv_deg, a)
    np.save(f"uptrace_{tag}.npy", trace)

#run_and_store_number_of_upspins_trajectory(1E-3, 100000, 12252,np.array([157, 360]), "higher_noise")
#run_and_store_number_of_upspins_trajectory(2E-5, 100000, 12252,np.array([157, 360]), "lower_noise")