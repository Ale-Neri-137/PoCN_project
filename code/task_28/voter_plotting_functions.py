def plot_p_sweep_and_BA(resultsSWSF, resultsSW, resultsBA,
                        sample_times, p_values, *,
                        dense0=None, cmap_name="viridis",
                        figsize=(14, 6)):
    # ---------------------------------------------------------------
    dense0 = dense0 or sample_times[0]
    times  = np.concatenate((np.arange(dense0), sample_times))

    cmap   = plt.get_cmap(cmap_name, len(p_values))
    colors = [cmap(i) for i in range(len(p_values))]

    # ---------------------------------------------------------------
    def _plot_one_dataset(ax, results):
        bucket = {p: [] for p in p_values}
        for _, params, μ, avg_t, _ in results:
            bucket[params["p"]].append(avg_t)

        handles = []
        for i, p in enumerate(p_values):
            runs = np.vstack(bucket[p])
            ax.plot(times, runs.T,
                    color=colors[i], alpha=0., lw=0.8, zorder=1)
            h, = ax.plot(times, runs.mean(0),
                         color=colors[i], alpha=0.9, lw=2,
                         label=f"p = {p:.4f}", zorder=30-i)   # top layer
            handles.append(h)
        return handles
    # ---------------------------------------------------------------

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    h_p = _plot_one_dataset(axes[0], resultsSWSF)   # SWSF
    _   = _plot_one_dataset(axes[1], resultsSW)     # SW

    # ---------------------------------------------------------------
    m_values  = [1, 2, 3]
    data_BA   = {m: [] for m in m_values}
    for _, params, μ, avg_t, _ in resultsBA:
        data_BA[params["m"]].append(avg_t)

    ba_colors = ['#215584', '#24768b', '#44948f', '#6db290']
    h_ba = []
    for i, m in enumerate(m_values):
        if m == 2:          # omitted per your note
            continue
        runs = np.vstack(data_BA[m])
        for j, ax in enumerate(axes):
            ax.plot(times, runs.T,
                    color=ba_colors[i], alpha=0., lw=0.8, zorder=1)
            if m==1:
                line, = ax.plot(times, runs.mean(0),
                                color=ba_colors[i], alpha=0.9, lw=2,
                                label=f"m = {m} (BA)", zorder=100,linestyle = "--")  # under p-curves
            if m==3:
                                line, = ax.plot(times, runs.mean(0),
                                color='#f4a261', alpha=0.9, lw=2,
                                label=f"m = {m} (BA)", zorder=2)  # under p-curves

            if j == 0:      # store handle once for the legend
                h_ba.append(line)

    # ---------------------------------------------------------------
    for ax, title in zip(axes, ["SWSF", "SW"]):
        ax.set_xlim(1, 1e5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time $t$ (sweeps)")
        if title == "SWSF":
            ax.set_ylabel(r"$\langle\rho(t)\rangle$")
        ax.set_title(title)
        ax.grid(False)

    # Legend inside first panel, vertical
    axes[0].legend(h_p + h_ba,
                   [h.get_label() for h in h_p + h_ba],
                   loc="lower left", bbox_to_anchor=(0.05, 0.05),
                   frameon=True, ncol=1)

    fig.tight_layout()
    plt.show()


def plot_noise_sweep_by_p(results,
                          sample_times,
                          p_values,
                          noise_values=None,
                          *,
                          dense0=None,
                          cmap_name="viridis",
                          figsize=(10, 4),      
                          sharey=True):
    """
    Disegna ⟨ρ(t)⟩ vs t per (p, noise) con:
      – p ∈ {0, 0.01}
      – legenda verticale nel primo pannello
      – etichette noise come potenze di dieci
    Resto dei parametri identico all’originale.
    """

    # ------------------------------------------------------------------
    # 0) quali p disegnare
    plotted_p = [p for p in p_values if p in (0, 0.01)]
    if not plotted_p:
        raise ValueError("Nessuno degli 'p_values' è 0 o 0.01.")
    
    # ------------------------------------------------------------------
    # 1) full time axis
    dense0 = dense0 or sample_times[0]
    times  = np.concatenate((np.arange(dense0), sample_times))

    # ------------------------------------------------------------------
    # 2) bucket the trajectories  →  data[p][noise] = [runs]
    if noise_values is None:
        noise_values = sorted({params['noise'] for _, params, _, _, _ in results})
    # escludi il noise bandito
    noise_values = [a for a in noise_values if not np.isclose(a, 1e-7)]
    data = {p: {a: [] for a in noise_values} for p in plotted_p}

    for _, params, _, avg_t, _ in results:
        p     = params['p']
        noise = params['noise']
        if p in data and noise in data[p]:
            data[p][noise].append(avg_t)

    # ------------------------------------------------------------------
    # 3) colori per i diversi noise
    cmap    = plt.get_cmap(cmap_name, len(noise_values))
    colours = [cmap(i) for i in range(len(noise_values))]

    # ------------------------------------------------------------------
    # 4) subplot grid fissa: 1×2
    n_p   = len(plotted_p)              # dovrebbe essere 2
    ncols = n_p
    nrows = 1
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=figsize,
                             sharey=sharey, sharex=True,
                             squeeze=False)
    axes = axes.ravel()

    # ------------------------------------------------------------------
    # funzione per formattare le etichette noise
    def fmt_noise(a):
        if a!=0:
            exp = int(np.floor(np.log10(a)))
            coeff = a / 10**exp
            if np.isclose(coeff, 1):
                return rf"$a = 10^{{{exp}}}$"
            else:
                return rf"$a = {coeff:g}\times10^{{{exp}}}$"
        else: return rf"$a = 0$"

    # ------------------------------------------------------------------
    for k, (ax, p) in enumerate(zip(axes, plotted_p)):
        for j, noise in enumerate(noise_values):
            runs = data[p][noise]
            if not runs:
                continue

            mean_run = np.vstack(runs).mean(axis=0)
            ax.plot(times, mean_run,
                    lw=1.5, color=colours[j],
                    label=fmt_noise(noise), alpha=0.8)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"$p = {p}$")
        ax.set_xlim(1, 1E5)
        ax.set_ylim(1E-4, 0.9)
        ax.grid(False)

        # legenda solo sul primo subplot
        if k == 0:
            h, l = ax.get_legend_handles_labels()
            ax.legend(h, l, title=r"noise $a$",
                      loc="lower left",
                      frameon=True,
                      ncol=1)          # verticale

    # etichette assi ----------------------------------------------------
    for ax in axes:                       # entrambi, perché sharex=True
        ax.set_xlabel(r"time $t$ (sweeps)")
    axes[0].set_ylabel(r"$\langle \rho(t)\rangle$")

    plt.tight_layout()
    plt.show()
