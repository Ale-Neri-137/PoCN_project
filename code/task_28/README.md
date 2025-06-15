


This code uses the `graph-tool` library, which is **not** published on PyPI, so `pip install graph-tool` will fail.  
Pick **one** of the methods below; Conda is by far the smoothest.

| Method | OS | Command(s) | Notes |
|--------|----|------------|-------|
| **Conda (recommended)** | Any platform with [Miniconda / Mamba](https://docs.conda.io/en/latest/miniconda.html) | ```bash\nconda install -c conda-forge graph-tool\n``` | Gives the latest stable build and handles the heavy C++/Boost stack automatically. |
| **APT repositories** | Ubuntu ≥ 22.04 / Debian Bookworm | ```bash\nsudo apt-get update\nsudo apt-get install python3-graph-tool\n``` | Versions track the distro’s release cycle and can lag the upstream project. |
| **Homebrew** | macOS (Apple Silicon & Intel) | ```bash\nbrew install graph-tool\n``` |  Homebrew sometimes removes the formula when it breaks—if that happens, fall back to Conda. |

Or dowload it from source with Linux/macOS

