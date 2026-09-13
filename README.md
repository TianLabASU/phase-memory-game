# Memory Under Growth

A self-contained, responsive interactive research figure based on Zhang et al., “Phase separation to buffer growth-mediated dilution in synthetic circuits,” *Cell* (2025). The trajectories are generated from the publication's stochastic model and parameters.

## Run locally

Open `index.html` directly, or serve the folder:

```bash
python3 -m http.server 8000 --directory phase-memory-game
```

Then visit `http://localhost:8000`.

## Website deployment

The project has no external dependencies. Upload the contents of this folder to any static host, including GitHub Pages, Netlify, Vercel, or an existing lab website.

## Simulation provenance

`data/generate_publication_data.py` is a seeded Python translation of the Figure 1 MATLAB simulation and its dependencies in the authors' public repository. It exports JSON plus a directly loadable browser data file, allowing `index.html` to work without a web server.

The animation communicates a reproducible stochastic simulation, not an experimental measurement or population average.
