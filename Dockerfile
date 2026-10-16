FROM julia:1.11

WORKDIR /usr/src/app

# Copy environment definitions first
COPY Project.toml Manifest.toml* ./
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy all your local .jl scripts into the container
COPY . .

# Execute the simulation runner when the container starts
CMD ["julia", "--project=.", "run_selectedE.jl"]

