# SpanSeq Developer Guide

## Architecture

SpanSeq follows a layered architecture:

```
CLI (cli.py)
  │
  ▼
Config (config.py)          ← maps CLI args → dataclass
  │
  ▼
Pipeline (pipeline.py)      ← orchestration: init, run_split, run_reduce
  │
  ├── Distance (distance.py) ← distance matrix computation + Hobohm1 (mixin)
  │
  ├── Output (output.py)     ← merge tables, FASTA partitions, class columns (mixin)
  │
  ├── FASTA (fasta.py)       ← parse_fasta / write_fasta (no BioPython)
  │
  └── Applications (applications/)
        ├── kma.py           ← KmaIndexApp, KmaDistApp
        ├── mash.py          ← MashApp
        ├── cdhit.py         ← CdHitApp
        ├── ccphylo.py       ← CCPhyloDbscanApp, CCPhyloMakespanApp, CCPhyloTreeApp
        ├── ggsearch.py      ← GGSearchApp
        └── mmseqs2.py       ← MMseqs2SearchApp
```

### Key design decisions

- **No BioPython dependency.** FASTA parsing uses a lightweight custom parser in `fasta.py`.
- **cgecore v3 base classes.** All application wrappers extend `ApplicationRunner` or `AlignerRunner` from `cgecore`. These handle executable validation, command building, and subprocess execution.
- **DistanceMixin.** Distance methods are separated into `distance.py` as a mixin class that `SpanSeqPipeline` inherits from. This keeps the pipeline file focused on orchestration.
- **OutputMixin.** Output formatting (table merging, FASTA partitioning, class columns) lives in `output.py` as a second mixin.
- **Dataclass config.** `SpanSeqConfig` is a `@dataclass` with computed properties (`distance_tool`, `kma_dist_flag`, etc.) that translate user-facing options into internal parameters.

## Module Reference

### `fasta.py`

```python
FastaRecord = namedtuple("FastaRecord", ["id", "description", "seq"])

def parse_fasta(path) -> List[FastaRecord]
def write_fasta(handle, record: FastaRecord) -> None
```

- `id`: first whitespace-delimited token after `>`
- `description`: everything after `>` (includes `id`)
- `seq`: concatenated sequence lines

### `config.py`

`SpanSeqConfig` dataclass. Key properties:

| Property | Returns |
|---|---|
| `distance_tool` | `"kma"`, `"mash"`, `"ggsearch36"`, `"mmseqs2"`, or `"mmseqs-fast"` based on `distance_method` |
| `kma_dist_flag` | Integer flag for KMA dist (64, 256, 2048, 32) or `None` |
| `effective_dist_value` | `min_dist / kma_dist_factor` |
| `needs_hobohm` | `True` if approach is `hobohm_reduce` or `hobohm_split` |
| `sample_name` | Stem of the input file path |

The `DISTANCE_TOOLS` dict maps CLI distance method names to tool names:

```python
DISTANCE_TOOLS = {
    "jaccard": "kma", "szymkiewicz_simpson": "kma",
    "cosine": "kma", "kmer_inv": "kma",
    "mash": "mash", "identity": "ggsearch36",
    "mmseqs2": "mmseqs2", "mmseqs-fast": "mmseqs-fast",
}
```

### `pipeline.py`

`SpanSeqPipeline(DistanceMixin, OutputMixin)` — the main entry point. Methods:

| Method | Role |
|---|---|
| `__init__(config)` | Sets up directories, resolves tool paths, initializes app runners |
| `run()` | Dispatches to `_run_split()` or `_run_reduce()` |
| `_run_split()` | Full split pipeline: hobohm1 → distance → (tree) → dbscan → makespan → output |
| `_run_reduce()` | Reduce pipeline: KMA index+hobohm1 → makespan |

### `output.py`

`OutputMixin` — inherited by `SpanSeqPipeline`. Methods:

| Method | Description |
|---|---|
| `_merge_tables()` | Joins cluster and makespan TSVs into a partition table |
| `_create_fasta_partitions()` | Splits a FASTA file into per-partition files |
| `_add_class_columns()` | Merges class labels for imbalance-aware makespan |

### `distance.py`

`DistanceMixin` — inherited by `SpanSeqPipeline`. Methods:

| Method | Tool | Description |
|---|---|---|
| `_compute_distance()` | — | Dispatcher based on `config.distance_tool` |
| `_distance_kma()` | KMA | Index + distance matrix |
| `_distance_mash()` | Mash | Sketch + triangle distance |
| `_distance_ggsearch()` | GGSearch36 | Incremental all-vs-all alignment → PHYLIP matrix |
| `_distance_mmseqs2()` | MMseqs2 | Easy-search self-vs-self → PHYLIP matrix |
| `_cluster_mmseqs2_fast()` | MMseqs2 | Direct clustering via easy-cluster (skips distance + DBSCAN) |
| `_run_hobohm1()` | CD-HIT or KMA | Pre-reduction of redundant sequences |

### `applications/`

Each application module wraps an external tool:

| Module | Class | Base class | Purpose |
|---|---|---|---|
| `kma.py` | `KmaIndexApp` | `ApplicationRunner` | Build KMA index |
| `kma.py` | `KmaDistApp` | `ApplicationRunner` | Compute KMA distance matrix |
| `mash.py` | `MashApp` | `ApplicationRunner` | Mash sketch + triangle |
| `cdhit.py` | `CdHitApp` | `ApplicationRunner` | CD-HIT clustering |
| `ccphylo.py` | `CCPhyloDbscanApp` | `ApplicationRunner` | CCPhylo DBSCAN |
| `ccphylo.py` | `CCPhyloMakespanApp` | `ApplicationRunner` | CCPhylo makespan |
| `ccphylo.py` | `CCPhyloTreeApp` | `ApplicationRunner` | CCPhylo tree (Newick from distance matrix) |
| `ggsearch.py` | `GGSearchApp` | `AlignerRunner` | GGSearch36 global alignment |
| `mmseqs2.py` | `MMseqs2SearchApp` | `AlignerRunner` | MMseqs2 easy-search |
| `mmseqs2.py` | `MMseqs2ClusterApp` | `ApplicationRunner` | MMseqs2 easy-cluster |

**ApplicationRunner** apps implement:
- `build_command(...)` → `List[str]` — builds the subprocess command
- `map_outputs(workdir, ...)` → `Dict[str, Path]` — maps output file paths

**AlignerRunner** apps additionally implement:
- `build_command(query, db, out_prefix, ...)` — takes `SampleSpec` and `DatabaseSpec`
- `parse_result_file(path)` → `dict` — parses tool output into structured results

## Adding a New Distance Method

To add a new distance method (e.g. `diamond`):

### 1. Create the application wrapper

Create `src/spanseq/applications/diamond.py`:

```python
from cgecore.applications.base import AlignerRunner  # or ApplicationRunner

class DiamondApp(AlignerRunner):
    def build_command(self, query, db, out_prefix, **kwargs):
        cmd = [str(self.exec_path), "blastp", ...]
        return cmd

    def map_outputs(self, workdir, out_prefix=None, **kwargs):
        if out_prefix is None:
            return {}
        return {"result": Path(f"{out_prefix}.tsv")}

    @staticmethod
    def parse_result_file(path):
        # Parse the output file into a dict
        ...
```

### 2. Register in config.py

Add to the `DISTANCE_TOOLS` dict:

```python
DISTANCE_TOOLS = {
    ...
    "diamond": "diamond",
}
```

### 3. Add distance method in distance.py

Add a `_distance_diamond` method to `DistanceMixin`:

```python
def _distance_diamond(self, input_file, output_file):
    cfg = self.config
    # Run diamond and build PHYLIP distance matrix
    ...
    return output_file
```

Update `_compute_distance` dispatcher:

```python
elif tool == "diamond":
    return self._distance_diamond(input_file, output_file)
```

### 4. Initialize in pipeline.py

In `_init_apps`, add the tool initialization:

```python
elif tool == "diamond":
    diamond = self._resolve_tool(cfg.diamond_path, "diamond")
    self.diamond = DiamondApp(exec_path=diamond)
```

Add `diamond_path` to `SpanSeqConfig` and the CLI parser.

### 5. Add to CLI

In `cli.py`, add `"diamond"` to the `--distanceMethod` choices.

### 6. Write tests

Create `tests/applications/test_diamond.py` following the pattern in existing test files.

## Testing

### Running tests

```bash
# cgecore v3 is not yet on PyPI — install from the development branch first
pip install git+https://bitbucket.org/genomicepidemiology/cgecore.git@v3_dev

pip install -e ".[test]"
pytest -v
```

### Coverage

```bash
pytest --cov=spanseq --cov-report=term-missing
```

### Test structure

```
tests/
├── conftest.py              # shared fixtures: sample_fasta
├── test_fasta.py            # parse_fasta / write_fasta
├── test_config.py           # SpanSeqConfig dataclass
├── test_cli.py              # CLI argument parsing + main()
├── test_main.py             # __main__ module
├── test_pipeline.py         # pipeline orchestration (all external calls mocked)
├── applications/
│   ├── conftest.py          # fake_exec fixture
│   ├── test_kma.py
│   ├── test_mash.py
│   ├── test_cdhit.py
│   ├── test_ccphylo.py
│   ├── test_ggsearch.py
│   └── test_mmseqs2.py
```

### Mocking patterns

**Fake executable fixture** — needed because `ApplicationRunner` validates the executable path:

```python
@pytest.fixture
def fake_exec(tmp_path):
    exe = tmp_path / "tool"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    return exe
```

**Mocking tool resolution** — avoids needing real tools installed:

```python
mocker.patch.object(SpanSeqPipeline, "_resolve_tool", return_value=fake_exe)
```

**Mocking subprocess.run** — for methods that call subprocess directly:

```python
mocker.patch("subprocess.run")
```

**Mocking application methods** — for pipeline integration tests:

```python
mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)
```

**Important:** When mocking with `side_effect` lambdas, the parameter names must match the keyword argument names used at the call site:

```python
# If the call site uses: self._compute_distance(input_file=..., output_file=...)
# Then the lambda must use the same kwarg names:
side_effect=lambda input_file, output_file: ...
```

### cgecore specifics

- `DatabaseSpec` requires `sequencetype` to be a non-None string (e.g. `"genes"`).
- `SampleSpec` requires `type` (e.g. `"assembled"`) and `files` (list of Paths).

## Project Layout

```
SpanSeq/
├── pyproject.toml           # build config, dependencies, entry point
├── src/
│   └── spanseq/
│       ├── __init__.py      # version
│       ├── __main__.py      # python -m spanseq
│       ├── cli.py           # argparse + main()
│       ├── config.py        # SpanSeqConfig dataclass
│       ├── fasta.py         # FASTA parser/writer
│       ├── distance.py      # DistanceMixin
│       ├── output.py        # OutputMixin
│       ├── pipeline.py      # SpanSeqPipeline
│       └── applications/
│           ├── __init__.py
│           ├── kma.py
│           ├── mash.py
│           ├── cdhit.py
│           ├── ccphylo.py
│           ├── ggsearch.py
│           └── mmseqs2.py
├── tests/
│   └── ...
├── docs/
│   ├── user_guide.md
│   └── developer_guide.md
└── data/
    └── envs/
        └── spanseqenv.yml
```
