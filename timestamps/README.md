# Timestamps

[OpenTimestamps](https://opentimestamps.org) proofs for the assets of release `v3`.
Each proof commits to the SHA-256 of the file below.

| file | SHA-256 |
|---|---|
| `li2-minus-half-irrationality-v3.pdf` | `6f111306eed20d3abc65230ce67604a5fe7cc3e1588d6f0d4a2a8c33acdfc8c2` |
| `li2-minus-half-irrationality-v3.tar.gz` | `40410d634aa1d5b761a2df59956fd3b382782444dc65313b1a2dc007ab5f4ee7` |

To verify, download the release asset into this directory and run

    pip install opentimestamps-client
    ots verify timestamps/<file>.ots
