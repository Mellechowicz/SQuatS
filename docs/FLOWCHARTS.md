# EXSQS / SQuatS — program flowcharts

Nine flowcharts describing how the engine works, from the component map down to the
per-candidate evaluation. Every node group is mapped to its implementation
(file : line, function/class; everything lives in `namespace exsqs` unless noted).
Line numbers refer to v1.7.1. SPEC tags `[Axx]`/`[Dxx]` refer to `docs/SPEC.md`.

> This file is generated documentation and deliberately **untracked** (the curated
> history carries no meta files). View with any Mermaid-capable renderer
> (GitHub, VS Code + Mermaid extension, mermaid.live).

Contents

1. [Component map](#1-component-map)
2. [CLI dispatch and run lifecycle](#2-cli-dispatch-and-run-lifecycle)
3. [Config loading and run setup](#3-config-loading-and-run-setup)
4. [Serial driver: lockstep islands](#4-serial-driver-lockstep-islands)
5. [One generation: `IslandEngine::advance`](#5-one-generation-islandengineadvance)
6. [Candidate pipeline and evaluation](#6-candidate-pipeline-and-evaluation)
7. [MPI driver: rank-distributed islands](#7-mpi-driver-rank-distributed-islands)
8. [Checkpoint and resume](#8-checkpoint-and-resume)
9. [`score` subcommand](#9-score-subcommand)

---

## 1. Component map

Orientation diagram (not a control flow): the two binaries, the `exsqs_core`
library modules, and who calls whom. `evolution` is the hub; `score` reuses the
same evaluation stack on external structures.

```mermaid
flowchart TD
    subgraph binaries["binaries"]
        MAIN["exsqs<br/>src/main.cpp"]
        MPIMAIN["exsqs_mpi<br/>src/mpi_main.cpp"]
    end
    subgraph core["library exsqs_core (namespace exsqs)"]
        CONFIG["config<br/>YAML → RunConfig"]
        EVO["evolution<br/>RunContext, IslandEngine,<br/>drivers, outputs"]
        SCORE["score<br/>external-structure scoring"]
        STRUCT["structure + lattice + linalg<br/>Structure, supercell, POSCAR I/O"]
        ZONES["zones<br/>ZoneTable coordination shells"]
        CORR["correlation<br/>count_pairs, E_pure, E_floor, weights"]
        SYM["symmetry<br/>spglib wrapper: SymmetryInfo"]
        DISP["displacements<br/>phonopy port: D count"]
        DEDUP["dedup<br/>canonical labels, hashes"]
        RNG["rng<br/>CounterRng (stateless, counter-keyed)"]
        SER["serialize<br/>ByteWriter/Reader, trajectory_signature"]
    end
    SPGLIB(["spglib 2.7.0<br/>(FetchContent)"])
    YAML(["yaml-cpp 0.8.0<br/>(FetchContent)"])

    MAIN --> CONFIG
    MAIN --> EVO
    MAIN --> SCORE
    MPIMAIN --> CONFIG
    MPIMAIN --> EVO
    MPIMAIN --> SER
    EVO --> STRUCT
    EVO --> ZONES
    EVO --> CORR
    EVO --> SYM
    EVO --> DISP
    EVO --> DEDUP
    EVO --> RNG
    EVO --> SER
    SCORE --> EVO
    SCORE --> STRUCT
    SCORE --> CORR
    SCORE --> SYM
    SCORE --> DISP
    CONFIG --> YAML
    SYM --> SPGLIB
    ZONES --> STRUCT
    CORR --> ZONES
    DISP --> SYM
    DEDUP --> STRUCT
```

Where in the code

| Node | Implementation |
|---|---|
| `exsqs` binary | `src/main.cpp` (`main` at `main.cpp:23`) |
| `exsqs_mpi` binary | `src/mpi_main.cpp` (`main` at `mpi_main.cpp:101`); built only when MPI is found (`CMakeLists.txt:74`) |
| config | `include/exsqs/config.hpp:14` (`struct RunConfig`), `src/config.cpp:71` (`load_config`) |
| evolution | `include/exsqs/evolution.hpp`, `src/evolution.cpp` (1213 lines; engine class `IslandEngine` at `evolution.cpp:359`) |
| score | `include/exsqs/score.hpp`, `src/score.cpp` |
| structure/lattice/linalg | `src/structure.cpp`, `src/lattice.cpp`, `include/exsqs/linalg.hpp` |
| zones | `src/zones.cpp:9` (`build_zones`) |
| correlation | `src/correlation.cpp` |
| symmetry | `src/symmetry.cpp:36` (`get_symmetry`, wraps `spg_get_dataset`) |
| displacements | `src/displacements.cpp` (phonopy 4.3.1 port, `[A15]`) |
| dedup | `src/dedup.cpp` |
| rng | `include/exsqs/rng.hpp:37` (`class CounterRng`), purposes enum `rng.hpp:17` |
| serialize | `src/serialize.cpp`, `include/exsqs/serialize.hpp` |

---

## 2. CLI dispatch and run lifecycle

What happens between `exsqs …` on the command line and the exit code.
Exit codes: **0** converged (min E_pure ≤ e_tol), **3** budget exhausted
(resumable, normal), **1** error.

```mermaid
flowchart TD
    START(["exsqs argv"]) --> SC{"argv[1] == 'score'?"}
    SC -- yes --> SCORE["score subcommand<br/>(flowchart 9)"]
    SC -- no --> GE{"argv[1] == 'geom'?"}
    GE -- yes --> GEOM["load_config → RunContext::build<br/>→ write_poscar(undecorated supercell)"]
    GEOM --> EX0(["exit 0 / 1"])
    GE -- no --> FLAGS["parse flags:<br/>--set key.path=value (accumulate)<br/>--out DIR ⇒ output.dir=DIR<br/>--resume DIR<br/>positional: config.yaml"]
    FLAGS --> RESQ{"--resume given<br/>and no --out?"}
    RESQ -- "yes (v1.7.1 fix)" --> RESFIX["append override<br/>output.dir = resume DIR<br/>(outputs default to the state dir)"]
    RESQ -- no --> LOAD
    RESFIX --> LOAD["load_config(path, overrides)<br/>(flowchart 3)"]
    LOAD --> RFC["run_from_config(cfg, resume_dir + '/state.ckpt')<br/>banner → RunContext::build → run_evolution<br/>→ write_outputs → summary print"]
    RFC --> OK{"min E_pure ≤ e_tol<br/>on any island?"}
    OK -- yes --> EXITOK(["exit 0 (SUCCESS)"])
    OK -- no --> EXIT3(["exit 3 (budget exhausted,<br/>resumable)"])
    RFC -. "std::exception" .-> EXIT1(["exit 1 (error)"])
    LOAD -. "std::exception" .-> EXIT1
```

Where in the code

| Node | Implementation |
|---|---|
| dispatch `score` | `src/main.cpp:24-50` → `run_score_cli` (`score.cpp:109`) |
| dispatch `geom` | `src/main.cpp:51-81` → `RunContext::build` + `write_poscar` (`structure.cpp:120`) |
| flag parsing | `src/main.cpp:85-115` (`--out` rewritten to `output.dir=` at `main.cpp:102`) |
| resume-output fix | `src/main.cpp:120-132` (the v1.7.1 foot-gun fix: resume defaults outputs to the state dir) |
| `load_config` | `src/config.cpp:71` |
| `run_from_config` | `src/evolution.cpp:1155-1211`; returns `out.success ? 0 : 3` at `evolution.cpp:1210` |
| exception → 1 | `src/main.cpp:137-140` |
| usage/exit-code banner | `src/main.cpp:9-21` |

---

## 3. Config loading and run setup

Two stages: `load_config` turns YAML + `--set` overrides into a validated
`RunConfig`; `RunContext::build` derives everything the run needs exactly once
(geometry, shells, weights, floor, symmetry groups).

```mermaid
flowchart TD
    subgraph LC["load_config — src/config.cpp:71"]
        Y["parse YAML (yaml-cpp)"] --> OV["apply --set overrides<br/>apply_override()"]
        OV --> LAT["lattice: sc|bcc|fcc|hcp (a, c)<br/>or POSCAR prototype file"]
        LAT --> SUP["supercell: diag [n1,n2,n3]<br/>or 3×3 integer matrix H, det(H) > 0"]
        SUP --> COMP["composition: fractions →<br/>integer counts (largest remainder)"]
        COMP --> CT{"each |x_achieved − x_target|<br/>≤ comp_tol = 0.02<br/>and every count ≥ 1?"}
        CT -- no --> FAIL(["throw → exit 1"])
        CT -- yes --> REST["zones / error(weights, mode, gamma) /<br/>evolution / symmetry / displacements /<br/>parallel / rng.seed / output"]
    end
    REST --> RCB
    subgraph RC["RunContext::build — src/evolution.cpp:198"]
        RCB["make_supercell(proto, H)<br/>undecorated geometry, N sites"] --> BZ["build_zones: shell radii +<br/>per-site neighbour lists<br/>(27-image background [A1])"]
        BZ --> MW["make_weights: w(n) =<br/>inv_n | inv_n_pow | inv_R | custom"]
        MW --> EF["E_floor: L1 quantization bound<br/>(diagonal or full mode)"]
        EF --> SP["site_permutations Π of the<br/>empty lattice [A9] (dedup group)"]
        SP --> SPS["seed_perms: permutations of<br/>non-identity ops (for constructive<br/>seeding [D4])"]
    end
    SPS --> ET{"e_tol < 0 (auto)?"}
    ET -- yes --> ET3["effective e_tol = 3 × E_floor"]
    ET -- no --> ETC["effective e_tol = cfg.e_tol"]
```

Where in the code

| Node | Implementation |
|---|---|
| override mechanics | `src/config.cpp:35` (`apply_override`), applied at `config.cpp:83` |
| lattice / supercell / composition sections | `src/config.cpp:89` / `:114` / `:131` |
| largest-remainder rounding | `src/config.cpp:50` (`largest_remainder`), counts at `:142` |
| comp_tol gate (fixed 2%) | `src/config.cpp:144-149` |
| zones/error/evolution/… sections | `src/config.cpp:153` / `:159` / `:192` / `:275` / `:283` / `:288` / `:311` / `:312` |
| `RunContext::build` | `src/evolution.cpp:198-217`; `struct RunContext` at `evolution.hpp:70` |
| `make_supercell` | `src/structure.cpp:27` |
| `build_zones` | `src/zones.cpp:9-80` (half-cell-width warning `[A3]` at `zones.cpp:78`) |
| `make_weights` | `src/correlation.cpp:8` |
| `e_floor_diagonal` / `e_floor_full` | `src/correlation.cpp:82` / `:95` |
| `site_permutations`, `permutation_of_op` | `src/dedup.cpp:55` / `:51` |
| `effective_e_tol` | `src/evolution.cpp:219-221` |

---

## 4. Serial driver: lockstep islands

`run_evolution` advances all islands one generation per round (lockstep) with
nested OpenMP: *outer* threads across islands, *inner* threads inside each
island's evaluation rounds. Thread counts change wall time only, never results
(counter-keyed RNG `[A14]`).

```mermaid
flowchart TD
    RFC["run_from_config: banner,<br/>E_floor / D(P1)=6N reference,<br/>checkpoint callback installed"] --> ENG["construct one EngineHandle<br/>per island (0 … islands−1)"]
    ENG --> THR["thread budget:<br/>outer = min(islands, T)<br/>inner = T / outer"]
    THR --> RES{"--resume?"}
    RES -- yes --> LRS["load_run_state(state.ckpt)<br/>(flowchart 8) → round"]
    RES -- no --> SEED["parallel over islands:<br/>seed_generation0()<br/>(flowchart 6, gen-0 seeding)"]
    LRS --> LOOP
    SEED --> LOOP{"any island<br/>not done?"}
    LOOP -- no --> FINAL["final save_run_state<br/>(enables resume chains [A13])"]
    LOOP -- yes --> ADV["parallel over islands (outer threads):<br/>engine.advance() — one generation<br/>(flowchart 5)"]
    ADV --> MIG{"round % migration_every == 0<br/>and islands > 1 and migrants > 0?"}
    MIG -- yes --> RING["synchronous ring migration [SPEC 8.1]:<br/>each active island sends its best-k pool<br/>copies to the next active island"]
    MIG -- no --> CKPT
    RING --> CKPT{"round %<br/>checkpoint_every == 0?"}
    CKPT -- yes --> SAVE["save_run_state<br/>→ outdir/state.ckpt"]
    CKPT -- no --> LOOP
    SAVE --> LOOP
    FINAL --> MERGE["merge_island_results [D8]:<br/>pool all islands, sort by<br/>(E_obj, E_pure, hash),<br/>dedup by canonical label, keep top-M"]
    MERGE --> WO["write_outputs:<br/>generations.csv,<br/>best_00.vasp … best_MM.vasp,<br/>summary.json (final) /<br/>checkpoint.json (mid-run)"]
    WO --> RC(["return success ? 0 : 3"])
```

Where in the code

| Node | Implementation |
|---|---|
| `run_from_config` (banner, callback, exit) | `src/evolution.cpp:1155-1211` (checkpoint callback lambda `:1183-1190`) |
| `run_evolution` | `src/evolution.cpp:954-1042` |
| engine construction | `src/evolution.cpp:957-959`; `class EngineHandle` (pimpl) `island_engine.hpp:19`, impl `evolution.cpp:820-846` |
| thread budget outer/inner | `src/evolution.cpp:963-969` |
| resume vs parallel seeding | `src/evolution.cpp:971-988` |
| lockstep round loop | `src/evolution.cpp:995-1033` |
| ring migration | `src/evolution.cpp:1016-1030` (`emigrants` `evolution.cpp:625`, `receive_migrants` nearby) |
| periodic + final state save | `src/evolution.cpp:1031-1034` |
| `merge_island_results` | `src/evolution.cpp:925-952` |
| `write_outputs` | `src/evolution.cpp:1044-1153` (`generations.csv` `:1049`, `best_%02zu.vasp` `:1075`, `summary.json`/`checkpoint.json` `:1083`) |
| single-island convenience path | `evolve_island`, `src/evolution.cpp:916-923` |

---

## 5. One generation: `IslandEngine::advance`

The v1.1 loop body: termination checks, extinction `[A11]` with elitism `[D5]`,
repopulation `[A12]` with a three-stage fallback ladder, stagnation stop `[A13]`.

```mermaid
flowchart TD
    A["advance() — no-op if done"] --> T1{"min E_pure ≤ e_tol?"}
    T1 -- yes --> S1(["stop 'e_tol'<br/>success = true"])
    T1 -- no --> T2{"gen ≥ max_generations?"}
    T2 -- yes --> S2(["stop 'max_generations'"])
    T2 -- no --> T3{"wall time > max_wall_s?"}
    T3 -- yes --> S3(["stop 'wall_time'"])
    T3 -- no --> SORT["rank population by<br/>(E_obj, E_pure, hash, index);<br/>mark top elitism_best as elite [D5]"]
    SORT --> BETA{"survival mode?"}
    BETA -- ratio --> PRAT["keep-prob p = min(1, e_min / e_i)"]
    BETA -- metropolis --> SCHED["β schedule [A11]:<br/>geometric β·growth^gen if schedule=1;<br/>auto β = ln2 / median(e_i − e_min) if β ≤ 0"]
    SCHED --> PMET["keep-prob p = min(1, exp(−β(e_i − e_min)))"]
    PRAT --> DRAW["per individual: elite ⇒ keep;<br/>else keep iff u < p,<br/>u from CounterRng(seed, island, gen, i,<br/>ExtinctionDraw)"]
    PMET --> DRAW
    DRAW --> KILL["killed = P − survivors"]
    KILL --> REPOP["repopulate killed slots [A12] —<br/>per-slot fallback ladder:<br/>stage 0: pick parent, mutate_sigma [D6]<br/>(≤ retry_budget attempts)<br/>stage 1: constructive seed<br/>(≤ retry_budget attempts)<br/>stage 2: plain shuffle<br/>(≤ 10×retry_budget, then throw)"]
    REPOP --> PIPE["per round: candidates through the<br/>canonical → dedup → evaluate → P1-commit<br/>pipeline (flowchart 6);<br/>P1-rejected slots resume next round"]
    PIPE --> INC["++gen; record GenStats;<br/>every checkpoint_every gens:<br/>checkpoint callback (mutex-guarded)"]
    INC --> STAG{"killed > 0 and<br/>fallback seeds ≥ killed<br/>for stagnation_stop<br/>consecutive gens? [A13]"}
    STAG -- yes --> S4(["stop 'stagnation'"])
    STAG -- no --> DONE(["generation complete;<br/>driver calls advance() again"])
```

Where in the code

| Node | Implementation |
|---|---|
| `IslandEngine::advance` | `src/evolution.cpp:436-622` (class at `:359`) |
| termination checks | `src/evolution.cpp:439-453` |
| ranking + elitism | `src/evolution.cpp:455-468` |
| β schedule + auto-β | `src/evolution.cpp:469-478` |
| survival draws (ratio `:490`, metropolis `:492`) | `src/evolution.cpp:482-503`; RNG purpose `ExtinctionDraw` `rng.hpp:17` |
| slot ladder (`SlotState`, stages 0/1/2) | `src/evolution.cpp:505-601` (`struct SlotState` `:516`; stage transitions `:546`, `:560`; exhaustion throw `:569-572`) |
| `mutate_sigma` `[D6]` | `src/evolution.cpp:267-346` (Poisson swap count `poisson_draw` `:254`) |
| constructive fallback | `constructive_from_rng`, `src/evolution.cpp:132` |
| record + checkpoint callback | `src/evolution.cpp:607-612` (`record` `:760`) |
| stagnation stop | `src/evolution.cpp:613-621` |

---

## 6. Candidate pipeline and evaluation

Every candidate — gen-0 seed or child — passes the same four-step pipeline.
`fill_canonical` and `evaluate_individual` are parallel-safe; `dup_or_insert`
and `commit_p1` run serially in slot order, which keeps results bit-identical
for any thread count.

```mermaid
flowchart TD
    subgraph SRC["σ sources"]
        SEED0["gen-0 seeding [D4] (mixed default):<br/>odd slots: constructive —<br/>invariant under a non-identity op<br/>⇒ guaranteed non-P1;<br/>even slots: rejection —<br/>shuffled composition vector"]
        MUT["mutation [D6] (mutate_sigma):<br/>k = 1 + Poisson(λ) swaps;<br/>pick ONE nontrivial parent op g,<br/>swap species between equal-size<br/>unlike-species cycles of ⟨g⟩<br/>(child keeps ⟨g⟩ ⇒ stays non-P1);<br/>fallback: plain unlike-species<br/>pair swap [A12]"]
    end
    SEED0 --> P1
    MUT --> P1
    subgraph PIPE["pipeline (per candidate)"]
        P1["1. fill_canonical:<br/>canonical_labels = lex-min<br/>relabelling over dedup group Π [A10];<br/>FNV-1a hash"]
        P1 --> P2{"2. dup_or_insert:<br/>canonical key already<br/>in archive? [A9]"}
        P2 -- "yes → duplicate" --> REJ1(["reject, next attempt"])
        P2 -- "no → insert" --> P3
        P3["3. evaluate_individual"] --> P4{"4. commit_p1 [A7][A8]:<br/>space group == 1 and<br/>P1 count ≥ p1_elite_quota?"}
        P4 -- yes --> REJ2(["reject (stays archived),<br/>slot resumes next round"])
        P4 -- no --> ACC(["accepted into population;<br/>pool_insert keeps island's<br/>best-M by (E_obj, E_pure, hash)"])
    end
    subgraph EVAL["evaluate_individual — one spglib call"]
        E1["decorate: species onto<br/>the fixed geometry"] --> E2["get_symmetry (spglib):<br/>SG number/symbol,<br/>equivalent atoms, ops;<br/>keep nontrivial ops as stab_ops"]
        E2 --> E3["D = displacement_count:<br/>phonopy 4.3.1 least-displacements<br/>port [A15] over independent atoms"]
        E3 --> E4["count_pairs over ZoneTable shells<br/>→ Π(t,t2 | shell n)"]
        E4 --> E5["E_pure = Σ w(n) · |Π − x| (L1)<br/>diagonal or full-pairs mode"]
        E5 --> E6(["E_obj = E_pure · D^γ  [D3]"])
    end
    P3 -.-> E1
```

Where in the code

| Node | Implementation |
|---|---|
| gen-0 seeding round loop | `IslandEngine::seed_generation0`, `src/evolution.cpp:379-432` (mixed odd/even split `:395`; exhaustion cap `:410-413`) |
| rejection / constructive seed σ | `src/evolution.cpp:230` / `:240` (+ `constructive_from_rng` `:132`) |
| `mutate_sigma` `[D6]` | `src/evolution.cpp:267-346` (cycle decomposition of ⟨g⟩ `:280-306`; cycle swap `:307-334`; plain-swap fallback `:335-343`) |
| `fill_canonical` | `src/evolution.cpp:711-714`; `canonical_labels` `dedup.cpp:67`, `hash_labels` `dedup.cpp:82` |
| `dup_or_insert` (archive) | `src/evolution.cpp:716-724` |
| `evaluate_individual` | `src/evolution.cpp:726-740` (`E_obj` formula at `:739`) |
| `commit_p1` (P1 quota) | `src/evolution.cpp:742-749` |
| `pool_insert` | `src/evolution.cpp:751-758` |
| `decorate` | `src/structure.cpp:95` |
| `get_symmetry` | `src/symmetry.cpp:36-77` (atomic numbers `:16`, `pointgroup_order` `:98`) |
| `displacement_count` | `src/displacements.cpp:89/:93` (`least_displacements` `:72/:76`; internals `disp_one` `:17`, `disp_two` `:34`, `minus_needed` `:53`) |
| `count_pairs` | `src/correlation.cpp:36` (`CorrData` `correlation.hpp:17`) |
| `e_pure_diagonal` / `e_pure_full` | `src/correlation.cpp:53` / `:66` |
| `Individual` (σ, canonical, sg, stab_ops, D, E_pure, E_obj) | `include/exsqs/evolution.hpp:18` |

---

## 7. MPI driver: rank-distributed islands

`exsqs_mpi` lifts the same lockstep schedule across processes. Islands are
owned round-robin (`island i → rank i mod R`); migration and checkpointing go
through collective byte-blob exchanges. Results are **bit-identical** to the
serial binary for any rank count (`[A14]`, gate T-MPI1).

```mermaid
flowchart TD
    INIT["MPI_Init; rank, size"] --> PARSE["parse flags (as serial);<br/>only rank 0 prints"]
    PARSE --> SETUP["ensure_little_endian →<br/>load_config (rank ≠ 0: logging off) →<br/>RunContext::build"]
    SETUP --> OWN["ownership: island i owned by<br/>rank i mod size;<br/>ALL ranks construct all EngineHandles,<br/>each advances only its own"]
    OWN --> THR["per-rank threads:<br/>outer = min(#owned, T), inner = T/outer"]
    THR --> RES{"--resume?"}
    RES -- yes --> LOADALL["every rank reads the shared<br/>state.ckpt, keeps owned islands"]
    RES -- no --> SEEDOWN["parallel seed_generation0()<br/>over owned islands"]
    LOADALL --> LOOP
    SEEDOWN --> LOOP["lockstep round:<br/>done-flags → MPI_Allreduce(MAX)"]
    LOOP --> ANY{"any island<br/>still active?"}
    ANY -- no --> FSAVE["final checkpoint gather<br/>(resubmission chains [A13])"]
    ANY -- yes --> ADV["parallel advance() over<br/>owned islands (flowchart 5)"]
    ADV --> MIG{"round % migration_every == 0?"}
    MIG -- yes --> XCHG["active flags → Allreduce;<br/>serialize emigrants as<br/>(island id, blob) records →<br/>MPI_Allgatherv to all ranks;<br/>ring: dst = next active island,<br/>owner rank deserializes + receives"]
    MIG -- no --> CK
    XCHG --> CK{"round %<br/>checkpoint_every == 0?"}
    CK -- yes --> GATH["owned engine blobs →<br/>MPI_Gatherv to rank 0 →<br/>rank 0 writes state.ckpt → Barrier"]
    CK -- no --> LOOP
    GATH --> LOOP
    FSAVE --> RESULTS["IslandResult blobs →<br/>Gatherv to rank 0;<br/>rank 0: merge_island_results →<br/>write_outputs → code = 0 or 3"]
    RESULTS --> BCAST["MPI_Bcast(code) → MPI_Finalize<br/>→ exit code on every rank"]
    SETUP -. "exception" .-> ABORT(["MPI_Abort(…, 1)"])
```

Where in the code

| Node | Implementation |
|---|---|
| `main` | `src/mpi_main.cpp:101-302` |
| flag parsing | `src/mpi_main.cpp:109-134` |
| setup + logging silencing | `src/mpi_main.cpp:136-141` (`cfg.log_info = false` for rank ≠ 0 at `:139`) |
| ownership round-robin | `src/mpi_main.cpp:143-149` |
| thread budget | `src/mpi_main.cpp:151-157` (`local_threads` `:41`) |
| resume (all ranks read shared file) | `src/mpi_main.cpp:167-170` → `load_run_state` (`evolution.cpp:883`) |
| seeding of owned islands | `src/mpi_main.cpp:171-185` |
| done-flag `MPI_Allreduce(MAX)` | `src/mpi_main.cpp:203-210` (`:207`) |
| parallel advance | `src/mpi_main.cpp:212-225` |
| migration via `MPI_Allgatherv` | `src/mpi_main.cpp:227-260` (helper `allgather_bytes` `:51`; record codec `parse_records` `:88`; individuals codec `put/get_individuals`, `serialize.cpp`) |
| checkpoint gather to rank 0 | `src/mpi_main.cpp:187-200` (`save_state_root` lambda; `gather_bytes_root` `:69`; `save_run_state_blobs` `evolution.cpp:848`), called at `:261`, final `:263` |
| result gather + outputs + Bcast | `src/mpi_main.cpp:265-296` |
| exception → `MPI_Abort` | `src/mpi_main.cpp:297-300` |

---

## 8. Checkpoint and resume

`state.ckpt` is a little-endian binary snapshot of every island engine (RNG
counters, population, archive, pool, generation), guarded by a **trajectory
signature** that rejects any config change except raising budget caps.

```mermaid
flowchart TD
    subgraph SAVE["save — evolution.cpp:848/:872"]
        SV1["serialize each EngineHandle<br/>→ per-island blob"] --> SV2["ByteWriter: magic 'EXSQSTAT',<br/>format version u32 = 1,<br/>trajectory_signature(cfg),<br/>island count, (id, blob) records"]
        SV2 --> SV3["atomic write:<br/>state.ckpt.tmp → fs::rename"]
    end
    subgraph LOADCK["load — evolution.cpp:883"]
        LD1["read file → ByteReader"] --> LD2{"magic == 'EXSQSTAT'?<br/>version == 1?"}
        LD2 -- no --> LERR(["throw → exit 1"])
        LD2 -- yes --> LD3{"trajectory_signature<br/>matches current config?"}
        LD3 -- no --> LERR2(["throw: trajectory signature<br/>mismatch (only budget caps<br/>may change on resume)"])
        LD3 -- yes --> LD4{"island count ==<br/>cfg.islands?"}
        LD4 -- no --> LERR
        LD4 -- yes --> LD5["deserialize each engine,<br/>reopen_for_resume()"]
        LD5 --> LD6(["resume round =<br/>max generation across islands"])
    end
    SIG["trajectory_signature hashes every<br/>trajectory-relevant field (lattice, H,<br/>composition, zones, weights, γ, survival,<br/>mutation, seeding, symprec, seed, islands,<br/>migration) and deliberately EXCLUDES<br/>max_generations, max_wall_s,<br/>output.*, threads"]
    SIG -.-> SV2
    SIG -.-> LD3
```

Where in the code

| Node | Implementation |
|---|---|
| `save_run_state_blobs` (format + atomic write) | `src/evolution.cpp:848-870` |
| `save_run_state` (engines → blobs) | `src/evolution.cpp:872-881` |
| `load_run_state` (all checks + resume round) | `src/evolution.cpp:883-914` |
| engine (de)serialization | `EngineHandle::serialize/deserialize`, `src/evolution.cpp:844-845`; `reopen_for_resume` via handle |
| `trajectory_signature` | `src/serialize.cpp:148` (exclusion rationale documented at `serialize.hpp:99`) |
| byte codecs (`ByteWriter`/`ByteReader`, individuals, results) | `include/exsqs/serialize.hpp:19-87`, `src/serialize.cpp:29-146` |
| CLI wiring (`--resume DIR` → `DIR/state.ckpt`) | `src/main.cpp:103-108`, `:135-136`; MPI: `mpi_main.cpp:167-170` |

---

## 9. `score` subcommand

Scores externally produced structures (POSCAR) on the config geometry with the
exact evaluation stack of the search — E_pure, E/E_floor, D, SG, E_obj.

```mermaid
flowchart TD
    CLI(["exsqs score config.yaml<br/>[--set …] POSCAR… [--json PATH]"]) --> LC["load_config"]
    LC --> CTX["RunContext::build;<br/>print header: E_floor,<br/>reference D(P1) = 6N"]
    CTX --> EACH["for each POSCAR file"]
    EACH --> RP["read_poscar"]
    RP --> MAP["sigma_on_geometry:<br/>1. lattice matches config geometry<br/>elementwise (orientation-strict)?<br/>2. same site count?<br/>3. species ⊆ config composition?<br/>4. each geometry site ↔ exactly one<br/>input site (wrapped frac dist ≤ 1e-4)?<br/>5. composition == config counts [A5]?"]
    MAP -- "any check fails" --> FERR["print error for this file;<br/>all_ok = false"]
    MAP -- ok --> EV["score_structure —<br/>same stack as the search:<br/>decorate → get_symmetry →<br/>displacement_count D →<br/>count_pairs → E_pure →<br/>E_obj = E_pure · D^γ"]
    EV --> ROW["print table row:<br/>E_pure, E/floor, D, D(P1)/D,<br/>SG (symbol), E_obj"]
    ROW --> EACH
    FERR --> EACH
    EACH -- "all files done" --> JS{"--json given?"}
    JS -- yes --> JSON["write JSON: e_floor, d_p1,<br/>scores[] per file"]
    JS -- no --> RET
    JSON --> RET(["exit 0 iff every<br/>file scored, else 1"])
```

Where in the code

| Node | Implementation |
|---|---|
| CLI dispatch | `src/main.cpp:24-50` |
| `run_score_cli` (header, loop, table, JSON, exit) | `src/score.cpp:109-166` |
| `sigma_on_geometry` (all five checks) | `src/score.cpp:33-91` (lattice tol `:40`, site match tol `:61`, composition check `:84-89`) |
| `score_structure` | `src/score.cpp:93-107` |
| `read_poscar` | `src/structure.cpp:149` |
| evaluation stack | same functions as flowchart 6: `symmetry.cpp:36`, `displacements.cpp:93`, `correlation.cpp:36/:53/:66` |
| `ScoreResult` | `include/exsqs/score.hpp:18` |
