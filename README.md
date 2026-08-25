# DNA Storage Noise-Guessing Decoder (MATLAB)

Proof-of-concept MATLAB implementation of a self-correcting algorithm for DNA data storage systems, based on a noise-guessing (GRAND-inspired) paradigm. The channel error model is derived from real Oxford Nanopore sequencing measurements published in Lopez et al. (_Nature Communications_, 2019).

---

## Overview

DNA storage introduces three types of errors during sequencing:

| Error type     | Symbol | Probability (measured) |
| -------------- | :----: | ---------------------: |
| Deletion       |   D    |                 2.95 % |
| Insertion      |   I    |                 2.16 % |
| Substitution   |   S    |                 1.75 % |
| Any impairment |   —    |                 6.87 % |

Binary payload is encoded into nucleotide sequences using the mapping `00→A, 01→C, 10→G, 11→T`. A CRC checksum is appended to each codeword. After the simulated channel corrupts the codeword, the decoder tries to recover the original sequence by guessing the most probable error patterns first, verifying each guess with a CRC check.

---

## Repository structure

```
Algorithm.m           – Main simulation entry point
RollingAlg.m          – DNA channel error model (deletions / insertions / substitutions)
SelfCorrectAlg.m      – Noise-guessing self-correcting decoder
CRCSender.m           – CRC encoding (sender side)
CRCDetector.m         – CRC verification (receiver side)
BinomialDis.m         – Binomial distribution for multi-error threshold filtering
BinomialDisSingle.m   – Binomial distribution for single-error probability
calcErrors.m          – Priority scoring for combined error-type classes (t1,t2,t3)
calcErrorPerm.m       – Unique nucleotide combination generator
calcSubComb.m         – Ordered substitution candidate generator (Bayes-ranked)
SubstitutionNuc.m     – Per-nucleotide substitution sampler
GenProbRandi.m        – Uniform random integer helper
MyError.m             – Class: single-type error object (prob, count, type)
MyErrors.m            – Class: combined error object (prob, errDel, errIns, errSub)
combProbs.m           – Class: nucleotide combination with probability
untitled.m            – Complexity measurement plot (k=32, 1/2/3 errors)
```

---

## Algorithm summary

### 1. Encoding (`CRCSender`)

- Random binary payload of length `N` bits is generated.
- Every two bits are mapped to one nucleotide → codeword of length `N/2` nucleotides.
- CRC remainder is appended. Polynomial used: `x^8 + x^7 + x^4 + x^3 + x + 1`.

### 2. Channel simulation (`RollingAlg`)

- Each nucleotide is independently subjected to deletion, insertion, or substitution with the probabilities measured from Nanopore data.
- Substitution targets are drawn using per-nucleotide conditional probabilities (e.g. A→G most likely at 70.9 %).
- Substitutions are applied in a second pass to maximise randomness of error positions.

### 3. Decoding (`SelfCorrectAlg`)

The decoder is a priority-guided exhaustive search:

1. **Outer loop** — `decisionTree`: binomial distribution over total error count. Classes above threshold `λ` are queued, ordered by probability. This decides _how many_ errors to attempt to correct.

2. **Single-error branch** — tries deletion, insertion and substitution (in order of decreasing per-type probability). Uses length-check (`mod`) to skip physically impossible error types before iterating positions.

3. **Multi-error branch** — `calcErrors` scores every $(t_1,t_2,t_3)$ combination by the ranking function
   $$s(t_1,t_2,t_3) = p_d^{\,t_1}\,p_i^{\,t_2}\,p_s^{\,t_3}$$
   (proportional to the probability of one specific pattern of that type). Classes are tried highest-score first. For each class, all position combinations (`nchoosek`) and nucleotide combinations are tested.

4. Each candidate codeword is verified by CRC. On success the corrected word is returned immediately.

---

## Parameters

| Parameter | Default | Description                                                                     |
| --------- | :-----: | ------------------------------------------------------------------------------- |
| `n`       |   24    | Codeword length in nucleotides                                                  |
| `N`       |   32    | Payload length in bits                                                          |
| `thresh`  |  0.10   | Minimum binomial probability to include an error count in the guessing list (λ) |
| `count`   |  1250   | Number of simulation repetitions                                                |
| `probD`   | 0.0295  | Deletion probability                                                            |
| `probI`   | 0.0216  | Insertion probability                                                           |
| `probS`   | 0.0175  | Substitution probability                                                        |

---

## Key results (from associated thesis measurements)

### BER/BLER vs CRC degree (k = 16, ~1 Mbit payload)

| CRC    | λ = 20 %               | λ = 10 %                | λ = 5 %                  |
| ------ | ---------------------- | ----------------------- | ------------------------ |
| CRC-6  | BER 0.032 / BLER 0.194 | BER 0.0285 / BLER 0.17  | BER 0.0282 / BLER 0.17   |
| CRC-8  | BER 0.033 / BLER 0.199 | BER 0.022 / BLER 0.1325 | BER 0.022 / BLER 0.1325  |
| CRC-10 | BER 0.034 / BLER 0.218 | BER 0.017 / BLER 0.1    | BER 0.017 / BLER 0.104   |
| CRC-16 | BER 0.019 / BLER 0.124 | BER 0.0195 / BLER 0.124 | BER 0.0141 / BLER 0.0935 |

- **CRC-10 at λ = 10 %** is the best practical trade-off: BLER = 0.1 (90 % correct blocks), 3× faster than CRC-16.

### Decoder complexity (k = 16, CRC-6/10/16)

| Errors corrected |      Time |  CRC checks |
| :--------------: | --------: | ----------: |
|        1         |  6 m 15 s |     629 026 |
|        2         | 57 m 23 s |   6 341 280 |
|        3         |   1 605 m | 142 230 810 |

Complexity grows roughly 10× per additional error for time and combinations.

---

## Requirements

- MATLAB R2020b or newer
- Communications Toolbox (for `comm.CRCGenerator` / `comm.CRCDetector`)

---

## References

1. R. Lopez et al., "DNA assembly for nanopore data storage readout," _Nature Communications_ 10, art. 2933 (2019). https://doi.org/10.1038/s41467-019-10978-4
2. K. R. Duffy, J. Li and M. Médard, "Capacity-Achieving Guessing Random Additive Noise Decoding," _IEEE Trans. Inf. Theory_, vol. 65, no. 7, pp. 4023–4040, July 2019.
3. MathWorks, `comm.CRCGenerator` documentation. https://www.mathworks.com/help/comm/ref/comm.crcgenerator-system-object.html

---

## License

See [LICENSE](../LICENSE) in the root of this repository.
