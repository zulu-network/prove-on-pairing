# prove-on-pairing




## Here to run
* run poc
```bash
cargo test test_dual_miller_loop_with_c_wi_fixed -- --nocapture
```

* run groth16's verifier leverage power of 
```bash
cargo test test_groth16_verifier_with_c_wi -- --nocapture
```



## Reference
* [On Proving Pairings](https://eprint.iacr.org/2024/640)
* [Optimal Pairings](https://eprint.iacr.org/2008/096)
* [A remark on the computation of cube roots in finite fields](https://eprint.iacr.org/2009/457)