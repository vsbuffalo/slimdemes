# slimdemes

This is a simple implementation of
[demes](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html#sec-intro)
demographic models for SLiM. This is heavily inspired
by Graham Gower's
[deme-slims](https://github.com/grahamgower/demes-slim)
which is a few years old and not compatible with SLiM
version 4. **Note that this implementation does not
support gene flow (e.g. the migrations and pulses
blocks in the deme specification).

Validation is done entirely through the
[demes](https://popsim-consortium.github.io/demes-docs/latest/introduction.html)
Python package. This is, in my view, a cleaner way to
interface SLiM with demes, since SLiM does not have
native YAML parsing and requires that demes YAML
files be convert to JSON anyways. The `slimdemes
convert` tool converts a YAML demes model to JSON,
removes the `migrations` and `pulses` blocks (since
these are not yet supported), and converts the demes
`time_units` to generations based on the generation
time.
