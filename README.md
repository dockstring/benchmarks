# Dockstring benchmarks

Official repository for dockstring-based benchmarks
(currently the regression, virtual screening, and optimization benchmarks from our paper).
This repository contains example code to help understand what the benchmarks
entail, and to help you to test new methods on the benchmarks.

## Original benchmarks

To showcase the benchmarks from the
[original `dockstring` paper](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01334),
we have created a set of jupyter notebooks which tests
a dummy model on all the benchmarks.
If you want to test a model on dockstring, we recommend that you
modify these notebooks and replace our dummy model with your model.

The following python packages are required to run the benchmarks:
- numpy
- pandas
- scikit-learn
- jupyter (to run the notebook)
- rdkit

The notebooks are linked here:
- [Regression benchmark](TODO)
- Virtual Screening benchmark: unfortunately we are still awaiting permission to re-distribute the version of the ZINC dataset which we used to run this benchmark. Stay tuned or email us for more information!
- [Optimization benchmark](TODO)

## New benchmarks?

Have you created a new benchmark that uses `dockstring`?
We would be happy to put it into this repository.
Feel free to submit a pull request to this repo with your benchmark
and we can add it to this repo!
