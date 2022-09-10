![JuDGE](docs/src/assets/judge-small.png)

| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-latest-img]][docs-latest-url] | [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]

JuDGE stands for: Julia Decomposition for Generalized Expansion. Functionally,
it is a solver which leverages the syntax of the JuMP modelling language to
solve complex multi-stage (and multi-horizon) capacity expansion problems.

JuDGE supports generalized set partitioning problems. This, however, is
experimental.

Please see the [documentation](https://EPOC-NZ.github.io/JuDGE.jl/)
for details about installing JuDGE.jl, and examples showing how to set up a
stochastic capacity example model using the JuDGE.jl package.

For more details see our working paper: [JuDGE.jl: a Julia package for optimizing capacity expansion](http://www.optimization-online.org/DB_HTML/2020/11/8086.html).

[build-img]: https://github.com/EPOC-NZ/JuDGE.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/EPOC-NZ/JuDGE.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/github/EPOC-NZ/JuDGE.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/EPOC-NZ/JuDGE.jl?branch=master

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://EPOC-NZ.github.io/JuDGE.jl
