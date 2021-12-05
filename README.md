### Status

![Build Status](https://github.com/jacobwilliams/carlson-elliptic-integrals/actions/workflows/CI.yml/badge.svg)

### Description

This is a modern Fortran implementation of the Carlson symmetric forms of elliptic integrals code from the [SLATEC library](http://www.netlib.org/slatec/src/). It has been extensively refactored.

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `carson-elliptic-intergrals` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
carson-elliptic-intergrals = { git="https://github.com/jacobwilliams/carson-elliptic-intergrals.git" }
```

### References

1. B. C. Carlson and E. M. Notis, [Algorithms for incomplete
   elliptic integrals](http://dl.acm.org/citation.cfm?id=355970),
   ACM Transactions on Mathematical
   Software 7, 3 (September 1981), pp. 398-403.
2. B. C. Carlson, [Computing elliptic integrals by
   duplication](http://link.springer.com/article/10.1007%2FBF01396491), Numerische Mathematik 33, (1979),
   pp. 1-16.
3. B. C. Carlson, [Elliptic integrals of the first kind](http://epubs.siam.org/doi/abs/10.1137/0508016),
   SIAM Journal of Mathematical Analysis 8, (1977),
   pp. 231-242.

### Documentation

The latest API documentation for the `master` branch can be found [here](http://jacobwilliams.github.io/carlson-elliptic-integrals/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### License

The carlson-elliptic-integrals source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/carlson-elliptic-integrals/blob/master/LICENSE) (BSD-style).  The original Fortran 77 code is [public domain](http://www.netlib.org/slatec/guide).
