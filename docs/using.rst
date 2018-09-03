Using Fundem
============


Function Signatures
-------------------

This library supports R, Python, and C++, so function signatures are different
from one language to the next. The functions are all pretty similar, so
they follow these rules.

 * Function names for R and Python use underscores between words,
   while C++ names are CamelCase.

 * `nx` is the set of age group interval widths, known as :math:`n_x`
   in demography, as an array whose length
   is the number of age groups. That's a `c()` array in R, a
   `numpy.ndarray` with `dtype=numpy.float` in Python, and
   a `double *` in C++.

   If the last age group is half-open, such as 95+, then its
   interval should be infinity, as a floating point value.

 * `mx`, `ax`, `px`, and friends are assumed to be two-dimensional
   arrays, so they are treated as a long array of
   of the same length as `nx`, multiplied by some :math:`N` that's
   called the population count, or `pop_cnt`. Therefore,
   `array[pop,age]` is a single string of length `pop * age`.

   In R, every `c()` array looks like a one-dimensional array of
   the correct type as long as the ages are the last dimension.

   In Python, every `ndarray` is the right shape as long as
   ages are the last dimension. An XArray's values are a suitable
   choice.

   In C++, the input array of type `double *` is assumed to be of
   length `age_cnt * pop_cnt`.

 * `age_cnt` is an integer. This is an `int` for R and Python.

 * `pop_cnt` is a size_t in C++, but this, too, looks like an
   `int` to R and Python.

 * Return values will be returned from the function for R and
   Python, but they are passed into C++ in order to simplify
   the memory allocation.



.. index:: first_moment, uniform_deaths, balducci, constant_mortality

Naming
------

There are sets of functions that go together, and these have
prefixed or suffixed names.

 * `first_moment` refers to using :math:`({}_nm_x,{}_na_x)` together
   in order to specify a life table. It is an exact specification.

 * `uniform_deaths` refers to an assumption that the continuous distribution of
   deaths over an interval is uniform. Therefore, :math:`{}_na_x=n_x/2`.

 * `balducci` uses the Balducci Hypothesis. Actually, nobody uses that
   any more, but it's a great name.

 * `constant_mortality` applies the assumption of a constant continuous
   mortality rate to derive the mean age of death.
